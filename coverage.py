#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os
import json
import csv
import re
import urllib.request
import tempfile
from pathlib import Path
from typing import List, Dict, Set, Optional

# --- Configuration ---
SAMBAMBA_OUTPUT_FILENAME = "sambamba_output.bed"
CHANJO_DB_NAME = "coverage.sqlite3"
CHANJO_CONFIG_FILENAME = "chanjo.yaml"
CHANJO_JSON_OUTPUT_FILENAME = "chanjo_raw_output.json"
GENE_SUMMARY_FILENAME = "gene_summary.txt"
EXON_REPORT_FILENAME = "exon_level_report.csv"
GENE_REPORT_FILENAME = "gene_level_report.csv"

# --- Panel Configuration URLs ---
PANEL_CONFIG_URL = "https://raw.githubusercontent.com/moka-guys/automate_demultiplex/main/config/panel_config.py"
DNANEXUS_PROJECT = "project-ByfFPz00jy1fk6PjpZ95F27J"


def extract_pan_number_from_bam(bam_filename: str) -> Optional[str]:
    """
    Extracts the Pan number from a BAM filename.
    
    Args:
        bam_filename: BAM filename like 'NGS698PoolA_15_368301_JB_F_CP2R207Via_Pan4150.bam'
    
    Returns:
        Pan number like 'Pan4150' or None if not found
    """
    # Look for Pan followed by digits
    match = re.search(r'(Pan\d+)', bam_filename)
    return match.group(1) if match else None


def download_panel_config() -> str:
    """
    Downloads the panel configuration file from the automate_demultiplex repository.
    
    Returns:
        Content of the panel configuration file
    """
    print(f"   - Downloading panel configuration from: {PANEL_CONFIG_URL}")
    with urllib.request.urlopen(PANEL_CONFIG_URL) as response:
        return response.read().decode('utf-8')


def extract_ed_cnvcalling_bedfile(panel_config_content: str, pan_number: str) -> Optional[str]:
    """
    Extracts the ed_cnvcalling_bedfile value for a specific Pan number.
    
    Args:
        panel_config_content: Content of the panel_config.py file
        pan_number: Pan number like 'Pan4150'
    
    Returns:
        ed_cnvcalling_bedfile value or None if not found
    """
    # Look for the Pan configuration block
    pan_pattern = rf'"{pan_number}":\s*\{{'
    match = re.search(pan_pattern, panel_config_content)
    
    if not match:
        return None
    
    # Find the ed_cnvcalling_bedfile line within this Pan block
    # Look from the Pan definition to the next closing brace
    start_pos = match.end()
    
    # Find the matching closing brace (simple approach - look for next },)
    brace_count = 1
    pos = start_pos
    while pos < len(panel_config_content) and brace_count > 0:
        if panel_config_content[pos] == '{':
            brace_count += 1
        elif panel_config_content[pos] == '}':
            brace_count -= 1
        pos += 1
    
    pan_block = panel_config_content[start_pos:pos-1]
    
    # Look for ed_cnvcalling_bedfile in this block
    bedfile_match = re.search(r'"ed_cnvcalling_bedfile":\s*"([^"]+)"', pan_block)
    return bedfile_match.group(1) if bedfile_match else None


def download_panel_bed_file(bedfile_name: str, output_dir: Path) -> Optional[Path]:
    """
    Downloads the panel-specific BED file from DNAnexus.
    
    Args:
        bedfile_name: Name like 'Pan5250'
        output_dir: Directory to save the file
    
    Returns:
        Path to downloaded BED file or None if failed
    """
    bed_filename = f"{bedfile_name}_CNV.bed"
    dx_path = f"{DNANEXUS_PROJECT}:/Data/BED/{bed_filename}"
    local_path = output_dir / bed_filename
    
    # Check if file already exists
    if local_path.exists():
        print(f"   - Panel BED file already exists: {local_path}")
        return local_path
    
    print(f"   - Downloading panel BED file: {dx_path}")
    
    try:
        # Use dx download command
        cmd = ["dx", "download", dx_path, "-o", str(local_path)]
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"   - Downloaded panel BED file to: {local_path}")
        return local_path
    except subprocess.CalledProcessError as e:
        print(f"   - Warning: Failed to download panel BED file: {e}")
        print(f"     Error: {e.stderr}")
        return None
    except FileNotFoundError:
        print(f"   - Warning: dx command not found. Cannot download panel BED file.")
        print(f"     Please ensure DNAnexus CLI is installed and configured.")
        return None


def extract_panel_genes(panel_bed_file: Path) -> Set[str]:
    """
    Extracts unique gene symbols from the panel BED file.
    
    Args:
        panel_bed_file: Path to the panel BED file
    
    Returns:
        Set of gene symbols found in the BED file
    """
    genes = set()
    
    with open(panel_bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                # Extract gene symbol from the 4th column (name field)
                # Format: MSH2_EPCAM_1;ENSG00000095002 or MSH2;NM_000251.2 or BRCA1_3UTR;transcript
                name_field = parts[3]
                # Take the part before the semicolon
                gene_region = name_field.split(';')[0]
                # Extract only the gene symbol before the first underscore
                # This removes user-defined region suffixes like _3UTR, _IN11_1, _PM_5_1, etc.
                gene_symbol = gene_region.split('_')[0]
                if gene_symbol:  # Only add non-empty gene symbols
                    genes.add(gene_symbol)
    
    return genes


def get_panel_genes(bam_file: Path, output_dir: Path) -> Optional[Set[str]]:
    """
    Gets the set of genes for the panel based on the BAM filename.
    
    Args:
        bam_file: Path to the BAM file
        output_dir: Directory for temporary files
    
    Returns:
        Set of gene symbols for this panel, or None if panel filtering unavailable
    """
    print("\n--- Panel Gene Filtering ---")
    
    # Extract Pan number from BAM filename
    pan_number = extract_pan_number_from_bam(bam_file.name)
    if not pan_number:
        print(f"   - No Pan number found in BAM filename: {bam_file.name}")
        return None
    
    print(f"   - Detected Pan number: {pan_number}")
    
    try:
        # Download panel configuration
        panel_config = download_panel_config()
        
        # Extract ed_cnvcalling_bedfile for this Pan
        bedfile_name = extract_ed_cnvcalling_bedfile(panel_config, pan_number)
        if not bedfile_name:
            print(f"   - No ed_cnvcalling_bedfile found for {pan_number}")
            return None
        
        print(f"   - Panel BED file: {bedfile_name}")
        
        # Download the panel BED file
        panel_bed_file = download_panel_bed_file(bedfile_name, output_dir)
        if not panel_bed_file:
            return None
        
        # Extract genes from the panel BED file
        panel_genes = extract_panel_genes(panel_bed_file)
        print(f"   - Found {len(panel_genes)} unique genes in panel: {', '.join(sorted(list(panel_genes))[:10])}{'...' if len(panel_genes) > 10 else ''}")
        
        return panel_genes
        
    except Exception as e:
        print(f"   - Error during panel gene extraction: {e}")
        return None


def run_command(command: List[str], cwd: Path = None):
    """Executes a command and raises an error if it fails."""
    print(f"🚀 Running command: {' '.join(command)}")
    try:
        process = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            cwd=cwd,
        )
        return process
    except FileNotFoundError as e:
        print(f"❌ Error: Command '{command[0]}' not found.")
        print("Please ensure sambamba and chanjo are installed and in your system's PATH.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"❌ Command failed: {' '.join(command)}")
        print(f"   Return Code: {e.returncode}")
        print(f"   Stdout: {e.stdout}")
        print(f"   Stderr: {e.stderr}")
        sys.exit(1)


def run_sambamba(args: argparse.Namespace, output_dir: Path) -> Path:
    """
    Constructs and executes the sambamba depth region command.
    """
    print("\n--- Running Sambamba ---")
    output_bed = output_dir / SAMBAMBA_OUTPUT_FILENAME

    # Build the filter command for the -F flag
    filter_conditions = [f"mapping_quality >= {args.min_mapping_qual}"]
    if args.exclude_failed_qc:
        filter_conditions.append("not failed_quality_control")
    if args.exclude_duplicates:
        filter_conditions.append("not duplicate")
    if args.additional_filter:
        filter_conditions.append(args.additional_filter)
    
    filter_command_str = " and ".join(filter_conditions)

    # Build the main sambamba command
    command = [
        "sambamba", "depth", "region",
        "-L", str(args.bed_file),
        "-t", str(args.threads),
        "-T", str(args.coverage_level),
        "--min-base-quality", str(args.min_base_qual),
        "-F", filter_command_str,
    ]
    
    if args.merge_overlapping_mates:
        command.append("-m")
        
    if args.additional_sambamba_flags:
        command.extend(args.additional_sambamba_flags.split())

    command.append(str(args.bam_file))

    # Execute and redirect output to a file
    result = run_command(command)
    with open(output_bed, "w") as f_out:
        f_out.write(result.stdout)
        
    print(f"✅ Sambamba analysis complete. Output saved to: {output_bed}")
    return output_bed


def run_chanjo(sambamba_output: Path, args: argparse.Namespace, output_dir: Path) -> Path:
    """
    Initializes a chanjo database, loads original BED data, and calculates coverage.
    """
    print("\n--- Running Chanjo ---")
    chanjo_db_path = output_dir / CHANJO_DB_NAME
    chanjo_config_path = output_dir / CHANJO_CONFIG_FILENAME
    chanjo_json_output = output_dir / CHANJO_JSON_OUTPUT_FILENAME
    
    # 1. Create chanjo config file
    config_content = (
        f"database: {chanjo_db_path.name}\n"
        f"sambamba:\n"
        f"  cov_threshold:\n"
        f"    - {args.coverage_level}"
    )
    with open(chanjo_config_path, "w") as f:
        f.write(config_content)
    print(f"   - Wrote chanjo config to: {chanjo_config_path}")

    # 2. Initialize and set up the database
    if chanjo_db_path.exists():
        print(f"   - Deleting existing database: {chanjo_db_path}")
        chanjo_db_path.unlink()
        
    run_command(["chanjo", "init", "-a"], cwd=output_dir)
    run_command(["chanjo", "db", "setup"], cwd=output_dir)
    print(f"   - Initialized new database: {chanjo_db_path}")

    # 3. Convert BED file to chanjo format and link/load
    # Chanjo expects a specific 7-column BED format, but our file has 7 columns in a different arrangement
    # Input format: chr, start, end, name, transcript, gene_id, gene_symbol  
    # Chanjo format: chr, start, end, name, transcript, gene_id, gene_symbol
    chanjo_bed_file = output_dir / "chanjo_format.bed"
    bed_file_abs = args.bed_file.resolve()
    
    # Convert BED file to chanjo format
    skipped_count = 0
    processed_count = 0
    with open(bed_file_abs, 'r') as f_in, open(chanjo_bed_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 7:  # 7-column format
                gene_id = parts[5]  # Column 6
                gene_symbol = parts[6]  # Column 7
                
                # Skip entries with non-numeric gene IDs (like Ensembl IDs)
                if not gene_id.isdigit():
                    skipped_count += 1
                    continue
                    
                # Output in chanjo format: chr, start, end, name, transcript, gene_id, gene_symbol
                # Use parts[4] as transcript since the input format seems to have transcript in column 5
                transcript = parts[4] if len(parts) > 4 else 'UNKNOWN'
                f_out.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\t{transcript}\t{gene_id}\t{gene_symbol}\n")
                processed_count += 1
    
    print(f"   - Converted BED file to chanjo format: {chanjo_bed_file}")
    print(f"   - Processed {processed_count} entries, skipped {skipped_count} entries with non-numeric gene IDs")
    
    run_command(["chanjo", "link", str(chanjo_bed_file)], cwd=output_dir)
    run_command(["chanjo", "load", str(sambamba_output)], cwd=output_dir)
    print(f"   - Linked transcript definitions from: {chanjo_bed_file}")
    print(f"   - Loaded coverage data from: {sambamba_output}")

    # 4. Get unique gene list from original BED file (column 6 - EntrezID) 
    gene_ids = set()
    with open(args.bed_file, "r") as f_in:
        for line in f_in:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 6:  # Changed from 8 to 6 for 7-column format
                gene_id = parts[5]  # Changed from parts[7] to parts[5] (column 6)
                if gene_id and gene_id != '0' and gene_id.isdigit():
                    gene_ids.add(gene_id)
    
    unique_genes = sorted(list(gene_ids))
    print(f"   - Found {len(unique_genes)} unique gene IDs to analyze.")

    # 5. Generate per-gene coverage statistics in legacy format
    # Instead of using chanjo calculate coverage, we'll process sambamba output directly
    # to match the legacy format: {"genes": {"gene_id": {"completeness_30": X, "mean_coverage": Y}}}
    
    print(f"   - Generating per-gene coverage statistics for {len(unique_genes)} genes...")
    
    # Calculate per-gene statistics from sambamba output
    gene_stats = {}
    sample_id = os.path.splitext(os.path.basename(sambamba_output))[0]
    
    with open(sambamba_output, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 11:
                # Extract gene info and coverage data
                gene_info = fields[6] if len(fields) > 6 else ""
                entrez_id = fields[7] if len(fields) > 7 else ""
                mean_coverage = float(fields[9]) if len(fields) > 9 and fields[9] != '.' else 0.0
                completeness = float(fields[10]) if len(fields) > 10 and fields[10] != '.' else 0.0
                
                # Only include numeric gene IDs
                if entrez_id and entrez_id.isdigit():
                    if entrez_id not in gene_stats:
                        gene_stats[entrez_id] = {'coverage_values': [], 'completeness_values': []}
                    
                    gene_stats[entrez_id]['coverage_values'].append(mean_coverage)
                    gene_stats[entrez_id]['completeness_values'].append(completeness)
    
    # Generate JSON output in legacy format - one line per gene
    with open(chanjo_json_output, "w") as f_out:
        for gene_id in sorted(gene_stats.keys(), key=int):
            stats = gene_stats[gene_id]
            avg_coverage = sum(stats['coverage_values']) / len(stats['coverage_values']) if stats['coverage_values'] else 0.0
            avg_completeness = sum(stats['completeness_values']) / len(stats['completeness_values']) if stats['completeness_values'] else 0.0
            
            # Create legacy format JSON object
            gene_json = {
                "genes": {
                    gene_id: {
                        f"completeness_{args.coverage_level}": avg_completeness,
                        "mean_coverage": avg_coverage
                    }
                },
                "sample_id": sample_id
            }
            
            f_out.write(json.dumps(gene_json) + "\n")

    print(f"✅ Chanjo analysis complete. Per-gene statistics saved to: {chanjo_json_output}")
    return chanjo_json_output


def process_results(sambamba_bed: Path, chanjo_json: Path, args: argparse.Namespace, output_dir: Path, panel_genes: Optional[Set[str]] = None):
    """
    Parses the outputs from sambamba and chanjo to create final reports.
    This function is a Python 3 rewrite of the logic in `read_chanjo.py`.
    
    Args:
        sambamba_bed: Path to sambamba output BED file
        chanjo_json: Path to chanjo JSON output file 
        args: Command line arguments
        output_dir: Output directory
        panel_genes: Optional set of gene symbols to filter results (panel-specific filtering)
    """
    print("\n--- Processing and Generating Final Reports ---")
    
    if panel_genes:
        print(f"   - Applying panel-specific filtering for {len(panel_genes)} genes")
    
    # --- Part 1: Create a simple summary from Chanjo's JSON output (overall stats) ---
    gene_summary_path = output_dir / GENE_SUMMARY_FILENAME
    
    # Read the line-delimited JSON and collect all genes
    all_genes = {}
    sample_id = "unknown"
    
    with open(chanjo_json, 'r') as f_in:
        for line in f_in:
            if line.strip():
                data = json.loads(line)
                sample_id = data.get('sample_id', sample_id)
                genes_data = data.get('genes', {})
                all_genes.update(genes_data)
    
    # Calculate overall statistics
    if all_genes:
        total_coverage = sum(gene_data.get('mean_coverage', 0) for gene_data in all_genes.values())
        avg_coverage = total_coverage / len(all_genes)
        
        # Get coverage level from first gene's completeness key
        coverage_level = args.coverage_level
        completeness_key = f"completeness_{coverage_level}"
        total_completeness = sum(gene_data.get(completeness_key, 0) for gene_data in all_genes.values())
        avg_completeness = total_completeness / len(all_genes)
    else:
        avg_coverage = 0
        avg_completeness = 0
    
    with open(gene_summary_path, 'w') as f_out:
        f_out.write(f"Sample ID: {sample_id}\n")
        f_out.write(f"Mean Coverage: {avg_coverage:.2f}\n")
        f_out.write(f"Completeness at {args.coverage_level}X: {avg_completeness:.2f}\n")
        f_out.write(f"Total Genes Analyzed: {len(all_genes)}\n")
        if panel_genes:
            f_out.write(f"Panel Genes: {len(panel_genes)} genes\n")
    
    print(f"   - Generated overall summary: {gene_summary_path}")

    # --- Part 2: Create an exon-level report for regions below 100% coverage ---
    exon_report_path = output_dir / EXON_REPORT_FILENAME
    exon_count = 0
    filtered_exon_count = 0
    
    with open(sambamba_bed, 'r') as f_in, open(exon_report_path, 'w', newline='') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["gene_symbol", "transcript", "entrez_id", "chr", "start", "stop", "mean_coverage", "completeness"])
        
        for line in f_in:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            # Expected format: chrom, start, end, name, score, strand, gene;transcript, entrezID, ...
            if len(fields) < 11:
                continue
                
            completeness = float(fields[10])

            if completeness < 100.0:
                exon_count += 1
                gene_info = fields[6]
                entrez_id = fields[7]
                coords = fields[3] # Expects format like "chr1-12345-67890"
                mean_coverage = fields[9]

                # Parse coordinates
                try:
                    chrom, start, stop = coords.split('-')
                except ValueError:
                    chrom, start, stop = "N/A", "N/A", "N/A" # Handle unexpected format

                # Parse gene symbol and transcript
                if ";" in gene_info:
                    gene_symbol, transcript = gene_info.split(";", 1)
                else:
                    gene_symbol, transcript = gene_info, "N/A"
                
                # Apply panel filtering if enabled
                if panel_genes and gene_symbol not in panel_genes:
                    filtered_exon_count += 1
                    continue
                
                writer.writerow([gene_symbol, transcript, entrez_id, chrom, start, stop, mean_coverage, completeness])
                
    if panel_genes and filtered_exon_count > 0:
        print(f"   - Generated exon-level report: {exon_report_path} (filtered {filtered_exon_count}/{exon_count} exons not in panel)")
    else:
        print(f"   - Generated exon-level report: {exon_report_path}")

    # --- Part 3: Create a gene-level report from sambamba data ---
    # Calculate average coverage and completeness per gene from sambamba output
    gene_stats: Dict[str, Dict[str, List[float]]] = {}
    
    with open(sambamba_bed, 'r') as f_in:
        for line in f_in:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 11:  # sambamba output has more columns
                gene_info = fields[6]
                entrez_id = fields[7]
                mean_coverage = float(fields[9]) if fields[9] != '.' else 0.0
                completeness = float(fields[10]) if fields[10] != '.' else 0.0
                
                # Skip non-numeric gene IDs
                if not entrez_id.isdigit():
                    continue
                
                # Extract gene symbol
                gene_symbol = gene_info.split(';')[0] if ';' in gene_info else gene_info
                
                # Special hardcoded case from original script
                if entrez_id == "11200":
                    gene_symbol = "CHEK2"
                
                # Apply panel filtering if enabled
                if panel_genes and gene_symbol not in panel_genes:
                    continue
                
                if gene_symbol not in gene_stats:
                    gene_stats[gene_symbol] = {'coverage': [], 'completeness': []}
                
                gene_stats[gene_symbol]['coverage'].append(mean_coverage)
                gene_stats[gene_symbol]['completeness'].append(completeness)

    # Generate gene-level report
    gene_report_path = output_dir / GENE_REPORT_FILENAME
    with open(gene_report_path, 'w', newline='') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["gene_symbol", f"mean_coverage", f"average_completeness_at_{args.coverage_level}X", "exon_count"])
        
        for gene_symbol, stats in gene_stats.items():
            avg_coverage = sum(stats['coverage']) / len(stats['coverage']) if stats['coverage'] else 0
            avg_completeness = sum(stats['completeness']) / len(stats['completeness']) if stats['completeness'] else 0
            exon_count = len(stats['coverage'])
            writer.writerow([gene_symbol, f"{avg_coverage:.2f}", f"{avg_completeness:.2f}", exon_count])

    total_genes = len(gene_stats)
    if panel_genes:
        print(f"   - Generated gene-level report: {gene_report_path} ({total_genes} panel genes)")
    else:
        print(f"   - Generated gene-level report: {gene_report_path} ({total_genes} genes)")
        
    print(f"\n✅ All reports generated successfully in: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="A Python 3 script to run a sambamba and chanjo coverage analysis workflow.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # --- Input Files ---
    parser.add_argument("--bam-file", type=Path, required=True, help="Path to the input BAM file.")
    parser.add_argument("--bed-file", type=Path, required=True, help="Path to the input BED file.")

    # --- Sambamba Settings ---
    parser.add_argument("-T", "--coverage-level", type=int, default=30, help="Sambamba: Minimum coverage level required.")
    parser.add_argument("--min-mapping-qual", type=int, default=20, help="Sambamba: Minimum mapping quality.")
    parser.add_argument("--min-base-qual", type=int, default=20, help="Sambamba: Minimum base quality.")
    parser.add_argument("--exclude-failed-qc", action='store_true', help="Sambamba: Exclude reads that failed quality control.")
    parser.add_argument("--exclude-duplicates", action='store_true', help="Sambamba: Exclude duplicate reads.")
    parser.add_argument("--merge-overlapping-mates", action='store_true', help="Sambamba: Count overlapping mate reads only once (-m flag).")
    parser.add_argument("--additional-filter", type=str, help="Sambamba: Additional filter commands (e.g., 'first_of_pair').")
    parser.add_argument("--additional-sambamba-flags", type=str, help="Sambamba: Other flags to pass, enclosed in quotes.")

    # --- General Settings ---
    parser.add_argument("-t", "--threads", type=int, default=os.cpu_count(), help="Number of threads to use.")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path.cwd() / "coverage_results", help="Directory to store all outputs.")

    # --- Panel Filtering ---
    parser.add_argument("--panel-filter", action='store_true', help="Enable panel-specific gene filtering based on BAM filename.")

    args = parser.parse_args()

    # --- Basic Input Validation ---
    if not args.bam_file.is_file():
        print(f"❌ Error: BAM file not found at '{args.bam_file}'")
        sys.exit(1)
        
    bai_file = args.bam_file.with_suffix(".bam.bai")
    if not bai_file.is_file():
        print(f"❌ Error: BAM index file (.bam.bai) not found for '{args.bam_file}'. Expected at '{bai_file}'")
        sys.exit(1)
        
    if not args.bed_file.is_file():
        print(f"❌ Error: BED file not found at '{args.bed_file}'")
        sys.exit(1)

    # --- Create Output Directory ---
    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"📂 All outputs will be saved in: {args.output_dir.resolve()}")

    # --- Execute Workflow ---
    if args.panel_filter:
        panel_genes = get_panel_genes(args.bam_file, args.output_dir)
    else:
        panel_genes = None

    sambamba_output = run_sambamba(args, args.output_dir)
    chanjo_output = run_chanjo(sambamba_output, args, args.output_dir)
    process_results(sambamba_output, chanjo_output, args, args.output_dir, panel_genes)


if __name__ == "__main__":
    main()