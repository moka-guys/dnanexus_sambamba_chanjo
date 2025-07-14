#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os
import json
import csv
from pathlib import Path
from typing import List, Dict

# --- Configuration ---
SAMBAMBA_OUTPUT_FILENAME = "sambamba_output.bed"
CHANJO_DB_NAME = "coverage.sqlite3"
CHANJO_CONFIG_FILENAME = "chanjo.yaml"
CHANJO_JSON_OUTPUT_FILENAME = "chanjo_raw_output.json"
GENE_SUMMARY_FILENAME = "gene_summary.txt"
EXON_REPORT_FILENAME = "exon_level_report.csv"
GENE_REPORT_FILENAME = "gene_level_report.csv"


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
        f"  cov_treshold:\n"
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
    # Chanjo expects a specific 7-column BED format, but our file has 8 columns
    # Convert: chr, start, end, name, transcript, gene_id, gene_symbol
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
            if len(parts) >= 8:
                gene_id = parts[7]
                # Skip entries with non-numeric gene IDs (like Ensembl IDs)
                if not gene_id.isdigit():
                    skipped_count += 1
                    continue
                    
                # Split gene;transcript field
                gene_parts = parts[6].split(';')
                gene_symbol = gene_parts[0] if gene_parts else 'UNKNOWN'
                transcript = gene_parts[1] if len(gene_parts) > 1 else 'UNKNOWN'
                # Output in chanjo format: chr, start, end, name, transcript, gene_id, gene_symbol
                f_out.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{parts[3]}\t{transcript}\t{gene_id}\t{gene_symbol}\n")
                processed_count += 1
    
    print(f"   - Converted BED file to chanjo format: {chanjo_bed_file}")
    print(f"   - Processed {processed_count} entries, skipped {skipped_count} entries with non-numeric gene IDs")
    
    run_command(["chanjo", "link", str(chanjo_bed_file)], cwd=output_dir)
    run_command(["chanjo", "load", str(sambamba_output)], cwd=output_dir)
    print(f"   - Linked transcript definitions from: {chanjo_bed_file}")
    print(f"   - Loaded coverage data from: {sambamba_output}")

    # 4. Get unique gene list from original BED file (column 7 - EntrezID)
    gene_ids = set()
    with open(args.bed_file, "r") as f_in:
        for line in f_in:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 8:
                gene_id = parts[7]
                if gene_id and gene_id != '0':
                    gene_ids.add(gene_id)
    
    unique_genes = sorted(list(gene_ids))
    print(f"   - Found {len(unique_genes)} unique gene IDs to analyze.")

    # 5. Get overall statistics from chanjo
    result = run_command(["chanjo", "calculate", "mean"], cwd=output_dir)
    with open(chanjo_json_output, "w") as f_out:
        f_out.write(result.stdout)

    print(f"✅ Chanjo analysis complete. Overall statistics saved to: {chanjo_json_output}")
    return chanjo_json_output


def process_results(sambamba_bed: Path, chanjo_json: Path, args: argparse.Namespace, output_dir: Path):
    """
    Parses the outputs from sambamba and chanjo to create final reports.
    This function is a Python 3 rewrite of the logic in `read_chanjo.py`.
    """
    print("\n--- Processing and Generating Final Reports ---")
    
    # --- Part 1: Create a simple summary from Chanjo's JSON output (overall stats) ---
    gene_summary_path = output_dir / GENE_SUMMARY_FILENAME
    with open(chanjo_json, 'r') as f_in, open(gene_summary_path, 'w') as f_out:
        data = json.load(f_in)
        sample_id = data.get('sample_id', 'unknown')
        mean_coverage = data.get('mean_coverage', 'N/A')
        f_out.write(f"Sample ID: {sample_id}\n")
        f_out.write(f"Mean Coverage: {mean_coverage}\n")
        for key, value in data.items():
            if key.startswith('completeness_'):
                level = key.split('_')[1]
                f_out.write(f"Completeness at {level}X: {value}\n")
    print(f"   - Generated overall summary: {gene_summary_path}")

    # --- Part 2: Create an exon-level report for regions below 100% coverage ---
    exon_report_path = output_dir / EXON_REPORT_FILENAME
    with open(sambamba_bed, 'r') as f_in, open(exon_report_path, 'w', newline='') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["gene_symbol", "transcript", "entrez_id", "chr", "start", "stop", "mean_coverage", "completeness"])
        
        for line in f_in:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            # Expected format: chrom, start, end, name, score, strand, gene;transcript, entrezID, ...
            completeness = float(fields[10])

            if completeness < 100.0:
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
                
                writer.writerow([gene_symbol, transcript, entrez_id, chrom, start, stop, mean_coverage, completeness])
                
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

    print(f"   - Generated gene-level report: {gene_report_path}")
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
    sambamba_output = run_sambamba(args, args.output_dir)
    chanjo_output = run_chanjo(sambamba_output, args, args.output_dir)
    process_results(sambamba_output, chanjo_output, args, args.output_dir)


if __name__ == "__main__":
    main()