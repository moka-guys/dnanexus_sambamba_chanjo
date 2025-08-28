#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os
from pathlib import Path
from typing import List, Dict

# --- Configuration ---
SAMBAMBA_OUTPUT_FILENAME = "sambamba_output.bed"
EXON_REPORT_FILENAME = "exon_level.txt"
GENE_REPORT_FILENAME = "gene_level.txt"


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

    # Execute sambamba with suppressed stderr to avoid verbose warnings
    try:
        stderr_setting = None if args.verbose else subprocess.DEVNULL
        
        process = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=stderr_setting,
            text=True,
        )
        
        with open(output_bed, "w") as f_out:
            f_out.write(process.stdout)
            
    except FileNotFoundError as e:
        print(f"Error: Command 'sambamba' not found.")
        print("Please ensure sambamba is installed and in your system's PATH.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Sambamba command failed with return code: {e.returncode}")
        if not args.verbose:
            print("Tip: Use --verbose flag to see sambamba's detailed error messages")
        print("This may indicate an issue with the input files or command parameters.")
        sys.exit(1)
        
    print(f"Sambamba analysis complete. Output saved to: {output_bed}")
    return output_bed


def process_results(sambamba_bed: Path, args: argparse.Namespace, output_dir: Path):
    """
    Parses the outputs from sambamba to create final reports.
    This function is a Python 3 rewrite of the logic in `read_chanjo.py`.
    
    Args:
        sambamba_bed: Path to sambamba output BED file
        args: Command line arguments
        output_dir: Output directory
    """
    print("\n--- Processing and Generating Final Reports ---")
    
    # --- Part 1: Create gene-level statistics from sambamba output ---
    gene_stats: Dict[str, Dict[str, List[float]]] = {}
    
    with open(sambamba_bed, 'r') as f_in:
        for line in f_in:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 11:
                # Extract gene info and coverage data
                gene_info = fields[6] # Column 7: gene;transcript format
                entrez_id = fields[7] # Column 8: Entrez ID 
                coords = fields[3] # Column 4: chr-start-stop format
                mean_coverage = float(fields[9]) if len(fields) > 9 and fields[9] != '.' else 0.0 # Column 10
                completeness = float(fields[10]) if len(fields) > 10 and fields[10] != '.' else 0.0 # Column 11
                
                # Skip non-numeric gene IDs
                if not entrez_id.isdigit():
                    continue
                
                # Extract gene symbol
                gene_symbol = gene_info.split(';')[0] if ';' in gene_info else gene_info
                
                # Special hardcoded case from original script
                if entrez_id == "11200":
                    gene_symbol = "CHEK2"
                
                # Calculate exon length from coordinates
                try:
                    chrom, start, stop = coords.split('-')
                    exon_length = int(stop) - int(start) + 1
                except (ValueError, IndexError):
                    print(f"   - Warning: Could not parse coordinates '{coords}' for {gene_symbol}, skipping")
                    continue
                
                # For gene-level statistics, only include coding regions (exclude UTR regions)
                # Extract base gene symbol (before any underscore) for gene-level aggregation
                base_gene_symbol = gene_symbol.split('_')[0]
                
                if base_gene_symbol not in gene_stats:
                    gene_stats[base_gene_symbol] = {
                        'coverage_values': [], 
                        'completeness_values': [], 
                        'exon_lengths': [],
                        'total_length': 0
                    }
                
                gene_stats[base_gene_symbol]['coverage_values'].append(mean_coverage)
                gene_stats[base_gene_symbol]['completeness_values'].append(completeness)
                gene_stats[base_gene_symbol]['exon_lengths'].append(exon_length)
                gene_stats[base_gene_symbol]['total_length'] += exon_length

    # Generate gene-level report with length-weighted averages
    gene_report_path = output_dir / GENE_REPORT_FILENAME
    with open(gene_report_path, 'w') as f_out:
        f_out.write("gene_symbol\tmean_coverage\taverage_completeness_at_{}X\texon_count\ttotal_length_bp\n".format(args.coverage_level))
        
        for gene_symbol in sorted(gene_stats.keys()):
            stats = gene_stats[gene_symbol]
            
            # Calculate length-weighted averages
            if stats['coverage_values'] and stats['total_length'] > 0:
                # Length-weighted mean coverage
                weighted_coverage_sum = sum(cov * length for cov, length in zip(stats['coverage_values'], stats['exon_lengths']))
                avg_coverage = weighted_coverage_sum / stats['total_length']
                
                # Length-weighted mean completeness
                weighted_completeness_sum = sum(comp * length for comp, length in zip(stats['completeness_values'], stats['exon_lengths']))
                avg_completeness = weighted_completeness_sum / stats['total_length']
            else:
                avg_coverage = 0
                avg_completeness = 0
            
            exon_count = len(stats['coverage_values'])
            total_length = stats['total_length']
            
            f_out.write(f"{gene_symbol}\t{avg_coverage:.2f}\t{avg_completeness:.2f}\t{exon_count}\t{total_length}\n")

    total_genes = len(gene_stats)
    print(f"   - Generated gene-level report: {gene_report_path} ({total_genes} genes)")
        
    # --- Part 2: Create an exon-level report for regions below 100% coverage ---
    exon_report_path = output_dir / EXON_REPORT_FILENAME
    coding_exon_data = []
    utr_exon_data = []
    total_exons_processed = 0
    
    # First pass: collect all exon data and count total exons
    with open(sambamba_bed, 'r') as f_in:
        for line in f_in:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
                
            total_exons_processed += 1
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
                
                # Store exon data for sorting, separating UTR from coding regions
                exon_record = {
                    'gene_symbol': gene_symbol,
                    'transcript': transcript,
                    'entrez_id': entrez_id,
                    'chrom': chrom,
                    'start': start,
                    'stop': stop,
                    'mean_coverage': mean_coverage,
                    'completeness': completeness
                }
                
                # Check if this is a UTR region (contains "_#UTR" pattern)
                if "_" in gene_symbol and "UTR" in gene_symbol:
                    utr_exon_data.append(exon_record)
                else:
                    coding_exon_data.append(exon_record)
    
    # Sort both datasets alphabetically by gene symbol, then by chromosome, then by start position
    coding_exon_data.sort(key=lambda x: (x['gene_symbol'], x['chrom'], int(x['start']) if x['start'].isdigit() else 0))
    utr_exon_data.sort(key=lambda x: (x['gene_symbol'], x['chrom'], int(x['start']) if x['start'].isdigit() else 0))
    
    # Write sorted data to file with separate sections
    coding_exon_count = len(coding_exon_data)
    utr_exon_count = len(utr_exon_data)
    total_under_coverage = coding_exon_count + utr_exon_count
    exons_at_100_percent = total_exons_processed - total_under_coverage
    
    with open(exon_report_path, 'w') as f_out:
        # Header
        f_out.write("gene_symbol\ttranscript\tentrez_id\tchr\tstart\tstop\tmean_coverage\tcompleteness\n")
        
        # Section 1: Coding regions with <100% coverage
        if coding_exon_count > 0:
            f_out.write(f"\n# CODING REGIONS WITH <100% COVERAGE AT {args.coverage_level}X ({coding_exon_count} exons)\n")
            for exon in coding_exon_data:
                f_out.write(f"{exon['gene_symbol']}\t{exon['transcript']}\t{exon['entrez_id']}\t{exon['chrom']}\t{exon['start']}\t{exon['stop']}\t{exon['mean_coverage']}\t{exon['completeness']}\n")
        
        # Section 2: UTR regions with <100% coverage
        if utr_exon_count > 0:
            f_out.write(f"\n# UTR REGIONS WITH <100% COVERAGE AT {args.coverage_level}X ({utr_exon_count} exons)\n")
            for exon in utr_exon_data:
                f_out.write(f"{exon['gene_symbol']}\t{exon['transcript']}\t{exon['entrez_id']}\t{exon['chrom']}\t{exon['start']}\t{exon['stop']}\t{exon['mean_coverage']}\t{exon['completeness']}\n")
        
        # Summary message at the end
        f_out.write(f"\n# Any exons not mentioned above are covered 100% at {args.coverage_level}X\n")

    print(f"   - Generated exon-level report: {exon_report_path}")
    print(f"     • {coding_exon_count} coding exons with <100% coverage")
    print(f"     • {utr_exon_count} UTR exons with <100% coverage") 
    print(f"     • {exons_at_100_percent} exons with 100% coverage at {args.coverage_level}X")

    print(f"\nAll reports generated successfully in: {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="A Python 3 script to run sambamba coverage analysis and generate detailed reports.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # --- Input Files ---
    parser.add_argument("--bam-file", type=Path, required=True, help="Path to the input BAM file.")
    parser.add_argument("--bed-file", type=Path, required=True, help="Path to the input BED file.")

    # --- Sambamba Settings ---
    parser.add_argument("-T", "--coverage-level", type=int, default=30, help="Sambamba: Minimum coverage level required.")
    parser.add_argument("--min-mapping-qual", type=int, default=20, help="Sambamba: Minimum mapping quality.")
    parser.add_argument("--min-base-qual", type=int, default=10, help="Sambamba: Minimum base quality.")
    parser.add_argument("--exclude-failed-qc", action='store_true', default=True, help="Sambamba: Exclude reads that failed quality control.")
    parser.add_argument("--exclude-duplicates", action='store_true', default=True, help="Sambamba: Exclude duplicate reads.")
    parser.add_argument("--merge-overlapping-mates", action='store_true', default=False, help="Sambamba: Count overlapping mate reads only once (-m flag).")
    parser.add_argument("--additional-filter", type=str, default="not (unmapped or secondary_alignment)", help="Sambaba: Additional filter commands (e.g., 'first_of_pair').")
    parser.add_argument("--additional-sambamba-flags", type=str, help="Sambamba: Other flags to pass, enclosed in quotes.")

    # --- General Settings ---
    parser.add_argument("-t", "--threads", type=int, default=os.cpu_count(), help="Number of threads to use.")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path.cwd() / "coverage_results", help="Directory to store all outputs.")

    # --- Verbosity ---
    parser.add_argument("--verbose", action='store_true', help="Show sambamba warnings and verbose output (default: suppressed).")

    args = parser.parse_args()

    # --- Basic Input Validation ---
    if not args.bam_file.is_file():
        print(f"Error: BAM file not found at '{args.bam_file}'")
        sys.exit(1)
        
    bai_file = args.bam_file.with_suffix(".bam.bai")
    if not bai_file.is_file():
        print(f"Error: BAM index file (.bam.bai) not found for '{args.bam_file}'. Expected at '{bai_file}'")
        sys.exit(1)
        
    if not args.bed_file.is_file():
        print(f"Error: BED file not found at '{args.bed_file}'")
        sys.exit(1)

    # --- Create Output Directory ---
    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"All outputs will be saved in: {args.output_dir.resolve()}")

    # --- Execute Workflow ---
    sambamba_output = run_sambamba(args, args.output_dir)
    process_results(sambamba_output, args, args.output_dir)


if __name__ == "__main__":
    main()