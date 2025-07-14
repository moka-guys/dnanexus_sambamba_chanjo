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
from datetime import datetime

# --- Configuration ---
SAMBAMBA_OUTPUT_FILENAME = "sambamba_output.bed"
EXON_REPORT_FILENAME = "exon_level.txt"
GENE_REPORT_FILENAME = "gene_level.txt"

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

    # Execute sambamba with suppressed stderr to avoid verbose warnings
    print(f"🚀 Running command: {' '.join(command)}")
    try:
        # Conditionally suppress stderr based on verbose flag
        stderr_setting = None if args.verbose else subprocess.DEVNULL
        
        process = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=stderr_setting,  # Show stderr if verbose, suppress otherwise
            text=True,
        )
        
        # Write stdout to output file
        with open(output_bed, "w") as f_out:
            f_out.write(process.stdout)
            
    except FileNotFoundError as e:
        print(f"❌ Error: Command 'sambamba' not found.")
        print("Please ensure sambamba is installed and in your system's PATH.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"❌ Sambamba command failed with return code: {e.returncode}")
        if not args.verbose:
            print("💡 Tip: Use --verbose flag to see sambamba's detailed error messages")
        print("This may indicate an issue with the input files or command parameters.")
        sys.exit(1)
        
    print(f"✅ Sambamba analysis complete. Output saved to: {output_bed}")
    return output_bed


def process_results(sambamba_bed: Path, args: argparse.Namespace, output_dir: Path, panel_genes: Optional[Set[str]] = None):
    """
    Parses the outputs from sambamba to create final reports.
    This function is a Python 3 rewrite of the logic in `read_chanjo.py`.
    
    Args:
        sambamba_bed: Path to sambamba output BED file
        args: Command line arguments
        output_dir: Output directory
        panel_genes: Optional set of gene symbols to filter results (only used for PDF, not text files)
    """
    print("\n--- Processing and Generating Final Reports ---")
    
    if panel_genes:
        print(f"   - Panel filtering will be applied to PDF report only ({len(panel_genes)} genes)")
    
    # --- Part 1: Create gene-level statistics from sambamba output (ALL GENES) ---
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
                
                # NO PANEL FILTERING HERE - include all genes in text report
                if gene_symbol not in gene_stats:
                    gene_stats[gene_symbol] = {
                        'coverage_values': [], 
                        'completeness_values': [], 
                        'exon_lengths': [],
                        'total_length': 0
                    }
                
                gene_stats[gene_symbol]['coverage_values'].append(mean_coverage)
                gene_stats[gene_symbol]['completeness_values'].append(completeness)
                gene_stats[gene_symbol]['exon_lengths'].append(exon_length)
                gene_stats[gene_symbol]['total_length'] += exon_length

    # Generate gene-level report (ALL GENES) with length-weighted averages
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
        
    # --- Part 2: Create an exon-level report for regions below 100% coverage (ALL EXONS) ---
    exon_report_path = output_dir / EXON_REPORT_FILENAME
    exon_count = 0
    
    with open(sambamba_bed, 'r') as f_in, open(exon_report_path, 'w') as f_out:
        f_out.write("gene_symbol\ttranscript\tentrez_id\tchr\tstart\tstop\tmean_coverage\tcompleteness\n")
        
        for line in f_in:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
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
                
                # NO PANEL FILTERING HERE - include all exons in text report
                f_out.write(f"{gene_symbol}\t{transcript}\t{entrez_id}\t{chrom}\t{start}\t{stop}\t{mean_coverage}\t{completeness}\n")
                
    print(f"   - Generated exon-level report: {exon_report_path} ({exon_count} exons with <100% coverage)")
 
    print(f"\n✅ All reports generated successfully in: {output_dir}")


def generate_pdf_report(sambamba_bed: Path, args: argparse.Namespace, output_dir: Path, panel_genes: Optional[Set[str]] = None):
    """
    Generates a comprehensive PDF coverage report.
    
    Args:
        sambamba_bed: Path to sambamba output BED file
        args: Command line arguments
        output_dir: Output directory
        panel_genes: Optional set of gene symbols for panel filtering
    """
    try:
        from reportlab.lib.pagesizes import letter, A4
        from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak, BaseDocTemplate, PageTemplate, Frame, KeepTogether
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.units import inch
        from reportlab.lib import colors
        from reportlab.lib.enums import TA_CENTER, TA_RIGHT, TA_LEFT
        from reportlab.graphics.shapes import Line, Drawing
        from reportlab.platypus.flowables import HRFlowable
    except ImportError:
        print("   - Warning: reportlab not installed. Install with: pip install reportlab")
        print("   - Skipping PDF report generation")
        return None

    print("\n--- Generating PDF Report ---")
    
    # Generate sample ID and report info
    sample_id = os.path.splitext(os.path.basename(args.bam_file))[0]
    current_date = datetime.now().strftime("%Y-%m-%d %H:%M")
    pdf_filename = f"{sample_id}_coverage_report.pdf"
    pdf_path = output_dir / pdf_filename
    
    # Custom document template class with header on every page
    class HeaderDocTemplate(BaseDocTemplate):
        def __init__(self, filename, sample_id, current_date, **kwargs):
            BaseDocTemplate.__init__(self, filename, **kwargs)
            self.sample_id = sample_id
            self.current_date = current_date
            
            # Define frame for content (leaving space for header)
            frame = Frame(
                self.leftMargin, self.bottomMargin,
                self.width, self.height - 50,  # 50 points for header space
                id='main_frame'
            )
            
            # Create page template with header
            template = PageTemplate(id='header_template', frames=[frame], onPage=self.add_header)
            self.addPageTemplates([template])
        
        def add_header(self, canvas, doc):
            """Add header to every page"""
            canvas.saveState()
            
            # Header text styles - smaller fonts
            canvas.setFont("Helvetica-Bold", 8)
            header_left = f"{self.sample_id} coverage report"
            canvas.drawString(self.leftMargin, self.height + self.topMargin - 25, header_left)
            
            canvas.setFont("Helvetica", 8)
            header_right = self.current_date
            canvas.drawRightString(self.width + self.leftMargin, self.height + self.topMargin - 25, header_right)
            
            # Header line
            canvas.setStrokeColor(colors.black)
            canvas.setLineWidth(0.5)
            canvas.line(self.leftMargin, self.height + self.topMargin - 35, 
                       self.width + self.leftMargin, self.height + self.topMargin - 35)
            
            canvas.restoreState()
    
    # Create PDF document with custom header template
    doc = HeaderDocTemplate(str(pdf_path), sample_id, current_date, pagesize=letter,
                           rightMargin=72, leftMargin=72, topMargin=100, bottomMargin=72)
    
    # Define styles
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        'CustomTitle',
        parent=styles['Heading1'],
        fontSize=16,
        spaceAfter=30,
        alignment=TA_CENTER
    )
    
    section_style = ParagraphStyle(
        'Section',
        parent=styles['Heading2'],
        fontSize=12,  # Reduced from 14
        spaceAfter=12,
        spaceBefore=20
    )
    
    # Style for boxed content
    boxed_style = ParagraphStyle(
        'Boxed',
        parent=styles['Normal'],
        fontSize=9,  # Reduced from 11
        leftIndent=10,
        rightIndent=10,
        spaceAfter=10,
        spaceBefore=10,
        borderWidth=1,
        borderColor=colors.black,
        borderPadding=8
    )
    
    # Collect gene and exon data
    gene_stats = {}
    exon_data = []
    
    with open(sambamba_bed, 'r') as f_in:
        for line in f_in:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 11:
                gene_info = fields[6]
                entrez_id = fields[7]
                coords = fields[3]
                mean_coverage = float(fields[9]) if fields[9] != '.' else 0.0
                completeness = float(fields[10]) if fields[10] != '.' else 0.0
                
                if not entrez_id.isdigit():
                    continue
                
                gene_symbol = gene_info.split(';')[0] if ';' in gene_info else gene_info
                if entrez_id == "11200":
                    gene_symbol = "CHEK2"
                
                # Calculate exon length
                try:
                    chrom, start, stop = coords.split('-')
                    exon_length = int(stop) - int(start) + 1
                except ValueError:
                    chrom, start, stop = "N/A", "N/A", "N/A"
                    exon_length = 0  # Skip if can't parse coordinates
                
                if ";" in gene_info:
                    transcript = gene_info.split(";", 1)[1]
                else:
                    transcript = "N/A"
                
                # Apply panel filtering if enabled
                if panel_genes and gene_symbol not in panel_genes:
                    continue
                
                exon_data.append([gene_symbol, transcript, entrez_id, chrom, start, stop, 
                                f"{mean_coverage:.1f}", f"{completeness:.1f}"])
                
                # Only collect length data if we have valid coordinates
                if exon_length > 0:
                    if gene_symbol not in gene_stats:
                        gene_stats[gene_symbol] = {
                            'coverage_values': [], 
                            'completeness_values': [], 
                            'exon_lengths': [],
                            'total_length': 0
                        }
                    
                    gene_stats[gene_symbol]['coverage_values'].append(mean_coverage)
                    gene_stats[gene_symbol]['completeness_values'].append(completeness)
                    gene_stats[gene_symbol]['exon_lengths'].append(exon_length)
                    gene_stats[gene_symbol]['total_length'] += exon_length
    
    # Prepare gene summary data with length-weighted averages
    gene_summary_data = [['Gene Symbol', 'Mean Coverage', f'Completeness at {args.coverage_level}X', 'Exon Count']]
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
        gene_summary_data.append([gene_symbol, f"{avg_coverage:.1f}", f"{avg_completeness:.1f}", str(exon_count)])
    
    # Prepare exon summary data (only exons with <100% coverage)
    exon_summary_data = [['Gene Symbol', 'Transcript', 'Entrez ID', 'Chr', 'Start', 'Stop', 'Mean Coverage', 'Completeness']]
    for exon in exon_data:
        if float(exon[7]) < 100.0:  # completeness < 100%
            exon_summary_data.append(exon)
    
    # Build PDF content
    content = []
    
    # Page 1: Title and Summary
    content.append(Paragraph("Coverage Analysis Report", title_style))
    content.append(Spacer(1, 20))
    
    if panel_genes:
        genes_text = f"Panel genes analysed: {', '.join(sorted(panel_genes)) if len(panel_genes) <= 10 else ', '.join(sorted(list(panel_genes))[:10]) + f' (and {len(panel_genes) - 10} more)'}"
    else:
        # Get sorted list of all genes for display
        gene_list = sorted(gene_stats.keys())
        if len(gene_list) <= 10:
            genes_text = f"Genes analysed: {', '.join(gene_list)}"
        else:
            genes_text = f"Genes analysed: {', '.join(gene_list[:10])} (and {len(gene_list) - 10} more)"
    
    content.append(Paragraph(genes_text, boxed_style))
    content.append(Spacer(1, 20))
    
    # Gene coverage summary
    content.append(Paragraph("Gene Coverage Summary:", section_style))
    content.append(HRFlowable(width="100%", thickness=1, lineCap='round', color=colors.black, spaceBefore=3, spaceAfter=10))
    
    if len(gene_summary_data) > 1:
        # Set column widths to make table span full width
        col_widths = [2*inch, 1.5*inch, 1.5*inch, 1*inch]  # Adjust widths for 4 columns
        gene_table = Table(gene_summary_data, colWidths=col_widths)
        gene_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 9),  # Match exon table
            ('FONTSIZE', (0, 1), (-1, -1), 8),  # Match exon table
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('LEFTPADDING', (0, 0), (-1, -1), 6),
            ('RIGHTPADDING', (0, 0), (-1, -1), 6)
        ]))
        content.append(gene_table)
    else:
        content.append(Paragraph("No gene data available.", styles['Normal']))
    
    content.append(Spacer(1, 20))
    
    # Exon coverage summary  
    content.append(Paragraph("Exons with < 100% coverage @ 30x:", section_style))
    content.append(HRFlowable(width="100%", thickness=1, lineCap='round', color=colors.black, spaceBefore=3, spaceAfter=10))
    
    if len(exon_summary_data) > 1:
        # Set column widths for 8 columns to span full width
        col_widths = [1*inch, 1.2*inch, 0.8*inch, 0.6*inch, 0.8*inch, 0.8*inch, 1*inch, 1*inch]
        exon_table = Table(exon_summary_data, colWidths=col_widths)
        
        # Build table style with conditional formatting for failed exons
        table_style_commands = [
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 9),
            ('FONTSIZE', (0, 1), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('LEFTPADDING', (0, 0), (-1, -1), 6),
            ('RIGHTPADDING', (0, 0), (-1, -1), 6)
        ]
        
        # Add red highlighting for rows with completeness < 90% (very poor coverage)
        for i, exon in enumerate(exon_summary_data[1:], 1):  # Skip header row
            if len(exon) > 7:  # Make sure we have completeness data
                try:
                    completeness = float(exon[7])  # Completeness is in column 7
                    if completeness < 90.0:  # Highlight very poor coverage in red
                        table_style_commands.append(('BACKGROUND', (0, i), (-1, i), colors.lightcoral))
                        table_style_commands.append(('TEXTCOLOR', (0, i), (-1, i), colors.black))
                except (ValueError, IndexError):
                    pass  # Skip if we can't parse the completeness value
        
        exon_table.setStyle(TableStyle(table_style_commands))
        content.append(exon_table)
    else:
        content.append(Paragraph("All exons have 100% coverage at a minimum of 30x.", styles['Normal']))
    
    # Page 2: All exons
    content.append(PageBreak())
    content.append(Paragraph("Complete Exon Coverage Report", section_style))
    content.append(HRFlowable(width="100%", thickness=1, lineCap='round', color=colors.black, spaceBefore=3, spaceAfter=10))
    
    all_exons_data = [['Gene Symbol', 'Transcript', 'Entrez ID', 'Chr', 'Start', 'Stop', 'Mean Coverage', 'Completeness']]
    all_exons_data.extend(exon_data)
    
    if len(all_exons_data) > 1:
        # Set column widths for 8 columns to span full width
        col_widths = [1*inch, 1.2*inch, 0.8*inch, 0.6*inch, 0.8*inch, 0.8*inch, 1*inch, 1*inch]
        all_exons_table = Table(all_exons_data, colWidths=col_widths)
        
        # Build table style with conditional formatting for failed exons
        table_style_commands = [
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 9),
            ('FONTSIZE', (0, 1), (-1, -1), 7),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('LEFTPADDING', (0, 0), (-1, -1), 4),
            ('RIGHTPADDING', (0, 0), (-1, -1), 4)
        ]
        
        # Add red highlighting for rows with completeness < 90% (very poor coverage)
        for i, exon in enumerate(all_exons_data[1:], 1):  # Skip header row
            if len(exon) > 7:  # Make sure we have completeness data
                try:
                    completeness = float(exon[7])  # Completeness is in column 7
                    if completeness < 90.0:  # Highlight very poor coverage in red
                        table_style_commands.append(('BACKGROUND', (0, i), (-1, i), colors.lightcoral))
                        table_style_commands.append(('TEXTCOLOR', (0, i), (-1, i), colors.black))
                except (ValueError, IndexError):
                    pass  # Skip if we can't parse the completeness value
        
        all_exons_table.setStyle(TableStyle(table_style_commands))
        content.append(all_exons_table)
    
    # Build PDF
    doc.build(content)
    print(f"   - Generated PDF report: {pdf_path}")
    return pdf_path


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
    parser.add_argument("--merge-overlapping-mates", action='store_true', default=True, help="Sambamba: Count overlapping mate reads only once (-m flag).")
    parser.add_argument("--additional-filter", type=str, default="not (unmapped or secondary_alignment)", help="Sambamba: Additional filter commands (e.g., 'first_of_pair').")
    parser.add_argument("--additional-sambamba-flags", type=str, help="Sambamba: Other flags to pass, enclosed in quotes.")

    # --- General Settings ---
    parser.add_argument("-t", "--threads", type=int, default=os.cpu_count(), help="Number of threads to use.")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path.cwd() / "coverage_results", help="Directory to store all outputs.")

    # --- Panel Filtering ---
    parser.add_argument("--panel-filter", action='store_true', help="Enable panel-specific gene filtering based on BAM filename.")

    # --- PDF Report ---
    parser.add_argument("--pdf-report", action='store_true', help="Generate a comprehensive PDF coverage report.")

    # --- Verbosity ---
    parser.add_argument("--verbose", action='store_true', help="Show sambamba warnings and verbose output (default: suppressed).")

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
    process_results(sambamba_output, args, args.output_dir, panel_genes)
    
    if args.pdf_report:
        generate_pdf_report(sambamba_output, args, args.output_dir, panel_genes)
    
    # Clean up downloaded panel BED file if it exists
    if panel_genes:
        panel_bed_files = list(args.output_dir.glob("Pan*_CNV.bed"))
        for panel_bed_file in panel_bed_files:
            panel_bed_file.unlink()
            print(f"   - Cleaned up downloaded panel BED file: {panel_bed_file.name}")


if __name__ == "__main__":
    main()