import datetime
import os
import pytest
import re
import requests
import subprocess
import time

working_project = "project-G3f0qj804fp7F74V9BVgJQZf"  # 003_Sambamba_Testing

MokaWES_config = {
    "benchmark_project": "project-G39450Q0X2181QPKKZZ5JYK4",
    "project_name": "002_210621_A01229_0024_AHFK5LDRXY_NGS407C2_WES70SKIN_WES70",
    "coverage_folder": "coverage",
    "bam_folder": "output",
    "sample_suffix": "_markdup.bam",
    "pipeline": "MokaWES_1.8",
    "sambamba_version": "sambamba and chanjo v1.9",
    "coverage_report_suffix": "chanjo_txt",
}


MokaAMP_config = {
    "benchmark_project": "project-G3Z06Z8008J23K9j3578q1PK",
    "project_name": "002_210708_M02353_0597_000000000-DCTW4_ONC21059",
    "coverage_folder": "MokaAMP_EGFR_trial/coverage",
    "bam_folder": "MokaAMP_EGFR_trial/output",
    "sample_suffix": "refined.primerclipped.bam",
    "pipeline": "MokaAMP_v1.6",
    "sambamba_version": "sambamba and chanjo v1.9",
    "coverage_report_suffix": "exon_level.txt",
}

MokaPipe_config = {
    "benchmark_project": "project-G3b0g7Q0B7Q27kK0444Vb6zQ",
    "project_name": "002_210709_NB552085_0135_AHVK5JAFX2_NGS417",
    "coverage_folder": "coverage",
    "bam_folder": "output",
    "sample_suffix": "refined.bam",
    "pipeline": "GATK3.5_2.12",
    "sambamba_version": "sambamba and chanjo v1.10",
    "coverage_report_suffix": "exon_level.txt",
}


def import_pan_numbers():
    """Import data from automate_demultiplex_config.py so that we have up-to-date info on the bedfiles to use"""
    url = "https://raw.githubusercontent.com/moka-guys/automate_demultiplex/master/automate_demultiplex_config.py"
    r = requests.get(url, allow_redirects=True)
    f = open("pan_number_lookup.py", "w")
    f.write(
        re.search(
            "#IMPORTANT: Lists below are used by the trend analysis scripts.*default_panel_properties = {",
            r.text,
            flags=re.DOTALL,
        ).group(0)
    )
    f.write("}\n")  # Close the matching dictionary
    f.close()
    f = open("pan_number_lookup.py", "a")
    f.write(
        re.search(
            "# override default panel settings.*# =====smartsheet API=====",
            r.text,
            flags=re.DOTALL,
        ).group(0)
    )
    f.close()
    from pan_number_lookup import panel_settings


def lookup_pan_number(sample_name):
    # Get Pan number from Sample name
    sample_pan = re.search("Pan[0-9]+", r.text, flags=re.DOTALL).group(0)
    # Lookup sambamba Pan number
    sambamba_pan_Num = panel_settings[sample_pan]["sambamba_bedfile"]
    return sambamba_pan_num


def generate_batch_cmd(config_dict):
    stream = subprocess.run(
        dx_batch_command=f"dx generate_batch_inputs \
        -ibamfile='(.*){config_dict[sample_suffix]}$' \
        -ibam_index='(.*){config_dict[sample_suffix]}.bai$' \
            --path='{config_dict[project_name]} --output_prefix {config_dict[pipeline]}"
    )

    stream = subprocess.run(
        dx_batch_command,
        check=True,
        shell=True,
        text=True,
        capture_output=True,
    )

    populate_batch_cmd_with_pan_numbers(config_dict[pipeline])

    """
    dx generate_batch_inputs -ibamfile='(.*)_markdup.bam$' -ibam_index='(.*)_markdup.bam.bai$' --path="002_210621_A01229_0024_AHFK5LDRXY_NGS407C2_WES70SKIN_WES70:output" --output_prefix "MokaWES"
    dx generate_batch_inputs -ibamfile='(.*)refined.primerclipped.bam$' -ibam_index='(.*)refined.primerclipped.bam.bai$' --path="002_210708_M02353_0597_000000000-DCTW4_ONC21059:MokaAMP_EGFR_trial/output" --output_prefix "MokaAMP"
    dx generate_batch_inputs -ibamfile='(.*)refined.bam$' -ibam_index='(.*)refined.bam.bai$' --path="002_210709_NB552085_0135_AHVK5JAFX2_NGS417:output" --output_prefix "MokaPIPE"
    """


def populate_batch_cmd_with_pan_numbers(dx_batch_file):
    amended_text = []
    for line in dx_batch_file:
        sambamba_bed_name = lookup_pan_number(line)
        # look up file ID
        sambamba_dx_id = f"dx find data --name {sambamba_bed_name} --path 001_ToolsReferenceData:Data/BED/"
        # dx find data --name "Pan468dataSambamba.bed" --path 001_ToolsReferenceData:Data/BED/
        amended_line = f""
        amended_text.append(amended_line)
    return amended_text


# TODO set variables for every input
def run_batch_command(sambamba_app, project_name, folder_path, sample_suffix):
    dx_batch_command = f"dx run --batch-tsv"


def download_benchmark_coverage(
    benchmark_project, benchmark_coverage_folder, destination_project
):
    f"dx cp"
    pass


def compare_coverage(benchmark_coverage_folder, test_coverage_folder):
    # Compare tested app coverage to benchmark app coverage
    f"diff -bur {benchmark_coverage_folder} {test_coverage_folder}"


# @pytest.fixture() TODO Uncomment
def compare_app_output(config_dict):
    import_pan_numbers()
    for pipeline, benchmark_config in config_dict():
        generate_batch_cmd(config_dict)


def create_working_directory(folder_path):
    now = datetime.datetime.today()
    nTime = now.strftime("%d-%m-%Y")
    dest = os.path.join(folder_path + nTime)
    if not os.path.exists(dest):
        os.makedirs(dest)  # creat dest dir
