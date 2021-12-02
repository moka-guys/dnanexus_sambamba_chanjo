import os
import pytest
import re
import subprocess
import time

# TODO split monolithic test into many

# TODO add to config file

app_version_to_test = "chanjo_sambamba_coverage_v1.13"


# TODO remove single necessity for single quotes

# TODO Replace list with dictionary

# TODO Refactor code to replace sleep function

# Combination of potential Sambamba settings
# test_name, merge_overlapping, exclude_failed_qc, exclude_duplicate, sambamba_flags, sambamba_filter
sambamba_settings = [
    [
        "test_111_standardSettings",
        "true",
        "true",
        "true",
        "",
        "",
    ],
    [
        "test_000_additionalFilter",
        "false",
        "false",
        "false",
        "",
        "first_of_pair",
    ],
    [
        "test_000_additionalFilter2",
        "false",
        "false",
        "false",
        "",
        '"not mate_is_unmapped and first_of_pair"',
    ],
    [
        "test_111_additionalFilter",
        "true",
        "true",
        "true",
        "",
        "first_of_pair",
    ],
    [
        "test_111_additionalFlags",
        "true",
        "true",
        "true",
        "--annotate",
        "",
    ],
    [
        "test_011_additionalFlags",
        "false",
        "true",
        "true",
        "-m",
        "",
    ],
    [
        "test_000_booleanFlags",
        "false",
        "false",
        "false",
        "",
        "",
    ],
    [
        "test_100_booleanFlags",
        "true",
        "false",
        "false",
        "",
        "",
    ],
    [
        "test_101_booleanFlags",
        "true",
        "false",
        "true",
        "",
        "",
    ],
    [
        "test_001_booleanFlags",
        "false",
        "false",
        "true",
        "",
        "",
    ],
    [
        "test_110_booleanFlags",
        "true",
        "true",
        "false",
        "",
        "",
    ],
    [
        "test_010_booleanFlags",
        "false",
        "true",
        "false",
        "",
        "",
    ],
]

def check_job_status(num_of_jobs_to_monitor):
    # Typically set num_of_jobs_to_monitor to number of conditions in 'sambamba_settings' defined above
    completed_jobs_ids = subprocess.run(
        "dx find jobs --num-results " + str(num_of_jobs_to_monitor) + " --project project-G3f0qj804fp7F74V9BVgJQZf --brief --state done",
        check=True,
        shell=True,
        stdout=subprocess.PIPE
    )
    return completed_jobs_ids.stdout.decode('utf-8').splitlines()


@pytest.fixture()
def run_dx_jobs():
    """Run sambamba app with various combination of inputs"""
    # Change to project 003_211202_Sambamba_Testing
    subprocess.run(
        "dx select project-G3f0qj804fp7F74V9BVgJQZf",
        check=True,
        shell=True,
    )
    dx_job_ids = []
    for opt_list in sambamba_settings:
        (
            test_name,
            merge_overlapping_boolean,
            exclude_failed_qc_boolean,
            exclude_duplicate_boolean,
            sambamba_flags,
            sambamba_filter,
        ) = opt_list
        stream = subprocess.run(
            f'dx run 003_211202_Sambamba_Testing:{app_version_to_test} \
        --yes --name {test_name} \
        -i sambamba_bed=project-ByfFPz00jy1fk6PjpZ95F27J:file-FgQ0G1Q0jy1j1221FzBqG089 \
        -i bamfile=project-G0FykV00GfYYJkF2B35Q27Yb:file-G0X31380zQkqYz293jk7VgJV \
        -i bam_index=project-G0FykV00GfYYJkF2B35Q27Yb:file-G0X31380zQkz3B403gvKK1y0 \
        -i merge_overlapping_mate_reads={merge_overlapping_boolean} \
        -i exclude_failed_quality_control={exclude_failed_qc_boolean} \
        -i exclude_duplicate_reads={exclude_duplicate_boolean} \
        -i coverage_level="150" \
        -i min_base_qual="25" \
        -i min_mapping_qual="20" \
        -i additional_sambamba_flags={sambamba_flags} \
        -i additional_filter_commands={sambamba_filter} \
        --brief',
            check=True,
            shell=True,
            text=True,
            capture_output=True,
        )
        dx_job_id = str(stream.stdout).rstrip()
        dx_job_ids.append(dx_job_id)
    """
    # Wait for all jobs to finish before running tests
    waiting_for_these_jobs = " ".join(dx_job_ids)
    process = subprocess.call(
        "dx wait {waiting_for_these_jobs}", shell=True, stdout=subprocess.PIPE
    )
    """

    status_cmd = ("dx find jobs --num-results 12 --project " + 'project-G3f0qj804fp7F74V9BVgJQZf'  " --brief --state done")
    time.sleep(10)
    
    return dx_job_ids


def test_dx_inputs(run_dx_jobs):
    """Test combination of inputs to DNANEXUS_SAMBAMBA_CHANJO app"""
    # Get job ids
    dx_job_ids = run_dx_jobs

    # test_111_standardSettings
    stream = subprocess.run(
        f"dx watch {dx_job_ids[0]}",
        check=True,
        shell=True,
        text=True,
        capture_output=True,
    )
    dx_job_log = str(stream.stdout)
    pattern = 'sambamba depth region -L bedfile -t 2 -T 150 -q 25  -m  -F "mapping_quality >= 20 and not failed_quality_control and not duplicate'
    assert re.search(pattern, dx_job_log)
    assert re.search("(done)", dx_job_log)

    # test_000_additionalFilter
    stream = subprocess.run(
        f"dx watch {dx_job_ids[1]}",
        check=True,
        shell=True,
        text=True,
        capture_output=True,
    )
    dx_job_log = str(stream.stdout)
    pattern = 'sambamba depth region -L bedfile -t 2 -T 150 -q 25  -F "mapping_quality >= 20 and first_of_pair'
    assert re.search(pattern, dx_job_log)
    assert re.search("(done)", dx_job_log)

    # test_000_additionalFilter2
    stream = subprocess.run(
        f"dx watch {dx_job_ids[2]}",
        check=True,
        shell=True,
        text=True,
        capture_output=True,
    )
    dx_job_log = str(stream.stdout)
    pattern = 'sambamba depth region -L bedfile -t 2 -T 150 -q 25  -F "mapping_quality >= 20 and not mate_is_unmapped and first_of_pair'
    assert re.search(pattern, dx_job_log)
    assert re.search("(done)", dx_job_log)

    # test_111_additionalFilter
    stream = subprocess.run(
        f"dx watch {dx_job_ids[3]}",
        check=True,
        shell=True,
        text=True,
        capture_output=True,
    )
    dx_job_log = str(stream.stdout)
    pattern = 'sambamba depth region -L bedfile -t 2 -T 150 -q 25  -m  -F "mapping_quality >= 20 and not failed_quality_control and not duplicate and first_of_pair'
    assert re.search(pattern, dx_job_log)
    assert re.search("(done)", dx_job_log)

    # test_111_additionalFlags
    stream = subprocess.run(
        f"dx watch {dx_job_ids[4]}",
        check=True,
        shell=True,
        text=True,
        capture_output=True,
    )
    dx_job_log = str(stream.stdout)
    pattern = 'sambamba depth region -L bedfile -t 2 -T 150 -q 25  -m  --annotate -F "mapping_quality >= 20 and not failed_quality_control and not duplicate'
    assert re.search(pattern, dx_job_log)
    assert re.search("(done)", dx_job_log)

    # test_011_additionalFlags
    stream = subprocess.run(
        f"dx watch {dx_job_ids[5]}",
        check=True,
        shell=True,
        text=True,
        capture_output=True,
    )
    dx_job_log = str(stream.stdout)
    pattern = 'sambamba depth region -L bedfile -t 2 -T 150 -q 25  -m -F "mapping_quality >= 20 and not failed_quality_control and not duplicate'
    assert re.search(pattern, dx_job_log)
    assert re.search("(done)", dx_job_log)

    # test_000_booleanFlags
    stream = subprocess.run(
        f"dx watch {dx_job_ids[6]}",
        check=True,
        shell=True,
        text=True,
        capture_output=True,
    )
    dx_job_log = str(stream.stdout)
    pattern = (
        'sambamba depth region -L bedfile -t 2 -T 150 -q 25  -F "mapping_quality >= 20'
    )
    assert re.search(pattern, dx_job_log)
    assert re.search("(done)", dx_job_log)

    # test_100_booleanFlags
    stream = subprocess.run(
        f"dx watch {dx_job_ids[7]}",
        check=True,
        shell=True,
        text=True,
        capture_output=True,
    )
    dx_job_log = str(stream.stdout)
    pattern = 'sambamba depth region -L bedfile -t 2 -T 150 -q 25  -m  -F "mapping_quality >= 20'
    assert re.search(pattern, dx_job_log)
    assert re.search("(done)", dx_job_log)

    # test_101_booleanFlags
    stream = subprocess.run(
        f"dx watch {dx_job_ids[8]}",
        check=True,
        shell=True,
        text=True,
        capture_output=True,
    )
    dx_job_log = str(stream.stdout)
    pattern = 'sambamba depth region -L bedfile -t 2 -T 150 -q 25  -m  -F "mapping_quality >= 20 and not duplicate'
    assert re.search(pattern, dx_job_log)
    assert re.search("(done)", dx_job_log)

    # test_001_booleanFlags
    stream = subprocess.run(
        f"dx watch {dx_job_ids[9]}",
        check=True,
        shell=True,
        text=True,
        capture_output=True,
    )
    dx_job_log = str(stream.stdout)
    pattern = 'sambamba depth region -L bedfile -t 2 -T 150 -q 25  -F "mapping_quality >= 20 and not duplicate'
    assert re.search(pattern, dx_job_log)
    assert re.search("(done)", dx_job_log)

    # test_110_booleanFlags
    stream = subprocess.run(
        f"dx watch {dx_job_ids[10]}",
        check=True,
        shell=True,
        text=True,
        capture_output=True,
    )
    dx_job_log = str(stream.stdout)
    pattern = 'sambamba depth region -L bedfile -t 2 -T 150 -q 25  -m  -F "mapping_quality >= 20 and not failed_quality_control'
    assert re.search(pattern, dx_job_log)
    assert re.search("(done)", dx_job_log)

    # test_010_booleanFlags
    stream = subprocess.run(
        f"dx watch {dx_job_ids[11]}",
        check=True,
        shell=True,
        text=True,
        capture_output=True,
    )
    dx_job_log = str(stream.stdout)
    pattern = 'sambamba depth region -L bedfile -t 2 -T 150 -q 25  -F "mapping_quality >= 20 and not failed_quality_control'
    assert re.search(pattern, dx_job_log)
    assert re.search("(done)", dx_job_log)
