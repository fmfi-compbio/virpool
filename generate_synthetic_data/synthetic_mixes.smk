import json
import subprocess
import sys
import time
from collections import defaultdict
from os.path import join

import pandas as pd
import pysam
import requests
import yaml

raw_bam_dir = "raw_bam"
filtered_bam_dir = "filtered_bam"
mixes_dir = "mixes"
raw_data_filename = "raw_data_info.tsv"
mixes_info_filename = "mixes.yaml"

raw_data_info = pd.read_csv(raw_data_filename, delimiter="\t")
print(raw_data_info)
with open(mixes_info_filename) as f:
    mixes_info = yaml.safe_load(f)

def get_ftp_link_by_ena_id_api(ena_id):
    time.sleep(0.5)
    url = f'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={ena_id}' \
          f'&download=false&fields=submitted_ftp&result=read_run&format=JSON'
    response = requests.get(url=url)
    if response.status_code != 200:
        raise Exception(f"Return code is not OK 200: {response.status_code}!")
    data = json.loads(response.text)
    if len(data) != 1:
        raise Exception(f"The ENA ID should be unique, got {len(data)} lines instead!")
    ftp_link = data[0]['submitted_ftp']
    return ftp_link

def get_ena_id_by_internal_id(internal_id, raw_data_info):
    result = raw_data_info[raw_data_info['internal_id'] == internal_id]['ena_id'].values[0]
    return result

def is_good_alignment(al: pysam.AlignedSegment):
    return not any((al.is_secondary,
                    al.is_supplementary,
                    al.is_unmapped,
                    al.mate_is_unmapped,))

rule download_ena_bam:
    output: temp(join(raw_bam_dir, "{sample}.bam"))
    resources:
        download=1
    run:
        ena_id = get_ena_id_by_internal_id(
            internal_id=wildcards.sample,
            raw_data_info=raw_data_info)

        ftp_link = get_ftp_link_by_ena_id_api(ena_id)

        return_code = subprocess.call(['wget', '-O', output[0], ftp_link])
        if return_code != 0:
            raise Exception(f"Something went wrong during download of '{ftp_link}'!")

rule filter_bam:
    input:
        bam=join(raw_bam_dir, "{sample}.bam"),
        genome="../genome/genome.fa"
    output: join(filtered_bam_dir, "{sample}.bam")
    shell: """samtools view -hb -F 2316 -T {input.genome} {input.bam} > {output}"""

rule create_mix:
    input:
        bams=lambda w: [join(filtered_bam_dir, f"{sample}.bam") for sample in mixes_info[f"{w.mix}"]]
    params:
        samples=lambda w: mixes_info[f"{w.mix}"]
    output:
        mix=join(mixes_dir, "{mix}.bam"),
        rcounts=join(mixes_dir, "{mix}.rcounts")
    run:
        print(f"Starting with mix {wildcards.mix}...")
        sys.stdout.flush()
        rcounts = defaultdict(int)
        with pysam.AlignmentFile(input.bams[0]) as template, \
                pysam.AlignmentFile(output.mix, "wb", template=template) as output_bam:
            for sample, bam_filename in zip(params.samples, input.bams):
                variant = raw_data_info[raw_data_info['internal_id'] == sample]['variant'].values[0]

                with pysam.AlignmentFile(bam_filename) as bam:
                    for alignment in filter(is_good_alignment, bam.fetch(until_eof=True)):
                        rcounts[variant] += 1
                        output_bam.write(alignment)

        with open(output.rcounts, "w") as f:
            for variant, count in rcounts.items():
                print(f"{variant}\t{count}", file=f)

rule all_mixes:
    input: [join(mixes_dir, f"{mix}.bam") for mix in mixes_info]

rule all:
    input:
        rules.all_mixes.input
