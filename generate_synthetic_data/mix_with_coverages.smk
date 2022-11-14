import os
from os.path import join

#minor_proportions = config.get('minor_proportions', [0.001, 0.01, 0.02, 0.05, 0.10, 0.20, 0.5])
minor_proportions = config.get('minor_proportions', [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20])
total_coverage = config.get('total_coverage', 5_000)
max_sample_number = config.get('max_sample_number', 10)

source_dir = config['source_dir']
target_dir = config.get('target_dir', "")
samples = config.get('samples', [f.removesuffix(".bam") for f in os.listdir(source_dir) if f.endswith(".bam")])

rule subsample_bam_simple:
    input: join(source_dir, "{sample}.bam")
    params:
        coverage=lambda w: w['coverage']
    output:
        bam=temp(join(target_dir, "subsampled", "{sample}_{coverage}_{snum}.bam")),
        written_alignment=temp(join(target_dir, "subsampled_count", "{sample}_{coverage}_{snum}.txt"))
    shell: """python ../../src/subsample_bam.py {input} {output.bam} {params.coverage} --sampled-reads-count-filename {output.written_alignment} """

rule mix_major_and_minor:
    input:
        major=lambda w: join(target_dir, "subsampled", f"{w['sample1']}_{total_coverage * (1 - float(w['minor_proportion']))}_{w['snum']}.bam"),
        minor=lambda w: join(target_dir, "subsampled", f"{w['sample2']}_{total_coverage * float(w['minor_proportion'])}_{w['snum']}.bam")
    output: join(target_dir, "{sample1}_{sample2}_{minor_proportion}_{snum}.bam")
    shell: """python ../../src/merge_bams.py {input} -o {output} """

rule count_true_minor_proportion:
    input:
        major=lambda w: join(target_dir, "subsampled_count", f"{w['sample1']}_{total_coverage * (1 - float(w['minor_proportion']))}_{w['snum']}.txt"),
        minor=lambda w: join(target_dir, "subsampled_count", f"{w['sample2']}_{total_coverage * float(w['minor_proportion'])}_{w['snum']}.txt")
    output: join(target_dir, "true_minor_proportion", "{sample1}_{sample2}_{minor_proportion}_{snum}.txt")
    run:
        with open(input.major) as f:
            major_count = int(f.__next__().strip())
        with open(input.minor) as f:
            minor_count = int(f.__next__().strip())

        true_minor_proportion = minor_count / (minor_count + major_count)
        with open(output[0], "w") as f:
            print(true_minor_proportion, file=f)

rule all_mixes:
    input:
            [
                join(target_dir, f"{s1}_{s2}_{mp}_{snum}.bam")
                for s1 in samples
                for s2 in samples
                for mp in minor_proportions
                for snum in range(max_sample_number)
                if s1 != s2]

rule merge_true_minor_proportions:
    input:
            [join(target_dir, "true_minor_proportion", f"{s1}_{s2}_{mp}_{snum}.txt")
                for s1 in samples
                for s2 in samples
                for mp in minor_proportions
                for snum in range(max_sample_number)
                if s1 != s2]
    output:
        join(target_dir, "merged_true_minor_proportions.tsv")
    run:
        with open(output[0], "w") as output_f:
            print("id\tproportion", file=output_f)
            for filename in input:
                basename = os.path.basename(filename)
                sample_id = basename.removesuffix(".txt")
                with open(filename) as f:
                    true_proportion = f.__next__().strip()
                print(f"{sample_id}\t{true_proportion}", file=output_f)

rule all:
    input:
        rules.all_mixes.input,
        rules.merge_true_minor_proportions.output