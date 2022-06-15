# VirPool

## Installation

```bash
conda install -c conda-forge mamba # mamba is used to speed up the installation
mamba env update -f requirements.yml -n virpool
```

## Usage

### Simplest case

```bash
covid activate virpool
python src/virpool -v variants.tsv -g genome.fa -m masking.vcf -o results reads.bam
```

This will determine the proportions of variants presented in file `variants.tsv` among reads in file `reads.bam`.
The output will be stored in directory `results`, with three files:

1. `estimated_weights.yaml` - a YAML file with estimated weights for each variant in the
   `variants.tsv` file
2. `posterior_coverage.svg` - a plot showing the coverage by reads with high posterior probability
   of belonging to a particular variant
3. `significant_positions.tsv` - a plot showing the coverage by reads on positions with singificant
   changes betweeen a given variant and the rest.


### All options available

```
usage: virpool [-h] [-o OUTPUT_DIR] [-v VARIANTS] [-g GENOME] [-s SUBST_RATE]
               [-t THREADS] [-n TRIES] [-cs CLIPPING_START]
               [-ce CLIPPING_END] [-m MASKING_VCF] [--percentile PERCENTILE]
               reads

positional arguments:
  reads                 SAM/BAM file with reads from a mixed sample

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory (default: 'results-YYMMDDThhmmss')
  -v VARIANTS, --variants VARIANTS
                        TSV file with variant profiles (default:
                        'resources/variants.tsv')
  -g GENOME, --genome GENOME
                        FASTA file with the virus genome (default:
                        'resources/genome.fa')
  -s SUBST_RATE, --subst-rate SUBST_RATE
                        Substitution error rate (suggested value: 0.05 for
                        ONT reads, 0.001 for Illumina reads) (default: 0.05)
  -t THREADS, --threads THREADS
                        Number of threads to use during posteriors
                        computation (default: 10)
  -n TRIES, --tries TRIES
                        Number of reruns of the optimisation (default: 5)
  -cs CLIPPING_START, --clipping-start CLIPPING_START
                        Number of ignored bases from the start of the ref.
                        sequence (default: 0)
  -ce CLIPPING_END, --clipping-end CLIPPING_END
                        Number of ignored bases from the end of the ref.
                        sequence (default: 0)
  -m MASKING_VCF, --masking-vcf MASKING_VCF
                        VCF file with positions to be masked (e.g.
                        homeoplasic sites or primers); the sixth column
                        should contain keyword 'mask' (default: -)
  --percentile PERCENTILE
                        Punish gaps with coverage < `percentile` quantile for
                        the given variant (default: 0.0)

```

See `examples/` for examples of usage.



## Profile estimation

TODO add description of already provided profiles.

### Download GISAID data

TODO into folder `profiles/`, files `metadata.tsv` and `gisaid.cram`

### Estimate a profile

First, put your preferred list of variants into a file, e.g. `profiles/variant_list.txt`:

```
B.1.1.7
B.1.617.2
B.1.1.529
```

Then, the profile is calculated by a following command:

```bash
conda activate virpool

cd profiles
python ../src/profile_estimation.py variant_list.txt aliases.json metadata.tsv gisaid.cram \
    --reference ../data/genome.fa \
    -o profile_name.tsv \
    --date-min 2020-01-01 \
    --date-max 2022-12-31
```
