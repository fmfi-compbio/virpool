# VirPool

## Installation

```bash
conda install -c conda-forge mamba
mamba env update -f requirements.yml -n virpool
```

## Usage

See `examples/` for examples of usage.

For a detailed interface description, use `virpool -h`.

TODO description of usage

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
