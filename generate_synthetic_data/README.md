## Data for synthetic experiments

### Folders:

- `synthetic-ont-long/`
- `synthetic-ont-short/`
- `synthetic-illumina/`

### Run:
```shell
conda activate covid-pooling2
cd <experiment_folder>
snakemake all -j4 --res download=2  --snakefile ../synthetic_mixes.smk
```

### Config files:

- `raw_data_info.tsv` - description of downloaded files (ENA run accession number, internal id, variant)
- `mixes.yaml` - description of desired mixes (mix name, list of internal ids)

### Prepare all data:
```shell
conda activate covid-pooling2
# cd problematic_sites
# wget -O 'masking.vcf' 'https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf'
cd ../synthetic-illumina
snakemake all -j4 --res download=2  --snakefile ../synthetic_mixes.smk
cd ../synthetic-ont-long
snakemake all -j4 --res download=2  --snakefile ../synthetic_mixes.smk
cd ../synthetic-ont-short
snakemake all -j4 --res download=2  --snakefile ../synthetic_mixes.smk
# cd ../synthetic-gisaid
# bash estimate_profile.sh # takes 1h 20min 
```

### Prepare single sample synthetic data
```shell
conda activate covid-pooling2
cd synthetic-ont-long-single
snakemake all -j4 --res download=2  --snakefile ../synthetic_mixes.smk
cd ../synthetic-ont-short-single
snakemake all -j4 --res download=2  --snakefile ../synthetic_mixes.smk
cd ../synthetic-illumina-single
snakemake all -j4 --res download=2  --snakefile ../synthetic_mixes.smk

```

### Prepare data for synthetic experiments of dominating variant setting
```shell
conda activate covid-pooling2
cd synth-ol-by-variant
snakemake all -j4 --res download=2  --snakefile ../synthetic_mixes.smk
cd ../synth-os-by-variant
snakemake all -j4 --res download=2  --snakefile ../synthetic_mixes.smk
cd ../synth-il-by-variant
snakemake all -j4 --res download=2  --snakefile ../synthetic_mixes.smk

cd synth-ol-dominating
snakemake all -j48 --snakefile ../mix_with_coverages.smk --configfile config.yaml -F --rerun-incomplete
cd ../synth-os-dominating;
snakemake all -j48 --snakefile ../mix_with_coverages.smk --configfile config.yaml -F --rerun-incomplete;
cd ../synth-il-dominating;
snakemake all -j48 --snakefile ../mix_with_coverages.smk --configfile config.yaml -F --rerun-incomplete

```

### Prepare data for synthetic experiments of domination variant setting alpha vs delta, delta vs omicron

```shell
conda activate covid-pooling2
cd synth-dom-il-alpha-delta-by-variant
snakemake all -j4 --res download=2  --snakefile ../synthetic_mixes.smk --keep-going --rerun-incomplete
cd ../synth-dom-os-alpha-delta-by-variant
snakemake all -j4 --res download=2  --snakefile ../synthetic_mixes.smk --keep-going --rerun-incomplete
cd ../synth-dom-il-delta-omicron-by-variant
snakemake all -j4 --res download=2  --snakefile ../synthetic_mixes.smk --keep-going --rerun-incomplete
cd ../synth-dom-os-delta-omicron-by-variant
snakemake all -j4 --res download=2  --snakefile ../synthetic_mixes.smk --keep-going --rerun-incomplete
cd ../

cd synth-dom-il-alpha-delta-by-variant
snakemake all -j48 --snakefile ../mix_with_coverages.smk --configfile config.yaml --rerun-incomplete
cd ../synth-dom-os-alpha-delta-by-variant
snakemake all -j48 --snakefile ../mix_with_coverages.smk --configfile config.yaml --rerun-incomplete
cd ../synth-dom-il-delta-omicron-by-variant
snakemake all -j48 --snakefile ../mix_with_coverages.smk --configfile config.yaml --rerun-incomplete
cd ../synth-dom-os-delta-omicron-by-variant
snakemake all -j48 --snakefile ../mix_with_coverages.smk --configfile config.yaml --rerun-incomplete
cd ../



```