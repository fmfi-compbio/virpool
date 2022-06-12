# serious example

```bash

# download reads from Nice, France by Rios et al. 2021 (FABRON, JAN 2021)
wget -O reads.fastq.gz --retry-connrefused ftp.sra.ebi.ac.uk/vol1/fastq/SRR152/052/SRR15275952/SRR15275952_1.fastq.gz

# align reads to genome using minimap2
minimap2 -ax map-ont  -t 4 ../../data/genome.fa reads.fastq.gz | samtools sort -o reads.bam

conda activate virpool
../../src/virpool -v ../../data/niw_2019-01-01_2021-03-31.tsv -g ../../data/genome.fa -m ../../data/ont-short.masking.vcf reads.bam
```
