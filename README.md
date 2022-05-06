# RNA-seq and differential expression 2022-03

Bash scripts to generate commands to submit to a slurm job scheduler for high-throughput differential gene expression analysis.

## Install

Tools used in this project:

BWA

BEDTools

SAMtools

HISAT2

STAR

## Usage

--- Creating Index files for alignments ---

```bash
sh hisat-star-indexing.sh
```
This only needs to be done once for each reference genome, reuse for each alignment.

--- Aligning rRNA reads ---

``bash
sh generate-bwa-rRNA-commands.sh
sbatch bwa-rRNA-job-script.sh
```

--- Converting non-rRNA bams to fastqs ---

```bash
sh generate-bamtofastq-commands.sh
sbatch bamtofastq-job-script.sh
```

--- Align filtered reads to genome with HISAT2 ---

```bash
sh generate-hisat2-commands.sh
sbatch hisat2-job-script.sh
```

--- Align filered reads to genome with STAR ---

```bash
sh generate-star-commands.sh
sbatch star-job-script.sh
```

--- featureCounts and DESeq2 in RStudio ---

`deseq2-edgeR-protocol.R`

This is not automated, change code as necessary.

## Contributing

Harry Pollitt, Rebekah White

## Acknowledgement

The authors would like to acknowledge the use of the University of Exeter High-Performance Computing (HPC) facility in
in carrying out this work.
