# RNA-seq and differential expression 2022-11

Bash scripts to generate commands to submitted via slurm job scheduler for high-throughput differential gene expression analysis. For unexpected error messages, check out handycode.txt and handycode2.txt for help.

## Tools Used

Cufflinks - gffread

BWA

BEDTools

SAMtools

HISAT2

MultiQC containing Picard, Preseq, Qualimap, RSeQC, and Samtools alignment metrics.

## Project Setup

Perform initial setup for project directories and files to make subsequent filtering and alignment run smoother.

Move to your Research_Project directory.

```bash
cd /lustre/projects/Research_Project-T110796
```


### Download seq data in new subproject folder
```bash
# Create new folder and move into folder
mkdir nf_10762/fastqs/
cd nf_10762/fastqs/

# Get QC Files
curl URL.tar | tar -xv

# Get Trimmed Reads
curl URL.tar | tar -xv
```


The reference genome for *Pristionchus pacificus* should already be in the el_paco directory. If a new version is released upload them
into this directory. You will then need to index the new version for HISAT2 or STAR using the scripts provided.



### Downloading Scripts 
(needs updating)
To obtain the scripts from the GitHub repository first move to diff_expr_scripts/ then load the git module and use the following command.

```bash
git clone https://github.com/Harry-Pollitt/RNA-Seq-Projects.git . # . is current working directory 

git pull origin  # should download updated scripts if needed
```



Set up the directories to look like this.

```bash
cd /lustre/projects/Research_Project-T110796/
|
└── nf_10762/
    ├── ages_nf/
    │   ├── samplesheet.csv
    │   └── ages_nf_jobs.sh
    |
    ├── dr_nf/
    │   ├── samplesheet.csv
    |   └── dr_nf_jobs.sh
    |
    ├── nys_nf/
    │   ├── samplesheet.csv
    │   └── nys_nf_jobs.sh
    |
    ├── el_paco/
    │   ├── El_Paco_genome_V3.gtf
    │   └── El_Paco_genome_V3.fa
    |
    ├── fastqs/
    │   ├── negativecontrol.fastq.gz
    │   ├── sample1.fastq.gz
    │   ├── sample2.fastq.gz
    │   └── Etc...      
    |
    └── generate_samplesheet.sh 


```

## Script Usage

Check each script and change the #SBATCH parameters and other lines as necessary.

In the directory of each study, run the job script, e.g., 

```bash
cd dr_nf/
sbatch dr_nf_jobs.sh
```

Nextflow will filter rRNA reads and generate aligned bam files.


--- featureCounts and DESeq2 in RStudio ---

```bash
deseq2-edgeR-protocol.R
```

To use this R script you will need to download your HISAT2 aligned bam files and the gtf annotation file and place them into your working directory for RStudio.

Then load this script and work through each step, modifying the script to your needs. 

Using your list of DEGs generated, find Ppa-specific functions using

```bash
degs_align_with_Sun.R
```
Then to get a summary of DEG biotypes


```bash
biotype_chart_dr.R
```



## Contributing

Harry Pollitt

Email: hap39@aber.ac.uk

Rebekah White

Email: rw617@exeter.ac.uk

Cameron Weadick (PI)

Email: c.weadick@exeter.ac.uk

## Acknowledgement

The authors would like to acknowledge the use of the University of Exeter High-Performance Computing (HPC) facility in
in carrying out this work.


## References

Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C (2017) Nextflow enables reproducible computational workflows. *Nature Biotechnology* 35(4):316-319.

Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S (2020) The nf-core framework for community-curated bioinformatics pipelines. *Nature Biotechnology* 38: 276–278.

See pipeline tool citations here: https://github.com/nf-core/rnaseq/blob/master/CITATIONS.md
