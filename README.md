# RNA-seq and differential expression 2022-10

Bash scripts to generate commands to submitted via slurm job scheduler for high-throughput differential gene expression analysis.

## Tools Used

Cufflinks - gffread

BWA

BEDTools

SAMtools

HISAT2

STAR

## Project Setup

Perform initial setup for project directories and files to make subsequent filtering and alignment run smoother.

Move to your Research_Project directory.

```bash
cd /lustre/projects/Research_Project-T110796
```


### Download seq data in new subproject folder
```bash
# Create new folder and move into folder
mkdir Project_10762
cd Project_10762

# Run V0268 Raw Reads
curl URL.tar | tar -xv

# Run V0268 QC Files
curl URL.tar | tar -xv

# Run V0268 Trimmed Reads
curl URL.tar | tar -xv
```


Make a directory tree that is compatible with job generation and job scripts. First you will need to rename the 
main project directory after your project.
This next command will make almost the entire tree. The -m 770 option will add full directory permissions to owner and group users.

```bash
mkdir -m 770 -p OCT22_RNA_seq/rRNA_filtering/{filtered_fastqs,STAR_alignment,HISAT2_alignment}
```

To make scripts look cleaner it is useful to link files that are stored else where into the project analysis directory.
This is called a symbolic link. To link to the fastp_trimmed directory that contains the trimmed fastq files received 
from the sequencing centre use this command. Change the directory names as required. I ran this from the Research_Project-T110796 directory.

```bash
ln -s /lustre/projects/Research_Project-T110796/Project_10762/V0268/11_fastp_trimmed/ /lustre/projects/Research_Project-T110796/Project_10762/ages_alignments/fastqs

/lustre/projects/Research_Project-T110796/Project_10762/ages_alignments

``` 

The reference genome for *Pristionchus pacificus* should already be in the el_paco_ref directory. If a new version is released upload them
into this directory. You will then need to index the new version for HISAT2 or STAR using the scripts provided.

### Downloading Scripts 
To obtain the scripts from the GitHub repository first move to diff_expr_scripts/ then load the git module and use the following command.

```bash
git clone https://github.com/Harry-Pollitt/RNA-Seq-Projects.git . # . is current working directory 

git pull origin  # should download updated scripts if needed
```



Once you have made the directory tree and linked the fastq directory, it should look something like this.
* is the symbolically linked directory

```bash
cd /lustre/projects/Research_Project-T110796/Project_10762
|
└── Project_10762/
    ├── V0268/
    │   ├── 01_raw_reads
    │   ├── 09_QC_reports
    │   └── 11_fastp_trimmed*
    |
    ├── ages_alignments/
    │   ├── fastqs*
    │   └── rRNA_filtering/
    │       ├── filtered_fastqs
    │       ├── STAR_alignment
    │       └── HISAT2_alignment
    |
    ├── el_paco_ref/
    │   ├── El_Paco_V3_gene_annotations.gff3
    │   ├── El_Paco_genome.fa 
    │   ├── El_Paco_V3_gene_annotations.gtf 
    │   ├── STAR_Index
    │   └── HISAT_Index
    |
    └── diff_expr_scripts/
        └── here be scripts...
```

## Script Usage

Check each script and change the #SBATCH parameters and other lines as necessary.

--- Creating Index files for alignments ---

```bash
sbatch hisat-indexing-job.sh
sbatch star-indexing-job.sh
```
This only needs to be done once for each reference genome, reuse the index for each alignment.

--- Aligning rRNA reads ---

# Create absolute symbolic link from diff_expr scripts in filtered_fastqs to run script in this directory

```bash
ln -s /lustre/projects/Research_Project-T110796/Project_10762/diff_expr_scripts/bwa-rRNA-job-script.sh rRNA-filtering.sh
```
In the job creation loop check that all file paths are correct. If in doubt use absolute path.

```bash
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

```bash
deseq2-edgeR-protocol.R
```

To use this R script you will need to download your HISAT2/STAR aligned bam files and the gtf annotation file and place them into your working directory for RStudio.
Then load this script and work through each step. Modifying the script to your needs. 

## Contributing

Harry Pollitt 

Email: hap39@aber.ac.uk

Rebekah White 

Email: rw617@exeter.ac.uk

Cameron Weadick

Email: c.weadick@exeter.ac.uk

## Acknowledgement

The authors would like to acknowledge the use of the University of Exeter High-Performance Computing (HPC) facility in
in carrying out this work.
