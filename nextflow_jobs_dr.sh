#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=30:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-******* # research project to submit under
#SBATCH --nodes=4 # specify number of nodes
#SBATCH --ntasks-per-node=12 # specify number of processors per node
#SBATCH --mem=90G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=*****@******.ac.uk # email address


## Install miniconda3 (will not need to be done again if re-running jobs)
#wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
#bash ~/miniconda.sh -b -p
#rm ~/miniconda.sh
source ~/miniconda3/bin/activate
#printf '\n# add path to conda\nexport PATH="$HOME/miniconda3/bin:$PATH"\n' >> ~/.bashrc
conda activate

## CONDA NOTE
#If above gives error, try restarting command line, then rerun conda activate



## Load latest Java module and unset options
module load Java/1.8.0_74
unset _JAVA_OPTIONS
JAVA_OPTIONS="-Djava.io.tmpdir=$HOME/tmp"

## Load Samtools
module load SAMtools/1.9-foss-2018b

## Load Nextflow module
module load Nextflow/22.04.0


## Run Nextflow RNAseq pipeline
nextflow run nf-core/rnaseq --input samplesheet.csv --outdir PpaDr --fasta /lustre/projects/Research_Project-T110796/nf_10762/el_paco/El_Paco_genome_v3.fa --gtf /lustre/projects/Research_Project-T110796/nf_10762/el_paco/El_Paco_genome_v3.gtf --email rw617@exeter.ac.uk --skip_dupradar true --skip_umi_extract true --remove_ribo_rna true --skip_fastqc true --skip_trimming true --skip_deseq2_qc false --aligner hisat2 --skip_bigwig true --skip_biotype_qc true -profile conda -resume
