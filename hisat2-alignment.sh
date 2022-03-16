#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D /lustre/projects/Research_Project-T110796/nys_project/RawReads # set working directory
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=08:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-XXXYYY # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=8 # specify number of processors per node
#SBATCH --mem=30G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=a.n.other@exeter.ac.uk # email address

# Commands you wish to run must go here, after the SLURM directives]
# load the modules required by the job
module load HISAT2/2.0.4-foss-2016b

# script to align RNA-seq to reference genome
# requires:
# HISAT2 version 2.0.4
# citation https://doi.org/10.1038/s41587-019-0201-4
# Samtools version 1.15
# citation https://doi.org/10.1093/bioinformatics/btp352 

# run this script in the folder containing your fastq files
# specify location of reference index
REF=/lustre/projects/Research_Project-T110796/el_paco_ref/El_Paco_Index

# align reads to reference
# loops command for each R1 file in the directory and matches it with the R2 file
for fname in ./nys_fastqs/*_R1_001_fastp.fastq.gz
do
    SAMPLE=${fname%_R1*}
	hisat2 -p 8 -x $REF --rna-strandness RF -1 "${SAMPLE}_R1_001_fastp.fastq.gz" -2 "${SAMPLE}_R2_001_fastp.fastq.gz" -S "${SAMPLE}_unsorted.sam" 2> "${SAMPLE}_unsorted_stats.txt"
done

