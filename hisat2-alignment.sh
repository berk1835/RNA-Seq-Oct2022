#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=HH:MM:SS # maximum walltime for the job
#SBATCH -A Research_Project-XXXYYY # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
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

# align reads to reference
destdir=/path/to/destination

for fname in *_R1.fastq.gz
do
    SAMPLE=${fname%_R1*}
	hisat2 -x PpacREF -1 "${SAMPLE}_R1.fastq.gz" -2 "${SAMPLE}_R2.fastq.gz" -S "$destdir/${SAMPLE}_alignment.sam"
done

