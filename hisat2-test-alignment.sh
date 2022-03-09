#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D /lustre/projects/Research_Project-T110796/nys_project/RawReads # set working directory
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

# run this script in the folder containing your fastq files
# specify destination directory for results
destdir=/lustre/projects/Research_Project-T110796/nys_project/analysis/alignments
# specify location of reference index
REF=/lustre/projects/Research_Project-T110796/nys_project/Reference

read1=$1
read2=$2
ID="awk '{split($1, array "_"); print array[2]}'"
hisat2 -x REF --rna-strandness RF -1 $1 -2 $2 -S "$destdir/${ID}_alignment.sam"


