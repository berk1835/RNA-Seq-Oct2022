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

# Index the reference genome
# usage: hisat2-build [options] reference outfile_prefix
# reference can be comma separated list of chromosome fasta files e.g. chr1.fa,chr2.fa,chr3.fa
hisat2-build chr1.fa,chr2.fa,chr3.fa PpacREF