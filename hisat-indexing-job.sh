#!/bin/bash

#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=00:30:00 # maximum walltime for the job
#SBATCH -A Research_Project-T110796 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=8 # specify number of processors per node
#SBATCH --mem=16G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=hp508@exeter.ac.uk # email address

### Run this script in the HISAT_Index directory within el_paco_ref/

### HISAT2 indexing
# load the HISAT2 module
module load HISAT2/2.0.4-foss-2016b
# run hisat2 index command with reference fasta and output prefix
hisat2-build ../El_Paco_genome.fa El_Paco_Index