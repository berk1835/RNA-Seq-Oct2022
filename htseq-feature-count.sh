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
module load HTSeq/0.9.1-foss-2016b-Python-2.7.12

# HTSeq - feature counting
# citation: https://doi.org/10.1093/bioinformatics/btu638

# make sure chromosome names match in SAM and GTF files.
#python -m HTSeq.scripts.count [options] <alignment_files> <gtf_file> 

for fname in *.sam
do
    SAMPLE=${fname%_alignment.sam}
	python -m HTSeq.scripts.count -c "${SAMPLE}_counts.csv" -t gene -i gene_id "${SAMPLE}_alignment.sam" ref_annotation.gtf
done


