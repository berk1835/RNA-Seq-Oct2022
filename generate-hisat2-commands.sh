#!/bin/bash

set -u
set -e
set -o pipefail

################################################################################
#                 Script for generating HISAT2 aligner commands                #
#                                                                              #
# This script will generate all the commands needed to align paired end fastqs #
# to a reference genome. Script may need changing based on file name changes.  #
#                                                                              #
# Author          : Harry Pollitt                                              #                
# Email           : hap39@aber.ac.uk    hp508@exeter.ac.uk                     #
# GitHub          : https://github.com/Harry-Pollitt                           #               
#                                                                              #                                           
#                                                                              #
# To use, copy or symbolic link this script into the directory where your      #
# BAM files will go. This should be nested within the fastq directory.         #
# From inside that directory run the script:                                   #
#                                                                              #
# sh generate-hisat2-commands.sh                                               #
#                                                                              #
# This will create a job script containing all the hisat2 commands ready to    #
# submit to the queue.                                                         #
#                                                                              #
################################################################################


cat <<\EOF > hisat2-job-script.sh
#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=03:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-T110796 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=8 # specify number of processors per node
#SBATCH --mem=30G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=hp508@exeter.ac.uk # email address

# Commands you wish to run must go here, after the SLURM directives]
# load the modules required by the job
module load HISAT2/2.0.4-foss-2016b
module load SAMtools/1.5-foss-2016b

# absolute path to reference files
ref=/lustre/projects/Research_Project-T110796/el_paco_ref/El_Paco_Index
EOF

################################################################################
# 
# Loop to generate commands to align reads to the reference genome
#
################################################################################

for fname in ../*Ppa*_R1.fastq
do
    sample=${fname%_R1*}
    output=$(echo $sample | cut -d '/' -f 3)
        printf "hisat2 -p 8 -x \$ref --rna-strandness RF -1 "${sample}_R1.fastq" -2 "${sample}_R2.fastq" | samtools sort -@ 8 -T "${output}" -o "${output}_hisat_sorted.bam" - \n" >> hisat2-job-script.sh
done


