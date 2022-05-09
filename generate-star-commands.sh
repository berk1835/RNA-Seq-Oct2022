#!/bin/bash

################################################################################
#          Script for generating the job script align reads wtih STAR          #
#                                                                              #
# This script will generate all the commands needed to align paired end fastqs #
# to a reference genome. Script may need changing based on file name changes.  #
# STAR is fairly slow, allow ~30 minutes per alignment with 8 threads          #
#                                                                              #
# Author          : Harry Pollitt                                              #
# Email           : hap39@aber.ac.uk    hp508@exeter.ac.uk                     #
# GitHub          : https://github.com/Harry-Pollitt                           #
#                                                                              #
#                                                                              #
# To use, copy or symbolic link this script into the directory where your      #
# filtered fastq files are.                                                    #
# From inside that directory run the script:                                   #
#                                                                              #
# sh generate-star-commands.sh                                                 #
#                                                                              #
# This will create a job script containing all the STAR commands ready to      #
# submit to the queue.                                                         #
#                                                                              #
################################################################################


cat <<\EOF > star-job-script.sh
#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=05:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-T110796 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=8 # specify number of processors per node
#SBATCH --mem=30G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=hp508@exeter.ac.uk # email address

# Commands you wish to run must go here, after the SLURM directives
# load the modules required by the job
module load STAR/2.7.3a-foss-2018b

# absolute path to reference directory
ref=/lustre/projects/Research_Project-T110796/el_paco_ref/STAR_Index/
EOF

################################################################################
# 
# Loop to generate commands to align reads to reference with STAR
#
################################################################################

for fname in ../fastqs/*Ppa*_R1.fastq
do
    sample=${fname%_R1*}
    output=$(echo $sample | cut -d '/' -f 3)_star_
    printf "STAR --runThreadN 8 --genomeDir \$ref --readFilesIn "${sample}_R1.fastq" "${sample}_R2.fastq" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $output\n" >> star-job-script.sh

done

