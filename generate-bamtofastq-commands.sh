#!/bin/bash

################################################################################
#       Script for generating the job script to extract non-rRNA reads         #
#                                                                              #
# This script will generate all the commands needed to align paired end fastqs #
# to a reference rRNA fasta. Script may need changing based on file            #
# name changes.                                                                #
#                                                                              #
# Author          : Harry Pollitt                                              #
# Email           : hap39@aber.ac.uk    hp508@exeter.ac.uk                     #
# GitHub          : https://github.com/Harry-Pollitt                           #
#                                                                              #
#                                                                              #
# To use, copy or symbolic link this script into the directory where your      #
# rRNA aligned BAM files are. Fastq files will be stored in a sub-directory.   #
# From inside that directory run the script:                                   #
#                                                                              #
# sh generate-bamtofastq-commands.sh                                           #
#                                                                              #
# This will create a job script containing all the SAMTools and BEDTools       #
# commands ready to submit to the queue.                                       #
#                                                                              #
################################################################################


cat <<\EOF > bamtofastq-job-script.sh
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
module load BEDTools/2.27.1-foss-2016b 
module load SAMtools/1.5-foss-2016b
EOF

###################################################################################
# 
# Loop to generate commands to convert rRNA aligned bam to fastqs of non-rRNA reads
#
###################################################################################

# samtools view options
# -f 4 # unmapped reads
# -F 4 # mapped reads

for fname in ./*Ppa*_rRNA_alignment.bam
do
    SAMPLE=${fname%_rRNA*}
        printf "samtools view -@ 8 -b -f 4 "${SAMPLE}_rRNA_alignment.bam" | bedtools bamtofastq -i - -fq "./filtered_fastqs/${SAMPLE}_filtered_R1.fastq" -fq2 "./filtered_fastqs/${SAMPLE}_filtered_R2.fastq" \n" >> bamtofastq-job-script.sh
done
