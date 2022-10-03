#!/bin/bash

################################################################################
#                 Script for generating rRNA alignment job                     #
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
# BAM files will go. This should be nested within the fastq directory.         #
# From inside that directory run the script:                                   #
#                                                                              #
# sh generate-bwa-rRNA-commands-SAM.sh                                         #
#                                                                              #
# This will create a job script containing all the bwa commands ready to       #
# submit to the queue.                                                         #
#                                                                              #
################################################################################


cat <<\EOF > bwa-rRNA-job-script-SAM.sh
#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-T110796 # research project to submit under
#SBATCH --nodes=2 # specify number of nodes
#SBATCH --ntasks-per-node=8 # specify number of processors per node
#SBATCH --mem=30G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=hp508@exeter.ac.uk # email address

# Commands you wish to run must go here, after the SLURM directives]
# load the modules required by the job
module load BWA/0.7.17-foss-2018a
module load SAMtools/1.7-foss-2018a


# absolute path to reference file
ref=/lustre/projects/Research_Project-T110796/Project_10762/el_paco_ref/el_paco_ref/rRNA-reference.fasta
EOF

################################################################################
#
# Loop to generate commands to align reads to the rRNA reference fasta and output in SAM format
#
################################################################################

for fname in /lustre/projects/Research_Project-T110796/Project_10762/ages_alignments/fastqs/*_R1_cuta.fastq.gz
do
    SAMPLE=${fname%_R1*}
    OUTPUT=$(echo $SAMPLE | cut -d '/' -f 3)
        printf "bwa mem /lustre/projects/Research_Project-T110796/Project_10762/el_paco_ref/el_paco_ref/rRNA-reference.fasta  "${SAMPLE}_R1_cuta.fa$
done




