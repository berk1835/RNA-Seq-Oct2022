
cat <<\EOF > sortedbamtofastq-job-script.sh
#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory
#SBATCH -p pq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-T110796 # research project to submit under
#SBATCH --nodes=2 # specify number of nodes
#SBATCH --ntasks-per-node=8 # specify number of processors per node
#SBATCH --mem=40G # specify bytes memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=rw617@exeter.ac.uk # email address

# Commands you wish to run must go here, after the SLURM directives
# load the modules required by the job
module load SAMtools/1.7-foss-2018a


ref=/lustre/projects/Research_Project-T110796/Project_10762/el_paco_ref/el_paco_ref/rRNA-reference.fasta
EOF

################################################################################################
#
# Loop to generate commands to convert sorted aligned bam to R1 and R2 fastqs of non-rRNA reads
#
################################################################################################


for fname in /lustre/projects/Research_Project-T110796/Project_10762/ages_alignments/sorted_bams/*_sorted.bam*
do
        SAMPLE=${fname%*_sorted.bam*}
        printf "samtools bam2fq "${SAMPLE}_sorted.bam" > "${SAMPLE}.fastq" | cat "${SAMPLE}.fastq" | grep '/1' -A 3 --no-group-separator > "${SAMPLE}_r1.fastq" | cat "${SAMPLE}.fastq" | grep '/2' -A 3 --no-group-separator > "${SAMPLE}_r2.fastq"\n" >> sortedbamtofastq-job-script.sh

done

