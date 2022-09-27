For ```sh generate-bwa-rRNA-commands.sh```, you would expect to see:


```bash
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
#SBATCH --mail-user=rw617@exeter.ac.uk # email address

# Commands you wish to run must go here, after the SLURM directives]
# load the modules required by the job
module load BWA/0.7.17-foss-2018a
module load SAMtools/1.7-foss-2018a

# absolute path to reference file
ref=/lustre/projects/Research_Project-T110796/el_paco_ref/rRNA-reference.fasta
bwa mem $ref ../../fastqs/10762_Ppa119_S32_R1_cuta.fastq.gz ../../fastqs/10762_Ppa119_S32_R2_cuta.fastq.gz | samtools view -b fastqs_rRNA_alignment.bam
bwa mem $ref ../../fastqs/10762_Ppa122_S33_R1_cuta.fastq.gz ../../fastqs/10762_Ppa122_S33_R2_cuta.fastq.gz | samtools view -b fastqs_rRNA_alignment.bam
bwa mem $ref ../../fastqs/10762_Ppa124_S34_R1_cuta.fastq.gz ../../fastqs/10762_Ppa124_S34_R2_cuta.fastq.gz | samtools view -b fastqs_rRNA_alignment.bam
bwa mem $ref ../../fastqs/10762_Ppa125_S35_R1_cuta.fastq.gz ../../fastqs/10762_Ppa125_S35_R2_cuta.fastq.gz | samtools view -b fastqs_rRN$
bwa mem $ref ../../fastqs/10762_Ppa128_S36_R1_cuta.fastq.gz ../../fastqs/10762_Ppa128_S36_R2_cuta.fastq.gz | samtools view -b fastqs_rRN$
bwa mem $ref ../../fastqs/10762_Ppa129_S37_R1_cuta.fastq.gz ../../fastqs/10762_Ppa129_S37_R2_cuta.fastq.gz | samtools view -b fastqs_rRN$
bwa mem $ref ../../fastqs/10762_Ppa132_S38_R1_cuta.fastq.gz ../../fastqs/10762_Ppa132_S38_R2_cuta.fastq.gz | samtools view -b fastqs_rRN$
bwa mem $ref ../../fastqs/10762_Ppa133_S39_R1_cuta.fastq.gz ../../fastqs/10762_Ppa133_S39_R2_cuta.fastq.gz | samtools view -b fastqs_rRN$
bwa mem $ref ../../fastqs/10762_Ppa140_S40_R1_cuta.fastq.gz ../../fastqs/10762_Ppa140_S40_R2_cuta.fastq.gz | samtools view -b fastqs_rRN$
bwa mem $ref ../../fastqs/10762_Ppa142_S41_R1_cuta.fastq.gz ../../fastqs/10762_Ppa142_S41_R2_cuta.fastq.gz | samtools view -b fastqs_rRN$
bwa mem $ref ../../fastqs/10762_Ppa164_S30_R1_cuta.fastq.gz ../../fastqs/10762_Ppa164_S30_R2_cuta.fastq.gz | samtools view -b fastqs_rRN$
bwa mem $ref ../../fastqs/10762_Ppa168_S31_R1_cuta.fastq.gz ../../fastqs/10762_Ppa168_S31_R2_cuta.fastq.gz | samtools view -b fastqs_rRN$
bwa mem $ref ../../fastqs/10762_Ppa51_S24_R1_cuta.fastq.gz ../../fastqs/10762_Ppa51_S24_R2_cuta.fastq.gz | samtools view -b fastqs_rRNA_$
bwa mem $ref ../../fastqs/10762_Ppa52_S25_R1_cuta.fastq.gz ../../fastqs/10762_Ppa52_S25_R2_cuta.fastq.gz | samtools view -b fastqs_rRNA_$
bwa mem $ref ../../fastqs/10762_Ppa53_S26_R1_cuta.fastq.gz ../../fastqs/10762_Ppa53_S26_R2_cuta.fastq.gz | samtools view -b fastqs_rRNA_$
bwa mem $ref ../../fastqs/10762_Ppa54_S27_R1_cuta.fastq.gz ../../fastqs/10762_Ppa54_S27_R2_cuta.fastq.gz | samtools view -b fastqs_rRNA_$
bwa mem $ref ../../fastqs/10762_Ppa55_S28_R1_cuta.fastq.gz ../../fastqs/10762_Ppa55_S28_R2_cuta.fastq.gz | samtools view -b fastqs_rRNA_$
bwa mem $ref ../../fastqs/10762_Ppa56_S29_R1_cuta.fastq.gz ../../fastqs/10762_Ppa56_S29_R2_cuta.fastq.gz | samtools view -b fastqs_rRNA_$
bwa mem $ref ../../fastqs/10762_PpaAge10_S22_R1_cuta.fastq.gz ../../fastqs/10762_PpaAge10_S22_R2_cuta.fastq.gz | samtools view -b fastqs$
bwa mem $ref ../../fastqs/10762_PpaAge11_S20_R1_cuta.fastq.gz ../../fastqs/10762_PpaAge11_S20_R2_cuta.fastq.gz | samtools view -b fastqs$
bwa mem $ref ../../fastqs/10762_PpaAge12_S23_R1_cuta.fastq.gz ../../fastqs/10762_PpaAge12_S23_R2_cuta.fastq.gz | samtools view -b fastqs$
bwa mem $ref ../../fastqs/10762_PpaAge1_S14_R1_cuta.fastq.gz ../../fastqs/10762_PpaAge1_S14_R2_cuta.fastq.gz | samtools view -b fastqs_r$
bwa mem $ref ../../fastqs/10762_PpaAge2_S15_R1_cuta.fastq.gz ../../fastqs/10762_PpaAge2_S15_R2_cuta.fastq.gz | samtools view -b fastqs_r$
bwa mem $ref ../../fastqs/10762_PpaAge3_S43_R1_cuta.fastq.gz ../../fastqs/10762_PpaAge3_S43_R2_cuta.fastq.gz | samtools view -b fastqs_r$
bwa mem $ref ../../fastqs/10762_PpaAge4_S42_R1_cuta.fastq.gz ../../fastqs/10762_PpaAge4_S42_R2_cuta.fastq.gz | samtools view -b fastqs_r$
bwa mem $ref ../../fastqs/10762_PpaAge5_S16_R1_cuta.fastq.gz ../../fastqs/10762_PpaAge5_S16_R2_cuta.fastq.gz | samtools view -b fastqs_r$
bwa mem $ref ../../fastqs/10762_PpaAge6_S17_R1_cuta.fastq.gz ../../fastqs/10762_PpaAge6_S17_R2_cuta.fastq.gz | samtools view -b fastqs_r$
bwa mem $ref ../../fastqs/10762_PpaAge7_S21_R1_cuta.fastq.gz ../../fastqs/10762_PpaAge7_S21_R2_cuta.fastq.gz | samtools view -b fastqs_r$
 ```
