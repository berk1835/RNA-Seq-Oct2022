
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

# run this script in the reference directory containing your fasta and annotation files

### gffread to convert annotation
# annotation needs to be in gtf format
# load Cufflinks
module load Cufflinks/2.2.1-foss-2016b 
# convert gff3 to gtf 
gffread El_Paco_V3_gene_annotations.gff3 -T -o El_Paco_V3_gene_annotations.gtf

### HISAT2 indexing
# load the HISAT2 module
module load HISAT2/2.0.4-foss-2016b
# run hisat2 index command with reference fasta and output prefix
hisat2-build El_Paco_genome.fa El_Paco_Index

### STAR indexing
# switch compiler version then load STAR
module switch foss/2016b foss/2018b
module switch zlib/1.2.8-foss-2016b zlib/1.2.11-GCCcore-7.3.0 
module load STAR/2.7.3a-foss-2018b

# make the output directory (genomeDir) first
# run STAR in index mode
# needs the reference fasta and annotation in gtf format
STAR --runThread 8 --runMode genomeGenerate --genomeDir ./STAR-Index \
--genomeFastaFiles El_Paco_genome.fa \
--sjdbGTFfile El_Paco_V3_gene_annotations.gtf 




