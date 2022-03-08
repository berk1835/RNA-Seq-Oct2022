# Title

RNA-seq and differential expression

## Install

Tools used in this project:

HISAT2/2.0.4-foss-2016b

HTSeq/0.9.1-foss-2016b-Python-2.7.12

## Usage

ALL scripts are to be submitted as jobs using "sbatch script.sh"

Check allocation of threads, time and node before submitting

hisat2-index.sh should be used first as alignment requires an indexed reference genome
 
hisat2-alignment.sh aligns all paired fastq sequence reads in a directory to the reference genome indexed previously

htseq-feature-count.sh produces a csv of feature/gene counts based on the hisat2 alignments

## Contributing

Harry Pollitt, Rebekah White

## License

CC BY-ND
