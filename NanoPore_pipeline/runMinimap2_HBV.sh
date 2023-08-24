#!/bin/bash -l

#$ -S /bin/bash
#$ -N Minimap2
#$ -j y
#$ -m bea
#$ -cwd
#$ -pe mpi 5
#$ -q development.q

module load samtools/1.9

Genome=/medstore/projects/B21-005/Data/Meta/db/HBV_D_HG19.fasta

echo "

~/Programs/minimap2/minimap2 -a -Y -t 5 $Genome $1 | samtools sort -o ${1%.fastq}.bam

"

~/Programs/minimap2/minimap2 -a -Y -t 5 $Genome $1 | samtools sort -o ${1%.fastq}.bam

echo "

samtools index ${1%.fastq}.bam

"

samtools index ${1%.fastq}.bam

echo "Finished"
