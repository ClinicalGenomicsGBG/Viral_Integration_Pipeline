#!/bin/bash -l

#$ -S /bin/bash
#$ -N Cutadapt
#$ -j y
#$ -M sanna.abrahamsson@gu.se
#$ -m bea
#$ -cwd
#$ -pe mpi 1
#$ -q bfxcore.q@node3-bfx.medair.lcl
#$ -q bfxcore.q@node4-bfx.medair.lcl

# Running Trimgalore, for removing illumina universal Adapter


echo "module load cutadapt/1.9"

module load cutadapt/1.9

commontag=GCCAGGTTCCAGTCAC

echo "cutadapt -g $commontag -o ${1%.fastq.gz}_Adapterfilt.fastq.gz $1"

cutadapt -g $commontag -o ${1%.fastq.gz}_Adapterfilt.fastq.gz $1

echo "Finished"
