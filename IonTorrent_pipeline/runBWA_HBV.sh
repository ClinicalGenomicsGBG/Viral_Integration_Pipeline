#!/bin/bash -l

#$ -S /bin/bash
#$ -N BWA
#$ -j y
#$ -M sanna.abrahamsson@gu.se
#$ -m bea
#$ -cwd
#$ -pe mpi 10
#$ -q bfx_short.q@node1-bfx.medair.l
#$ -q bfxcore.q@node4-bfx.medair.lcl


#genomefasta=/jumbo/WorkingDir/B21-005/Data/Meta/db/HBV_D.fasta
genomefasta=/jumbo/WorkingDir/B21-005/Data/Meta/db/HBV_D_HG19.fasta
#humanfasta=/jumbo/db/Homo_sapiens/Ensembl/GRCh38.90/BwaIndex/Homo_sapiens.GRCh38.dna.toplevel.canonical.fa


echo "
module load bwa/0.7.5a
module load samtools/1.6
"

module load bwa/0.7.5a
module load samtools/1.6


echo "
bwa mem -t 10 $genomefasta $1 | samtools sort -@15 -O BAM -o ${1%_Adapterfilt.fastq.gz}.Sorted.bam - 
"

bwa mem -t 10 $genomefasta $1 | samtools sort -@15 -O BAM -o ${1%_Adapterfilt.fastq.gz}.Sorted.bam - 

echo "Finished"
