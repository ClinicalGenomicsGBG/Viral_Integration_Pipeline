#!/bin/bash -l

#$ -S /bin/bash
#$ -N VirusPipeline_Nanopore
#$ -j y
#$ -m bea
#$ -cwd
#$ -pe mpi 1
#$ -q development.q

echo "

module load miniconda/4.14.0
source activate py2

module load samtools/1.9

"

module load miniconda/4.14.0
source activate py2

module load samtools/1.9

echo "

python /medstore/projects/P23-044/Code/Viral_Integration_Pipeline/NanoPore_pipeline/ScreeningFusionsForSE_Pipeline_vs4.py -b $1 -g HBV_D

"


python /medstore/projects/P23-044/Code/Viral_Integration_Pipeline/NanoPore_pipeline/ScreeningFusionsForSE_Pipeline_vs4.py -b $1 -g HBV_D


echo "

samtools view $1 -h | grep -e $'\tSA:' -e ^@ | grep -e HBV_D -e ^@ | samtools view -Sb - > ${1%.bam}_Chimeric.bam

"

samtools view $1 -h | grep -e $'\tSA:' -e ^@ | grep -e HBV_D -e ^@ | samtools view -Sb - > ${1%.bam}_Chimeric.bam

samtools index ${1%.bam}_Chimeric.bam

echo "Fin"
