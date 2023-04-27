#!/bin/bash -l

#$ -S /bin/bash
#$ -N VirusPipeline_Nanopore
#$ -j y
#$ -m bea
#$ -cwd
#$ -pe mpi 1
#$ -q bfxcore.q@node2-bfx.medair.lcl
##$ -q bfxcore.q@node6-bfx.medair.lcl
#$ -q bfxcore.q@node3-bfx.medair.lcl

echo "
module load anaconda2
"

module load anaconda2
source activate Sanna

module load samtools/1.9

echo "

python /jumbo/WorkingDir/B21-005/Code/NanoPore_pipeline/ScreeningFusionsForSE_Pipeline_vs4.py -b $1 -g HBV_D

"

python /jumbo/WorkingDir/B21-005/Code/NanoPore_pipeline/ScreeningFusionsForNanoPore_Pipeline_vs4.py -b $1 -g HBV_D


echo "

samtools view $1 -h | grep -e $'\tSA:' -e ^@ | grep -e HBV_D -e ^@ | samtools view -Sb - > ${1%.bam}_Chimeric.bam

"

samtools view $1 -h | grep -e $'\tSA:' -e ^@ | grep -e HBV_D -e ^@ | samtools view -Sb - > ${1%.bam}_Chimeric.bam



echo "Fin"
