#!/bin/bash -l

#$ -N Porechops
#$ -j y
#$ -m bea
#$ -cwd
#$ -pe mpi 1
#$ -q development.q

module load miniconda/4.14.0

source activate py3

echo "

/home/xabras/.conda/envs/py3/bin/porechop -i $1 -b $1_PorechopsOut

"

/home/xabras/.conda/envs/py3/bin/porechop -i $1 -b $1_PorechopsOut


echo "Fin"
