#!/bin/bash
#PBS -l walltime=1000:00:00
cd $PBS_O_WORKDIR

#echo $PBS_NODEFILE :
#cat $PBS_NODEFILE

source ~/.bashrc

make exp2a
