#!/bin/bash
#PBS -N linear-tau=1.00e4-Le=3.0
#PBS -l nodes=1:ppn=24
#PBS -l walltime=30:00:00
#PBS -m bea -M luis.silva@glasgow.ac.uk
#

cd $PBS_O_WORKDIR

pbsdsh -v $HOME/bin/glo_wrap.sh
