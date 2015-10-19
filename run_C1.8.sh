#! /bin/bash
#
#$ -S /bin/bash -cwd
#$ -j y
#$ -o nohup.sge
#$ -V
#$ -M jb3816@cumc.columbia.edu
#$ -m bea
# 

qmake -cwd -v PATH -- -j 8
