#!/bin/bash

Nx=225
Ny=220
nb_steps=200000
t_end=1
a=1
b=1.5
kappa=0.75
CFL_max=0.25
T_max=1
sigma=0.5
plot_every=0
calc_err=1
#
#  Name of the job (used to build the name of the standard output stream)
#$ -N seq_gin
#
#  The job is located in the current working directory
#$ -cwd
#
#  Merge standard error and standard output streams
#$ -j y
#
./seq $Nx $Ny $nb_steps $t_end $a $b $kappa $CFL_max $T_max $sigma $plot_every $calc_err

