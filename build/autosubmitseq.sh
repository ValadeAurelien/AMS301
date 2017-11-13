#!/bin/bash

Nx=225
Ny=225
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

OUTPUT_NAME=output_gin_seq_${Nx}_${Ny}_${nb_steps}
ARGS_LIST="-v Nx=$Nx -v Ny=$Ny -v nb_steps=$nb_steps -v t_end=$t_end -v a=$a -v b=$b -v kappa=$kappa -v CFL_max=$CFL_max -v T_max=$T_max -v sigma=$sigma -v plot_every=$plot_every -v calc_err=$calc_err"

qsub -N $OUTPUT_NAME -cwd -j y $ARGS_LIST auxseq.sh
