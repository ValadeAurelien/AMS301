#!/bin/bash

for np in $(seq $min $step $max)
do
    mpirun -np $np ./parallel $Nx $Ny $nb_steps $t_end $a $b $kappa $CFL_max $T_max $sigma $plot_every $calc_err
done 
