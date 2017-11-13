#!/bin/bash

for npsqr in $(seq $min $step $max)
do
    np=$(echo "$npsqr^2" | bc)
    Nx=$(echo "$npsqr*$mNx" | bc)
    Ny=$(echo "$npsqr*$mNy" | bc)
    mpirun -np $np ./parallel $Nx $Ny $nb_steps $t_end $a $b $kappa $CFL_max $T_max $sigma $plot_every $calc_err
done 
