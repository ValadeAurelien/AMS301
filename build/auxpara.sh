#! /bin/sh

mpirun ./parallel $Nx $Ny $nb_steps $t_end $a $b $kappa $CFL_max $T_max $sigma $plot_every $calc_err 


