#! /bin/sh

Nx=150
Ny=150
nb_steps=100000
t_end=1
a=1
b=1.5
kappa=0.75
CFL_max=0.25
T_max=1
sigma=0.5
plot_every=0
calc_err=1
nb_cores=24

make && mpirun -np $nb_cores ./para $Nx $Ny $nb_steps $t_end $a $b $kappa $CFL_max$ $T_max $sigma $plot_every $calc_err 


