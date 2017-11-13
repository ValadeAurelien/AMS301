#! /bin/sh

Nx=50
Ny=50
nb_steps=10000
t_end=1
a=1
b=1.5
kappa=0.75
CFL_max=0.25
T_max=1
sigma=0.5
plot_every=1
calc_err=1
nb_cores=16

mpirun -np $nb_cores ./parallel $Nx $Ny $nb_steps $t_end $a $b $kappa $CFL_max $T_max $sigma $plot_every $calc_err 


