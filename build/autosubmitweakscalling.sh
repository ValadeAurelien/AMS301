#!/bin/bash

mNx=10  #par cellule
mNy=10  #par cellule (garder Ny>Nx)
t_end=1 #garder=1
a=1
b=1
kappa=0.75   #garder <1
CFL_max=0.25
T_max=1      
sigma=0.5
plot_every=0
calc_err=0

#on fait que des carr�s pour être sur du découpage
NPminsqr=2      #coté ducarré minimal
NPmaxsqr=11     #coté max
NPstepsqr=1     #pas 

NPmax=$(echo "$NPmaxsqr^2" | bc)            #nombre de coeurs néc�essais a l'étude
nb_steps=$(echo "4*($NPmax*$mNy)^2" | bc)   #nombre de pas pour etre en dessous de la CFL


OUTPUT_NAME=output_gin_weak_${mNx}_${mNy}_${nb_steps}_${NPminsqr}_${NPstepsqr}_${NPmaxsqr}
ARGS_LIST="-v mNx=$mNx -v mNy=$mNy -v nb_steps=$nb_steps -v t_end=$t_end -v a=$a -v b=$b -v kappa=$kappa -v CFL_max=$CFL_max -v T_max=$T_max -v sigma=$sigma -v plot_every=$plot_every -v calc_err=$calc_err -v min=$NPminsqr -v step=$NPstepsqr -v max=$NPmaxsqr"

qsub -N $OUTPUT_NAME -cwd -pe orte $NPmax $ARGS_LIST auxweakscalling.sh

