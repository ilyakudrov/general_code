#!/bin/bash
t_size=32
n_eig=60
for mu in 3 5
do
#mu=3
w=$(($mu/10))
q=$(($mu - $w*10))
for((i = 1;i <= 2;i++))
do
a=$(($i/1000))
b=$((($i-$a*1000)/100))
c=$((($i-$a*1000-$b*100)/10))
d=$(($i-$a*1000-$b*100-$c*10))
scp kudrov@rrcmpi-a.itep.ru:/home/clusters/rrcmpi/kudrov/eigenvalues/Cuda-Arnoldi/result/time_$t_size/mu0.$w$q/eigenvalues_neig="$n_eig"_$a$b$c$d /home/ilya/lattice/general_code/tests/eigenvectors/time_$t_size/mu0.$w$q
scp kudrov@rrcmpi-a.itep.ru:/home/clusters/rrcmpi/kudrov/eigenvalues/Cuda-Arnoldi/result/time_$t_size/mu0.$w$q/eigenvectors_neig="$n_eig"_$a$b$c$d /home/ilya/lattice/general_code/tests/eigenvectors/time_$t_size/mu0.$w$q
done
done
