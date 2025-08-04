#!/bin/bash
path_conf="/home/ilya/soft/lattice/general_code/apps/mag_su2/test/result/conf_gaugefixed"
conf_format="lexicographical"
file_precision="double"
bytes_skip=0
T_step=0.1
T_init=10
T_final=0.1
OR_steps=4
thermalization_steps=20
tolerance_maximal=1e-10
tolerance_average=1e-14
x_size=64
y_size=64
z_size=64
t_size=64

../Landau_U1_fixation_test --path_conf ${path_conf} --conf_format ${conf_format} --file_precision ${file_precision} --bytes_skip ${bytes_skip}\
 --T_step ${T_step} --T_init ${T_init} --T_final ${T_final} --OR_steps ${OR_steps} --thermalization_steps ${thermalization_steps}\
 --tolerance_maximal ${tolerance_maximal} --tolerance_average ${tolerance_average} --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}
