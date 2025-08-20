#!/bin/bash
path_conf="../../../tests/confs/su2/qc2dstag/40^4/mu0.00/CONF0201"
conf_format="qcdstag"
file_precision="double"
bytes_skip=0
T_step=0.01
T_init=5
T_final=0.01
OR_steps=4
thermalization_steps=20
tolerance_maximal=1e-14
tolerance_average=1e-15
x_size=40
y_size=40
z_size=40
t_size=40

../Landau_su2_fixation_test --path_conf ${path_conf} --conf_format ${conf_format} --file_precision ${file_precision} --bytes_skip ${bytes_skip}\
 --T_step ${T_step} --T_init ${T_init} --T_final ${T_final} --OR_steps ${OR_steps} --thermalization_steps ${thermalization_steps}\
 --tolerance_maximal ${tolerance_maximal} --tolerance_average ${tolerance_average} --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}
