#!/bin/bash
path_conf="../../../tests/confs/su2/gluodynamics/64^4/beta2.9/CONF0002"
conf_format="qcdstag"
file_precision="double"
bytes_skip=0
path_functional_output="result/functional_thermalization_test.csv"
T_step=0.1
T_init=2.5
T_final=0.1
OR_steps=6
thermalization_steps=20
local_thermalization_steps=3
x_size=64
y_size=64
z_size=64
t_size=64

../mag_thermalization_test_test --path_conf ${path_conf} --conf_format ${conf_format} --file_precision ${file_precision} --bytes_skip ${bytes_skip}\
 --path_functional_output ${path_functional_output} --local_thermalization_steps ${local_thermalization_steps}\
 --T_step ${T_step} --T_init ${T_init} --T_final ${T_final} --OR_steps ${OR_steps} --thermalization_steps ${thermalization_steps}\
 --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}
