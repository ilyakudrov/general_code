#!/bin/bash
path_conf="../../../tests/confs/Coulomb_su3/QCD/140MeV/nt4/conf_Coulomb_gaugefixed_0501"
conf_format=ildg
matrix_type="su3"
path_output_correlator="./result/polyakov_correlator_0501"
bytes_skip=0
x_size=64
y_size=64
z_size=64
t_size=4
D_max=32

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_output_correlator ${path_output_correlator}\
    -bytes_skip ${bytes_skip} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size} -D_max ${D_max}"

../polyakov_loop_correlator_${matrix_type}_test ${parameters}