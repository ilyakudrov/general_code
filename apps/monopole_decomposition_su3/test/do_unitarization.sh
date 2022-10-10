#!/bin/bash
path_conf="../../../tests/confs/Landau_U1xU1/gluodynamics/36^4/beta6.3/conf_Landau_gaugefixed_0001"
conf_format=double_qc2dstag
path_conf_monopole="./result/conf_monopole_0001"
path_conf_monopoless="./result/conf_monopoless_36_0001"
bytes_skip=0
x_size=36
y_size=36
z_size=36
t_size=36

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_conf_monopole ${path_conf_monopole} -path_conf_monopoless ${path_conf_monopoless} \
    -bytes_skip ${bytes_skip} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size}"

../monopole_unitarity_test ${parameters}