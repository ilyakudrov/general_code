#!/bin/bash
path_conf="../../../tests/confs/Landau_U1xU1/gluodynamics/36^4/beta6.3/conf_Landau_gaugefixed_0001"
conf_format=double_qc2dstag
path_conf_monopole="./result/conf_monopole_36_0001"
path_conf_monopoless="./result/conf_monopoless_36_0001"
path_inverse_laplacian="../../../tests/confs/inverse_laplacian/ALPHA36x36_d.LAT"
bytes_skip=0
x_size=36
y_size=36
z_size=36
t_size=36
parallel=true

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_conf_monopole ${path_conf_monopole} -path_conf_monopoless ${path_conf_monopoless} \
    -path_inverse_laplacian ${path_inverse_laplacian} -bytes_skip ${bytes_skip} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size} \
    -parallel ${parallel}"

../decomposition_su3_test ${parameters}