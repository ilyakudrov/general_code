#!/bin/bash
path_conf="../../../tests/confs/su3/mag/gluodynamics/24^4/beta6.0/CONFDP_gaugefixed_0001"
conf_format=double_qc2dstag
path_conf_monopole="./result/conf_monopole_24_0001"
path_conf_monopoless="./result/conf_monopoless_24_0001"
path_inverse_laplacian="../../../tests/confs/su2/inverse_laplacian/ALPHA24x24_d.LAT"
bytes_skip=0
x_size=24
y_size=24
z_size=24
t_size=24
parallel=true

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_conf_monopole ${path_conf_monopole} -path_conf_monopoless ${path_conf_monopoless} \
    -path_inverse_laplacian ${path_inverse_laplacian} -bytes_skip ${bytes_skip} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size} \
    -parallel ${parallel}"

../decomposition_su3_test ${parameters}