#!/bin/bash
path_conf="../../../tests/confs/su3/Landau_U1/gluodynamics/32^4/beta6.2/conf_Landau_gaugefixed_0001"
conf_format=double_qc2dstag
path_conf_monopole="./result/conf_monopole_32_0001"
path_conf_monopoless="./result/conf_monopoless_32_0001"
path_inverse_laplacian="../../../tests/confs/su2/inverse_laplacian/ALPHA32x32_d.LAT"
bytes_skip=0
x_size=32
y_size=32
z_size=32
t_size=32
parallel=false

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_conf_monopole ${path_conf_monopole} -path_conf_monopoless ${path_conf_monopoless} \
    -path_inverse_laplacian ${path_inverse_laplacian} -bytes_skip ${bytes_skip} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size} \
    -parallel ${parallel}"

../decomposition_su3_test ${parameters}