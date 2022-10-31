#!/bin/bash
path_conf="../../../tests/confs/SU3_conf/gluodynamics/16^4/su3_mag_u1.01001.lat"
conf_format=double_vitaly
path_conf_monopole="./result/conf_monopole_16_1001"
path_conf_monopoless="./result/conf_monopoless_16_1001"
path_inverse_laplacian="../../../tests/confs/inverse_laplacian/ALPHA16x16_d.LAT"
bytes_skip=4
x_size=16
y_size=16
z_size=16
t_size=16
parallel=true

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_conf_monopole ${path_conf_monopole} -path_conf_monopoless ${path_conf_monopoless} \
    -path_inverse_laplacian ${path_inverse_laplacian} -bytes_skip ${bytes_skip} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size} \
    -parallel ${parallel}"

../decomposition_su3_test ${parameters}