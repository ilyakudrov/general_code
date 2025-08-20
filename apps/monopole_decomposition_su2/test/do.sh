#!/bin/bash
# path_conf="../../../tests/confs/MA_gauge/su2/su2_suzuki/24^4/beta2.6/T_step=0.0001/conf_0001"
path_conf="../../../tests/confs/MA_gauge/su2/gluodynamics/32^3x8/beta2.542/T_step=0.001/s0/conf_0001"
# path_conf="../../../tests/confs/su2/gluodynamics/32^3x8/beta2.542/s0/CONF0001"
conf_format="lexicographical"
file_precision=double
path_inverse_laplacian="../../../tests/confs/inverse_laplacian/ALPHA32x8_d.LAT"
path_conf_monopole="./result/conf_monopole_0001"
path_conf_monopoless="./result/conf_monopoless_0001"
bytes_skip=0
x_size=32
y_size=32
z_size=32
t_size=8

parameters="--conf_format ${conf_format} --path_conf $path_conf --bytes_skip ${bytes_skip} --file_precision ${file_precision}\
    --path_conf_monopole ${path_conf_monopole} --path_conf_monopoless ${path_conf_monopoless} --path_inverse_laplacian ${path_inverse_laplacian}\
    --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}"

../decomposition_test ${parameters}