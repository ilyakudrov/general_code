#!/bin/bash
path_conf="../../../tests/confs/Landau_U1xU1/48^3x12/beta6.257/steps_0/conf_Landau_gaugefixed_0263"
conf_format=qcdstag
file_precision=double
path_conf_monopole="./result/conf_monopole"
path_conf_monopoless="./result/conf_monopoless"
path_inverse_laplacian="../../../tests/confs/inverse_laplacian/inverse_laplacian_48x12"
bytes_skip=0
x_size=48
y_size=48
z_size=48
t_size=12

parameters="--path_conf ${path_conf} --file_precision ${file_precision} --conf_format ${conf_format} --path_conf_monopole ${path_conf_monopole} --path_conf_monopoless ${path_conf_monopoless} \
    --path_inverse_laplacian ${path_inverse_laplacian} --bytes_skip ${bytes_skip} --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}"

../decomposition_su3_test ${parameters}