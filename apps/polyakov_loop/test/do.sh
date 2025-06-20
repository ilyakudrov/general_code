#!/bin/bash
path_conf="../../../tests/confs/Coulomb_su3/QCD/140MeV/nt4/conf_Coulomb_gaugefixed_0501"
file_precision=double
conf_format=ildg
matrix_type="su3"
path_output="./result/polyakov_0501"
bytes_skip=0
x_size=64
y_size=64
z_size=64
t_size=4

parameters="--path_conf ${path_conf} --file_precision ${file_precision} --conf_format ${conf_format} --path_output ${path_output}\
    --bytes_skip ${bytes_skip} --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}"

../polyakov_loop_${matrix_type}_test ${parameters}