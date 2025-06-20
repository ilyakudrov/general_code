#!/bin/bash
path_conf="../../../tests/confs/su3/gluodynamics/24^4/beta6.0/CONF0001"
conf_format="qcdstag"
matrix_type="su3"
file_precision="double"
path_output_correlator="./result/polyakov_correlator_color_average_0001"
bytes_skip=0
x_size=24
y_size=24
z_size=24
t_size=24
D_max=24
correlator_type="color_average"
# correlator_type="singlet"

parameters="--path_conf ${path_conf} --file_precision ${file_precision} --conf_format ${conf_format} --path_output_correlator ${path_output_correlator}\
    --bytes_skip ${bytes_skip} --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size} --D_max ${D_max} --correlator_type ${correlator_type}"

../polyakov_loop_correlator_${matrix_type}_test ${parameters}