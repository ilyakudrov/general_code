#!/bin/bash
path_conf="../../../tests/confs/smeared/qc2dstag/40^4/mu0.00/HYP0_alpha=1_1_0.5_APE_alpha=0.5/smeared_0201"
conf_format="double"
matrix_type="su2"
path_output_correlator="./result/polyakov_correlator_0201"
bytes_skip=0
x_size=40
y_size=40
z_size=40
t_size=40
D_max=20
correlator_type="color_average"

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_output_correlator ${path_output_correlator}\
    -bytes_skip ${bytes_skip} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size} -D_max ${D_max} -correlator_type ${correlator_type}"

../polyakov_loop_correlator_${matrix_type}_test ${parameters}