#!/bin/bash
conf_format="lexicographical"
file_precision="float"
# conf_format="double_qc2dstag"
conf_path="../../../../tests/confs/su2/gluodynamics/66^3x8/beta2.701/CONF_B2701L66T08"
# conf_path="../../../../tests/confs/su2/gluodynamics/32^3x8/beta2.779/CONF0001"
bytes_skip=76
# bytes_skip=0
matrix_type="su2"
beta=2.701
convert=0
L_spat=66
L_time=8
output_path=./result/gluon_propagator
parameters="-conf_format ${conf_format} -conf_path ${conf_path} -bytes_skip ${bytes_skip} -convert ${convert}\
    -L_spat ${L_spat} -L_time ${L_time} -output_path ${output_path} -beta ${beta}"

../gluon_propagator_${matrix_type}_test ${parameters}