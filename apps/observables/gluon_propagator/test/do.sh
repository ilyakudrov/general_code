#!/bin/bash
conf_format="float"
conf_path="../../../../tests/confs/su2/gluodynamics/66^3x8/beta2.701/CONF_B2701L66T08"
bytes_skip=80
matrix_type="su2"
convert=0
L_spat=66
L_time=8
path=./result/gluon_propagator
parameters="-conf_format ${conf_format} -conf_path ${conf_path} -bytes_skip ${bytes_skip} -convert ${convert}\
    -L_spat ${L_spat} -L_time ${L_time} -path ${path}"

../gluon_propagator_${matrix_type}_test ${parameters}