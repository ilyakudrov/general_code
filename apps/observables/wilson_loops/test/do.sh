#!/bin/bash
conf_format="double_qc2dstag"
# conf_path="../../../../tests/confs/su2/gluodynamics/32^3x8/beta2.542/CONF0001"
conf_path="../../../../tests/confs/su3/gluodynamics/32^3x64/beta6.2/s1/CONF0001"
bytes_skip=0
matrix_type="su3"
L_spat=32
L_time=64
# representation="adjoint"
representation="fundamental"
path_wilson=./result/wilson_loops_${representation}
T_min=1
T_max=8
R_min=1
R_max=8
parameters="-conf_format ${conf_format} -conf_path ${conf_path} -bytes_skip ${bytes_skip}\
    -L_spat ${L_spat} -L_time ${L_time} -path_wilson ${path_wilson} -representation ${representation}\
    -T_min ${T_min} -T_max ${T_max} -R_min ${R_min} -R_max ${R_max}"

../wilson_loops_${matrix_type}_test ${parameters}