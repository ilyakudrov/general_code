#!/bin/bash
conf_format="double"
conf_path="../../../../tests/confs/smeared/qc2dstag/40^4/mu0.00/HYP0_alpha=1_1_0.5_APE_alpha=0.5/smeared_0201"
bytes_skip=0
matrix_type="su2"
L_spat=40
L_time=40
representation="adjoint"
# representation="fundamental"
path_wilson=./result/wilson_loops_${representation}
T_min=1
T_max=4
R_min=1
R_max=4
parameters="-conf_format ${conf_format} -conf_path ${conf_path} -bytes_skip ${bytes_skip}\
    -L_spat ${L_spat} -L_time ${L_time} -path_wilson ${path_wilson} -representation ${representation}\
    -T_min ${T_min} -T_max ${T_max} -R_min ${R_min} -R_max ${R_max}"

../wilson_loops_${matrix_type}_test ${parameters}