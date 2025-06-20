#!/bin/bash
conf_format="qcdstag"
file_precision="double"
conf_path="../../../../tests/confs/su3/gluodynamics/24^4/beta6.0/CONF0001"
# conf_path="../../../../tests/confs/su3/QCD/140MeV/nt20/conf.0501"
bytes_skip=0
matrix_type="su3"
L_spat=24
L_time=24
# representation="adjoint"
representation="fundamental"
path_wilson=./result/wilson_loops_${representation}
T_min=1
T_max=12
R_min=1
R_max=12
axis="off-axis"
# axis="on-axis"
parameters="--conf_format ${conf_format} --file_precision ${file_precision} --conf_path ${conf_path} --bytes_skip ${bytes_skip}\
    --L_spat ${L_spat} --L_time ${L_time} --path_wilson ${path_wilson} --representation ${representation}\
    --axis ${axis} --T_min ${T_min} --T_max ${T_max} --R_min ${R_min} --R_max ${R_max}"

../wilson_loops_${matrix_type}_test ${parameters}