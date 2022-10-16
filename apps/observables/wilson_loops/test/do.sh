#!/bin/bash
conf_format="ildg"
conf_path="../../../../tests/confs/Coulomb_su3/QCD/140MeV/nt16/conf_Coulomb_gaugefixed_0501"
bytes_skip=0
matrix_type="su3"
L_spat=64
L_time=16
path_wilson=./result/wilson_loops
T_min=1
T_max=14
R_min=1
R_max=32
parameters="-conf_format ${conf_format} -conf_path ${conf_path} -bytes_skip ${bytes_skip}\
    -L_spat ${L_spat} -L_time ${L_time} -path_wilson ${path_wilson}\
    -T_min ${T_min} -T_max ${T_max} -R_min ${R_min} -R_max ${R_max}"

../wilson_loops_${matrix_type}_test ${parameters}