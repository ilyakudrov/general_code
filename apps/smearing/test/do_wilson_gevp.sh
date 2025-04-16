#!/bin/bash
conf_format_wilson="double_qc2dstag"
conf_path_wilson="../../../tests/confs/MAG/su3/gluodynamics/36^4/beta6.1/HYP0_alpha=1_1_0.5_APE_alpha=0.5/steps_0/copies=10/s1/conf_gaugefixed_0001_0"
bytes_skip_wilson=0
wilson_type="su3"
convert_wilson=0
HYP_alpha1=1
HYP_alpha2=1
HYP_alpha3=0.5
APE_alpha=0.5
HYP_enabled=1
APE_steps=2
calculation_step_APE=1
calculation_APE_start=1
HYP_steps=0
L_spat=36
L_time=36
path_wilson=./result/wilson_loops_gevp
T_min=1
T_max=18
R_min=1
R_max=18
N_dir=1
representation="adjoint"
parameters="-conf_format_wilson ${conf_format_wilson} -conf_path_wilson ${conf_path_wilson} -bytes_skip_wilson ${bytes_skip_wilson} -convert_wilson ${convert_wilson}\
    -HYP_alpha1 ${HYP_alpha1} -HYP_alpha2 ${HYP_alpha2} -HYP_alpha3 ${HYP_alpha3} -representation ${representation}\
    -APE_alpha ${APE_alpha} -HYP_enabled ${HYP_enabled} -APE_steps ${APE_steps} -HYP_steps ${HYP_steps}\
    -L_spat ${L_spat} -L_time ${L_time} -path_wilson ${path_wilson} -N_dir ${N_dir}\
    -T_min ${T_min} -T_max ${T_max} -R_min ${R_min} -R_max ${R_max} -calculation_step_APE ${calculation_step_APE}\
    -calculation_APE_start ${calculation_APE_start}"

../smearing_wilson_gevp_${wilson_type}_test ${parameters}
