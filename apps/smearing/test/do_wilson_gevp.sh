#!/bin/bash
conf_format_wilson="ildg"
conf_path_wilson="../../../tests/confs/su3/gluodynamics/16^4/beta6.0/b6p00_L16x16x16x16.01001.lime"
bytes_skip_wilson=0
wilson_type="su3"
convert_wilson=0
HYP_alpha1=1
HYP_alpha2=1
HYP_alpha3=0.5
APE_alpha=0.6
HYP_enabled=1
APE_steps=2
calculation_step_APE=1
calculation_APE_start=1
HYP_steps=0
L_spat=16
L_time=16
path_wilson=./result/wilson_loops_gevp
T_min=1
T_max=4
R_min=1
R_max=4
N_dir=4
parameters="-conf_format_wilson ${conf_format_wilson} -conf_path_wilson ${conf_path_wilson} -bytes_skip_wilson ${bytes_skip_wilson} -convert_wilson ${convert_wilson}\
    -HYP_alpha1 ${HYP_alpha1} -HYP_alpha2 ${HYP_alpha2} -HYP_alpha3 ${HYP_alpha3}\
    -APE_alpha ${APE_alpha} -HYP_enabled ${HYP_enabled} -APE_steps ${APE_steps} -HYP_steps ${HYP_steps}\
    -L_spat ${L_spat} -L_time ${L_time} -path_wilson ${path_wilson} -N_dir ${N_dir}\
    -T_min ${T_min} -T_max ${T_max} -R_min ${R_min} -R_max ${R_max} -calculation_step_APE ${calculation_step_APE}\
    -calculation_APE_start ${calculation_APE_start}"

../smearing_wilson_gevp_${wilson_type}_test ${parameters}