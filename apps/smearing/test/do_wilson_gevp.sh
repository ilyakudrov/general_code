#!/bin/bash
conf_format="qcdstag"
file_precision="double"
conf_path="../../../tests/confs/su2/qc2dstag/40^4/mu0.00/CONF0201"
bytes_skip=0
matrix_type="su2"
convert=0
HYP_alpha1=1
HYP_alpha2=1
HYP_alpha3=0.5
APE_alpha=0.5
HYP_enabled=1
APE_steps=71
calculation_step_APE=20
calculation_APE_start=31
HYP_steps=2
L_spat=40
L_time=40
path_wilson=./result/wilson_loops_gevp
T_min=1
T_max=20
R_min=1
R_max=20
N_dir=1
representation="fundamental"
parameters="--conf_format ${conf_format} --file_precision ${file_precision} --conf_path ${conf_path} --bytes_skip ${bytes_skip} --convert ${convert}\
    --HYP_alpha1 ${HYP_alpha1} --HYP_alpha2 ${HYP_alpha2} --HYP_alpha3 ${HYP_alpha3} --representation ${representation}\
    --APE_alpha ${APE_alpha} --HYP_enabled ${HYP_enabled} --APE_steps ${APE_steps} --HYP_steps ${HYP_steps}\
    --L_spat ${L_spat} --L_time ${L_time} --path_wilson ${path_wilson} --N_dir ${N_dir}\
    --T_min ${T_min} --T_max ${T_max} --R_min ${R_min} --R_max ${R_max} --calculation_step_APE ${calculation_step_APE}\
    --calculation_APE_start ${calculation_APE_start}"

../smearing_wilson_gevp_${matrix_type}_test ${parameters}
