#!/bin/bash
conf_format="qcdstag"
file_precision="double"
conf_path="../../../tests/confs/su3/gluodynamics/24^4/beta6.0/CONF0175"
bytes_skip=0
matrix_type="su3"
convert=0
HYP_alpha1=1
HYP_alpha2=1
HYP_alpha3=0.5
APE_alpha=0.5
APE_enabled=1
HYP_enabled=0
APE_steps=11
calculation_step_APE=10
calculation_APE_start=1
calculation_step_HYP=1
calculation_HYP_start=1
HYP_steps=1
L_spat=24
L_time=24
path_wilson=./result/wilson_loops
path_polyakov_correlator=./result/polyakov_correlator
path_polyakov_loop=./result/polyakov_loop
wilson_enabled=1
polyakov_correlator_enabled=0
polyakov_loop_enabled=0
correlator_type="color_average"
T_min=1
T_max=12
R_min=1
R_max=12
polyakov_correlator_D=31
save_conf=0
conf_path_output="./result/smeared_01001"
parameters="--conf_format ${conf_format} --file_precision ${file_precision} --conf_path ${conf_path} --bytes_skip ${bytes_skip} --convert ${convert}\
    --conf_format_plaket ${conf_format_plaket} --conf_path_plaket ${conf_path_plaket} --bytes_skip_plaket ${bytes_skip_plaket} --convert_plaket ${convert_plaket}\
    --HYP_alpha1 ${HYP_alpha1} --HYP_alpha2 ${HYP_alpha2} --HYP_alpha3 ${HYP_alpha3}\
    --APE_alpha ${APE_alpha} --APE_enabled ${APE_enabled} --HYP_enabled ${HYP_enabled}\
    --APE_steps ${APE_steps} --HYP_steps ${HYP_steps} --L_spat ${L_spat} --L_time ${L_time} --conf_path_output ${conf_path_output}\
    --path_wilson ${path_wilson} --path_flux ${path_flux} --wilson_enabled ${wilson_enabled} --flux_enabled ${flux_enabled} --save_conf ${save_conf}\
    --path_polyakov_correlator ${path_polyakov_correlator} --polyakov_correlator_D ${polyakov_correlator_D} --path_polyakov_loop ${path_polyakov_loop}\
    --T_min ${T_min} --T_max ${T_max} --R_min ${R_min} --R_max ${R_max} --calculation_step_APE ${calculation_step_APE}\
    --polyakov_correlator_enabled ${polyakov_correlator_enabled} --polyakov_loop_enabled ${polyakov_loop_enabled} --correlator_type ${correlator_type}\
    --calculation_APE_start ${calculation_APE_start} --calculation_step_HYP ${calculation_step_HYP} --calculation_HYP_start ${calculation_HYP_start}"

../smearing_${matrix_type}_test ${parameters}
