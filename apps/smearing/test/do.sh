#!/bin/bash
conf_format_wilson="double_qc2dstag"
conf_path_wilson="../../../tests/confs/SU3_conf/gluodynamics/24^4/beta6.0/CONF0001"
bytes_skip_wilson=0
wilson_type="su3"
conf_format_plaket="double_qc2dstag"
# conf_path_plaket="../../../tests/confs/su2/su2_suzuki/24^4/beta2.4/CON_MC_001.LAT"
conf_path_plaket="../../../tests/confs/SU3_conf/gluodynamics/24^4/beta6.0/CONF0001"
bytes_skip_plaket=0
plaket_type="su3"
HYP_alpha1=1
HYP_alpha2=1
HYP_alpha3=0.5
APE_alpha=0.5
APE_enabled=1
HYP_enabled=1
APE_steps=10
calculation_step_APE=10
calculation_APE_start=10
HYP_steps=1
L_spat=24
L_time=24
path_wilson=./result/wilson_loops
path_flux=./result/flux_tube
wilson_enabled=0
flux_enabled=0
T_min=1
T_max=24
R_min=1
R_max=24
save_conf=1
conf_path_output="/home/ilya/soft/lattice/general_code/apps/smearing/test/result/smeared_0001"
parameters="-conf_format_wilson ${conf_format_wilson} -conf_path_wilson ${conf_path_wilson} -bytes_skip_wilson ${bytes_skip_wilson}\
    -conf_format_plaket ${conf_format_plaket} -conf_path_plaket ${conf_path_plaket} -bytes_skip_plaket ${bytes_skip_plaket}\
    -HYP_alpha1 ${HYP_alpha1} -HYP_alpha2 ${HYP_alpha2} -HYP_alpha3 ${HYP_alpha3}\
    -APE_alpha ${APE_alpha} -APE_enabled ${APE_enabled} -HYP_enabled ${HYP_enabled}\
    -APE_steps ${APE_steps} -HYP_steps ${HYP_steps} -L_spat ${L_spat} -L_time ${L_time} -conf_path_output ${conf_path_output}\
    -path_wilson ${path_wilson} -path_flux ${path_flux} -wilson_enabled ${wilson_enabled} -flux_enabled ${flux_enabled} -save_conf ${save_conf}\
    -T_min ${T_min} -T_max ${T_max} -R_min ${R_min} -R_max ${R_max} -calculation_step_APE ${calculation_step_APE}\
    -calculation_APE_start ${calculation_APE_start}"

../smearing_${wilson_type}_${plaket_type}_test ${parameters}