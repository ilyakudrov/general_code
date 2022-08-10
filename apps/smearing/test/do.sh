#!/bin/bash
conf_format_wilson=double
conf_path_wilson="../../../tests/confs/su2/su2_suzuki/48^4/beta2.8/CON_fxd_MAG_035.LAT"
bytes_skip_wilson=4
conf_format_plaket=double
conf_path_plaket="../../../tests/confs/su2/su2_suzuki/48^4/beta2.8/CON_fxd_MAG_035.LAT"
bytes_skip_plaket=4
HYP_alpha1=0.75
HYP_alpha2=0.6
HYP_alpha3=0.3
APE_alpha=0.5
APE_enabled=1
HYP_enabled=1
APE_steps=16
HYP_steps=1
L_spat=48
L_time=48
path_wilson=./result/wilson_loops
path_flux=./result/flux_tube
wilson_enabled=1
flux_enabled=1
T_min=1
T_max=12
R_min=1
R_max=12
calculation_step_APE=2
calculation_APE_start=6
parameters="-conf_format_wilson ${conf_format_wilson} -conf_path_wilson ${conf_path_wilson} -bytes_skip_wilson ${bytes_skip_wilson}\
    -conf_format_plaket ${conf_format_plaket} -conf_path_plaket ${conf_path_plaket} -bytes_skip_plaket ${bytes_skip_plaket}\
    -HYP_alpha1 ${HYP_alpha1} -HYP_alpha2 ${HYP_alpha2} -HYP_alpha3 ${HYP_alpha3}\
    -APE_alpha ${APE_alpha} -APE_enabled ${APE_enabled} -HYP_enabled ${HYP_enabled}\
    -APE_steps ${APE_steps} -HYP_steps ${HYP_steps} -L_spat ${L_spat} -L_time ${L_time}\
    -path_wilson ${path_wilson} -path_flux ${path_flux} -wilson_enabled ${wilson_enabled} -flux_enabled ${flux_enabled}\
    -T_min ${T_min} -T_max ${T_max} -R_min ${R_min} -R_max ${R_max} -calculation_step_APE ${calculation_step_APE}\
    -calculation_APE_start ${calculation_APE_start}"

../smearing_su2_su2_test ${parameters}