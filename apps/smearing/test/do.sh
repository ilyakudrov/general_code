#!/bin/bash
conf_format_wilson="ildg"
# conf_format_wilson="double"
# conf_format_wilson="double_qc2dstag"
# conf_path_wilson="../../../tests/confs/MAG/su2/qc2dstag/40^4/mu0.00/conf_abelian_0201"
# conf_path_wilson="../../../tests/confs/su2/qc2dstag/40^4/mu0.00/CONF0201"
#conf_path_wilson="../../../tests/confs/su3/QCD/140MeV/nt4/conf.0501"
# conf_path_wilson="../../../tests/confs/su3/gluodynamics/16^4/beta6.0/b6p00_L16x16x16x16.01001.lime"
# conf_path_wilson="../../../tests/confs/MAG/su3/gluodynamics/16^4/beta6.0/steps_0/copies=20/conf_gaugefixed_01001.lime_1"
# conf_path_wilson="../../../tests/confs/monopoless/su3/gluodynamics/16^4/beta6.0/steps_0/copies=20/conf_monopoless_1001_1"
conf_path_wilson="../../../tests/confs/su3/QCD/140MeV/nt20/conf.0501"
bytes_skip_wilson=0
wilson_type="su3"
convert_wilson=0
conf_format_plaket="double_qc2dstag"
conf_path_plaket="../../../tests/confs/su3/gluodynamics/24^4/beta6.0/CONF0001"
# conf_path_plaket="../../../tests/confs/su3/gluodynamics/16^4/beta6.0/b6p00_L16x16x16x16.01001.lime"
bytes_skip_plaket=0
plaket_type="su3"
convert_plaket=0
HYP_alpha1=1
HYP_alpha2=1
HYP_alpha3=0.5
APE_alpha=0.6
APE_enabled=0
HYP_enabled=1
APE_steps=11
calculation_step_APE=10
calculation_APE_start=1
calculation_step_HYP=1
calculation_HYP_start=1
HYP_steps=5
L_spat=64
L_time=20
path_wilson=./result/wilson_loops
path_flux=./result/flux_tube
path_polyakov_correlator=./result/polyakov_correlator
path_polyakov_loop=./result/polyakov_loop
wilson_enabled=0
flux_enabled=0
polyakov_correlator_enabled=0
polyakov_loop_enabled=1
T_min=1
T_max=8
R_min=1
R_max=8
polyakov_correlator_D=31
save_conf=0
conf_path_output="./result/smeared_01001"
parameters="-conf_format_wilson ${conf_format_wilson} -conf_path_wilson ${conf_path_wilson} -bytes_skip_wilson ${bytes_skip_wilson} -convert_wilson ${convert_wilson}\
    -conf_format_plaket ${conf_format_plaket} -conf_path_plaket ${conf_path_plaket} -bytes_skip_plaket ${bytes_skip_plaket} -convert_plaket ${convert_plaket}\
    -HYP_alpha1 ${HYP_alpha1} -HYP_alpha2 ${HYP_alpha2} -HYP_alpha3 ${HYP_alpha3}\
    -APE_alpha ${APE_alpha} -APE_enabled ${APE_enabled} -HYP_enabled ${HYP_enabled}\
    -APE_steps ${APE_steps} -HYP_steps ${HYP_steps} -L_spat ${L_spat} -L_time ${L_time} -conf_path_output ${conf_path_output}\
    -path_wilson ${path_wilson} -path_flux ${path_flux} -wilson_enabled ${wilson_enabled} -flux_enabled ${flux_enabled} -save_conf ${save_conf}\
    -path_polyakov_correlator ${path_polyakov_correlator} -polyakov_correlator_D ${polyakov_correlator_D} -path_polyakov_loop ${path_polyakov_loop}\
    -T_min ${T_min} -T_max ${T_max} -R_min ${R_min} -R_max ${R_max} -calculation_step_APE ${calculation_step_APE}\
    -polyakov_correlator_enabled ${polyakov_correlator_enabled} -polyakov_loop_enabled ${polyakov_loop_enabled}\
    -calculation_APE_start ${calculation_APE_start} -calculation_step_HYP ${calculation_step_HYP} -calculation_HYP_start ${calculation_HYP_start}"

../smearing_${wilson_type}_${plaket_type}_test ${parameters}
