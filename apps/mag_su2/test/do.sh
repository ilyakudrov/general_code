#!/bin/bash
path_conf="../../../tests/confs/su2/su2_suzuki/48^4/beta2.7/CON_MC_001.LAT"
conf_format="double"
bytes_skip=8
path_conf_output="result/conf_gaugefixed_0001"
path_spins_output="result/spins_0001"
path_previous="../../../tests/confs/spins/su2/su2_suzuki/48^4/beta2.7/T_step=0.0001/spins_0001"
path_functional_output="result/functional_0001"
T_step=0.01
T_init=2.5
T_final=0.5
OR_steps=4
thermalization_steps=50
tolerance_maximal=1e-8
tolerance_average=1e-12
tolerance_digits=7
gauge_copies=2
# is_new_trial=1
# is_final=0
# is_compare=1
# is_compare_spins=1
# is_functional_save=1
is_new_trial=0
is_final=1
is_compare=0
is_compare_spins=0
is_functional_save=0
x_size=48
y_size=48
z_size=48
t_size=48

../mag_fixation_test -path_conf ${path_conf} -conf_format ${conf_format} -bytes_skip ${bytes_skip}\
 -path_spins_output ${path_spins_output} -path_conf_output ${path_conf_output} -path_previous ${path_previous} -path_functional_output ${path_functional_output}\
 -T_step ${T_step} -T_init ${T_init} -T_final ${T_final} -OR_steps ${OR_steps} -thermalization_steps ${thermalization_steps}\
 -tolerance_maximal ${tolerance_maximal} -tolerance_average ${tolerance_average} -tolerance_digits ${tolerance_digits} -gauge_copies ${gauge_copies} -is_new_trial ${is_new_trial}\
 -is_final ${is_final} -is_compare ${is_compare} -is_compare_spins ${is_compare_spins} -is_functional_save ${is_functional_save} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size}
