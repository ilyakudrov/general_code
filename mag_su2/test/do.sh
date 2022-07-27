#!/bin/bash
path_conf="/home/ilya/soft/lattice/general_code/tests/confs/su2_suzuki/24^4/beta2.4/CON_MC_001.LAT"
conf_format="double"
bytes_skip=8
path_conf_output="/home/ilya/soft/lattice/general_code/tests/confs/MA_gauge/su2_suzuki/conf_gaugefixed/24^4/beta2.4/conf_gaugefixed_0001"
path_spins_output="/home/ilya/soft/lattice/general_code/mag_su2/test/spins"
path_previous="/home/ilya/soft/lattice/general_code/tests/confs/MA_gauge/su2_suzuki/spins/24^4/beta2.4/spins_0001"
path_functional_output="/home/ilya/soft/lattice/general_code/mag_su2/test/functional"
T_step=0.1
T_init=2.5
T_final=0.4
OR_steps=6
thermalization_steps=10
tolerance_maximal=1e-11
tolerance_average=1e-15
tolerance_digits=7
gauge_copies=5
is_new_trial=0
is_final=1
is_compare=1
is_compare_spins=1
is_functional_save=1
x_size=24
y_size=24
z_size=24
t_size=24

/home/ilya/soft/lattice/general_code/mag_su2/mag_fixation_test -path_conf ${path_conf} -conf_format ${conf_format} -bytes_skip ${bytes_skip}\
 -path_spins_output ${path_spins_output} -path_conf_output ${path_conf_output} -path_previous ${path_previous} -path_functional_output ${path_functional_output}\
 -T_step ${T_step} -T_init ${T_init} -T_final ${T_final} -OR_steps ${OR_steps} -thermalization_steps ${thermalization_steps}\
 -tolerance_maximal ${tolerance_maximal} -tolerance_average ${tolerance_average} -tolerance_digits ${tolerance_digits} -gauge_copies ${gauge_copies} -is_new_trial ${is_new_trial}\
 -is_final ${is_final} -is_compare ${is_compare} -is_compare_spins ${is_compare_spins} -is_functional_save ${is_functional_save} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size}
