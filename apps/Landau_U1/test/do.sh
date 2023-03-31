#!/bin/bash
path_conf="../../../tests/confs/su2/su2_suzuki/24^4/beta2.5/conf_0001"
conf_format="double"
bytes_skip=0
path_conf_output="./result/conf_gaugefixed_0001"
path_functional_output="./result/functional_0001"
x_size=24
y_size=24
z_size=24
t_size=24

../Landau_U1_fixation_test -path_conf ${path_conf} -conf_format ${conf_format} -bytes_skip ${bytes_skip}\
 -path_conf_output ${path_conf_output} -path_functional_output ${path_functional_output}\
 -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size}
