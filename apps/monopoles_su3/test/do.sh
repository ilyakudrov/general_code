#!/bin/bash
# path_conf="../../../tests/confs/MA_gauge/su3/gluodynamics/24^4/beta6.0/CONFDP_gaugefixed_0001"
path_conf="../../../tests/confs/decomposed/monopole/gluodynamics/36^4/beta6.3/steps_25/copies=4/conf_monopole_0001"
conf_format=doble_abelian_su3
path_output_clusters_unwrapped="./result/clusters_unwrapped_0001"
path_output_clusters_wrapped="./result/clusters_wrapped_0001"
path_output_windings="./result/windings_0001"
path_output_monopoles="./result/monopoles_0001"
bytes_skip=0
x_size=36
y_size=36
z_size=36
t_size=36

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_output_clusters_unwrapped ${path_output_clusters_unwrapped}\
    -path_output_clusters_wrapped ${path_output_clusters_wrapped} -path_output_windings ${path_output_windings} -path_output_monopoles ${path_output_monopoles} \
    -bytes_skip ${bytes_skip} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size}"

../monopoles_su3_test ${parameters}