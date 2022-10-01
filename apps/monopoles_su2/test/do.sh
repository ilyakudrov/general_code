#!/bin/bash
path_conf="../../../tests/confs/MA_gauge/su2/su2_suzuki/48^4/beta2.8/T_step=0.0001/T_final=0.5/OR_steps=4/conf_0001"
conf_format=double
path_output_clusters_unwrapped="./result/clusters_unwrapped_0001"
path_output_clusters_wrapped="./result/clusters_wrapped_0001"
path_output_windings="./result/windings_0001"
path_output_monopoles="./result/monopoles_0001"
bytes_skip=0
x_size=48
y_size=48
z_size=48
t_size=48

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_output_clusters_unwrapped ${path_output_clusters_unwrapped}\
    -path_output_clusters_wrapped ${path_output_clusters_wrapped} -path_output_windings ${path_output_windings} -path_output_monopoles ${path_output_monopoles} \
    -bytes_skip ${bytes_skip} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size}"

../monopoles_su2_test ${parameters}