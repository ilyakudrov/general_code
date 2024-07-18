#!/bin/bash
path_conf="../../../tests/confs/MAG/su3/gluodynamics/40^4/beta6.4/steps_0/copies=20/s1/conf_gaugefixed_0002_1"
conf_format=double_qc2dstag
convert=1
bytes_skip=0
path_output_clusters_unwrapped="./result/clusters_unwrapped_0001"
path_output_clusters_wrapped="./result/clusters_wrapped_0001"
path_output_windings="./result/windings_0001"
path_output_monopoles="./result/monopoles_0001"
x_size=40
y_size=40
z_size=40
t_size=40

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_output_clusters_unwrapped ${path_output_clusters_unwrapped}\
    -path_output_clusters_wrapped ${path_output_clusters_wrapped} -path_output_windings ${path_output_windings} -path_output_monopoles ${path_output_monopoles} \
    -bytes_skip ${bytes_skip} -convert ${convert} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size}"

../monopoles_su3_test ${parameters}