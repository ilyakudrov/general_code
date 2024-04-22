#!/bin/bash
path_conf="../../../tests/confs/su3/gluodynamics/abelian/24^4/beta6.0/steps_100/copies=20/0.01/conf_gaugefixed_0001_1"
conf_format=double_qc2dstag
convert=1
bytes_skip=0
path_output_clusters_unwrapped="./result/clusters_unwrapped_0001"
path_output_clusters_wrapped="./result/clusters_wrapped_0001"
path_output_windings="./result/windings_0001"
path_output_monopoles="./result/monopoles_0001"
x_size=24
y_size=24
z_size=24
t_size=24

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_output_clusters_unwrapped ${path_output_clusters_unwrapped}\
    -path_output_clusters_wrapped ${path_output_clusters_wrapped} -path_output_windings ${path_output_windings} -path_output_monopoles ${path_output_monopoles} \
    -bytes_skip ${bytes_skip} -convert ${convert} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size}"

../monopoles_su3_test ${parameters}