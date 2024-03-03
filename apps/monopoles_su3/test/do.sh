#!/bin/bash
path_conf="../../../tests/confs/MAG/su3/gluodynamics/16^4/beta6.0/steps_4000/copies=16/0.1/conf_gaugefixed_01001.lime"
conf_format=ildg
convert=1
bytes_skip=0
path_output_clusters_unwrapped="./result/clusters_unwrapped_01001"
path_output_clusters_wrapped="./result/clusters_wrapped_01001"
path_output_windings="./result/windings_01001"
path_output_monopoles="./result/monopoles_01001"
x_size=16
y_size=16
z_size=16
t_size=16

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_output_clusters_unwrapped ${path_output_clusters_unwrapped}\
    -path_output_clusters_wrapped ${path_output_clusters_wrapped} -path_output_windings ${path_output_windings} -path_output_monopoles ${path_output_monopoles} \
    -bytes_skip ${bytes_skip} -convert ${convert} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size}"

../monopoles_su3_test ${parameters}