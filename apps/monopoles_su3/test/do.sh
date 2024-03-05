#!/bin/bash
path_conf="../../../tests/confs/MAG/su3/gluodynamics/24^4/beta6.0/steps_500/copies=4/conf_gaugefixed_1486"
conf_format=double_qc2dstag
convert=1
bytes_skip=0
path_output_clusters_unwrapped="./result/clusters_unwrapped_1486"
path_output_clusters_wrapped="./result/clusters_wrapped_1486"
path_output_windings="./result/windings_1486"
path_output_monopoles="./result/monopoles_1486"
x_size=24
y_size=24
z_size=24
t_size=24

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_output_clusters_unwrapped ${path_output_clusters_unwrapped}\
    -path_output_clusters_wrapped ${path_output_clusters_wrapped} -path_output_windings ${path_output_windings} -path_output_monopoles ${path_output_monopoles} \
    -bytes_skip ${bytes_skip} -convert ${convert} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size}"

../monopoles_su3_test ${parameters}