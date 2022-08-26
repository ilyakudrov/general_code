#!/bin/bash
path_conf="../../../tests/confs/su3/mag/gluodynamics/36^4/beta6.3/CONFDP_gaugefixed_0001"
conf_format=double_qc2dstag
path_output_clusters="./result/clusters_36_0001"
path_output_windings="./result/windings_36_0001"
path_output_monopoles="./result/monopoles_36_0001"
bytes_skip=0
x_size=36
y_size=36
z_size=36
t_size=36

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_output_clusters ${path_output_clusters} \
    -path_output_windings ${path_output_windings} -path_output_monopoles ${path_output_monopoles} \
    -bytes_skip ${bytes_skip} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size}"

../monopoles_su3_test ${parameters}