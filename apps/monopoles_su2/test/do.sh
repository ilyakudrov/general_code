#!/bin/bash
# path_conf="../../../tests/confs/MAG/su2/qc2dstag/40^4/mu0.00/conf_abelian_0201"
# path_conf="../../../tests/confs/monopole/su2/qc2dstag/40^4/mu0.00/conf_monopole_0201"
path_conf="../../../tests/confs/monopoless/su2/qc2dstag/40^4/mu0.00/conf_monopoless_0201"
conf_format=lexicographical
file_precision=double
convert=1
path_output_clusters_unwrapped="./result/clusters_unwrapped2"
path_output_clusters_wrapped="./result/clusters_wrapped2"
path_output_windings="./result/windings2"
path_output_monopoles="./result/monopoles2"
bytes_skip=0
x_size=40
y_size=40
z_size=40
t_size=40

parameters="--path_conf ${path_conf} --file_precision ${file_precision} --conf_format ${conf_format} --path_output_clusters_unwrapped ${path_output_clusters_unwrapped}\
    --path_output_clusters_wrapped ${path_output_clusters_wrapped} --path_output_windings ${path_output_windings} --path_output_monopoles ${path_output_monopoles} \
    --bytes_skip ${bytes_skip} --convert ${convert} --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}"

../monopoles_su2_test ${parameters}