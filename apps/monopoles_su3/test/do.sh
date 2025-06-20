#!/bin/bash
path_conf="../../../tests/confs/monopoless/su3/gluodynamics/24^4/beta6.0/steps_0/copies=20/conf_monopoless_0001_0"
# path_conf="/home/ilya/soft/lattice/general_code/apps/monopole_decomposition_su3/test/result/conf_monopoless"
conf_format=lexicographical
file_precision="double"
convert=1
bytes_skip=0
path_output_clusters_unwrapped="./result/clusters_unwrapped_2"
path_output_clusters_wrapped="./result/clusters_wrapped_2"
path_output_windings="./result/windings_2"
path_output_monopoles="./result/monopoles_2"
x_size=24
y_size=24
z_size=24
t_size=24
parameters="--path_conf ${path_conf} --file_precision ${file_precision} --conf_format ${conf_format} --path_output_clusters_unwrapped ${path_output_clusters_unwrapped}\
    --path_output_clusters_wrapped ${path_output_clusters_wrapped} --path_output_windings ${path_output_windings} --path_output_monopoles ${path_output_monopoles} \
    --bytes_skip ${bytes_skip} --convert ${convert} --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}"

../monopoles_su3_test ${parameters}