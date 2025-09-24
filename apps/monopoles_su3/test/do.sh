#!/bin/bash
path_conf="../../../tests/confs/MAG/su3/gluodynamics/48^3x12/beta6.257/steps_0/conf_gaugefixed_0001"
# path_conf="/home/ilya/soft/lattice/general_code/apps/monopole_decomposition_su3/test/result/conf_monopoless"
conf_format=qcdstag
file_precision="double"
convert=1
bytes_skip=0
path_output_clusters_unwrapped="./result/clusters_unwrapped"
path_output_clusters_wrapped="./result/clusters_wrapped"
path_output_windings="./result/windings"
path_output_monopoles="./result/monopoles"
x_size=48
y_size=48
z_size=48
t_size=12
parameters="--path_conf ${path_conf} --file_precision ${file_precision} --conf_format ${conf_format} --path_output_clusters_unwrapped ${path_output_clusters_unwrapped}\
    --path_output_clusters_wrapped ${path_output_clusters_wrapped} --path_output_windings ${path_output_windings} --path_output_monopoles ${path_output_monopoles} \
    --bytes_skip ${bytes_skip} --convert ${convert} --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}"

../monopoles_su3_test ${parameters}