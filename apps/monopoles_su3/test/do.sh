#!/bin/bash
# path_conf="../../../tests/confs/MA_gauge/su3/gluodynamics/24^4/beta6.0/CONFDP_gaugefixed_0001"
# path_conf="../../../tests/confs/monopole/su3/QCD/140MeV/nt18/steps_500/copies=1/conf_monopole_0501"
# path_conf="../../../tests/confs/monopoless/su3/QCD/140MeV/nt18/steps_500/copies=1/conf_monopoless_0501"
# path_conf="../../../tests/confs/MA_gauge/su3/gluodynamics/24^4/beta6.0/steps_500/copies=4/conf_gaugefixed_0001"
# path_conf="../../../tests/confs/monopole/su3/gluodynamics/24^4/beta6.0/steps_500/copies=4/conf_monopole_0001"
path_conf="../../../tests/confs/monopoless/su3/monopoless/su3/gluodynamics/24^4/beta6.0/steps_500/copies=4/conf_monopoless_0001"
# path_conf="../../../tests/confs/MA_gauge/su3/QCD/140MeV/nt18/steps_500/copies=1/conf_gaugefixed_0501"
conf_format=double
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