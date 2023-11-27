#!/bin/bash
# path_conf="../../../tests/confs/MA_gauge/su2/su2_suzuki/24^4/beta2.6/T_step=0.0001/conf_0001"
# path_conf="../../../tests/confs/monopole/su2/su2_suzuki/24^4/beta2.6/T_step=0.0001/conf_monopole_0001"
# path_conf="../../../tests/confs/monopoless/su2/su2_suzuki/24^4/beta2.6/T_step=0.0001/conf_monopoless_0001"
# path_conf="../../../tests/confs/MA_gauge/su2/gluodynamics/24^3x6/beta2.35/CON_fxd_MAG_U1_001.LAT"
# path_conf="../../../tests/confs/monopole/su2/gluodynamics/24^3x6/beta2.35/CON_MON_MAG_001.LAT"
path_conf="../../../tests/confs/monopoless/su2/gluodynamics/24^3x6/beta2.35/CON_OFF_MAG_001.LAT"
# conf_format=double
conf_format=double
convert=1
path_output_clusters_unwrapped="./result/clusters_unwrapped_0001"
path_output_clusters_wrapped="./result/clusters_wrapped_0001"
path_output_windings="./result/windings_0001"
path_output_monopoles="./result/monopoles_0001"
bytes_skip=4
x_size=24
y_size=24
z_size=24
t_size=6

parameters="-path_conf ${path_conf} -conf_format ${conf_format} -path_output_clusters_unwrapped ${path_output_clusters_unwrapped}\
    -path_output_clusters_wrapped ${path_output_clusters_wrapped} -path_output_windings ${path_output_windings} -path_output_monopoles ${path_output_monopoles} \
    -bytes_skip ${bytes_skip} -convert ${convert} -x_size ${x_size} -y_size ${y_size} -z_size ${z_size} -t_size ${t_size}"

../monopoles_su2_test ${parameters}