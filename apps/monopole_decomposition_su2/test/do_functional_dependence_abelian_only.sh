#!/bin/bash
# path_conf="../../../tests/confs/MAG/su2/su2_suzuki/24^4/beta2.4/steps=0/conf_0001"
path_conf="../../../tests/confs/su2/gluodynamics/64^4/beta2.9/CONF0002"
# path_conf="../../../tests/confs/su2/gluodynamics/48^4/beta2.7/CON_MC_001.LAT"
conf_format="qcdstag"
file_precision="double"
bytes_skip=0
x_size=64
y_size=64
z_size=64
t_size=64
copies_required=2
mag_steps=100
path_functional_output="./result/functional"
path_clusters_unwrapped_abelian_output="./result/clusters_unwrapped_abelian"
path_clusters_wrapped_abelian_output="./result/clusters_wrapped_abelian"
path_windings_abelian_output="./result/windings_abelian"

../functional_dependence_abelian_only_test --path_conf ${path_conf} --conf_format ${conf_format} --file_precision ${file_precision} --bytes_skip ${bytes_skip}\
 --copies_required ${copies_required} --path_functional_output ${path_functional_output} --mag_steps ${mag_steps} \
 --path_clusters_unwrapped_abelian_output ${path_clusters_unwrapped_abelian_output} \
 --path_clusters_wrapped_abelian_output ${path_clusters_wrapped_abelian_output} \
 --path_windings_abelian_output ${path_windings_abelian_output} \
 --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}