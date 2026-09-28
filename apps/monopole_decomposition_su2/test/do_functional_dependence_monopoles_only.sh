#!/bin/bash
# path_conf="../../../tests/confs/MAG/su2/su2_suzuki/24^4/beta2.4/steps=0/conf_0001"
# path_conf="../../../tests/confs/su2/gluodynamics/64^4/beta2.9/CONF0002"
path_conf="../../../tests/confs/su2/gluodynamics/48^4/beta2.7/CON_MC_001.LAT"
conf_format="lexicographical"
file_precision="double"
bytes_skip=8
path_inverse_laplacian="./result/inverse_laplacian_48x48"
x_size=48
y_size=48
z_size=48
t_size=48
copies_required=2
mag_steps=0
path_functional_output="./result/functional"
path_clusters_unwrapped_abelian_output="./result/clusters_unwrapped_abelian"
path_clusters_unwrapped_monopole_output="./result/clusters_unwrapped_monopole"
path_clusters_unwrapped_monopoless_output="./result/clusters_unwrapped_monopoless"
path_clusters_wrapped_abelian_output="./result/clusters_wrapped_abelian"
path_clusters_wrapped_monopole_output="./result/clusters_wrapped_monopole"
path_clusters_wrapped_monopoless_output="./result/clusters_wrapped_monopoless"
path_windings_abelian_output="./result/windings_abelian"
path_windings_monopole_output="./result/windings_monopole"
path_windings_monopoless_output="./result/windings_monopoless"

../functional_dependence_monopoles_only_test --path_conf ${path_conf} --conf_format ${conf_format} --file_precision ${file_precision} --bytes_skip ${bytes_skip}\
 --path_inverse_laplacian ${path_inverse_laplacian}  \
 --copies_required ${copies_required} --path_functional_output ${path_functional_output} --mag_steps ${mag_steps} \
 --path_clusters_unwrapped_abelian_output ${path_clusters_unwrapped_abelian_output} --path_clusters_unwrapped_monopole_output ${path_clusters_unwrapped_monopole_output} --path_clusters_unwrapped_monopoless_output ${path_clusters_unwrapped_monopoless_output} \
 --path_clusters_wrapped_abelian_output ${path_clusters_wrapped_abelian_output} --path_clusters_wrapped_monopole_output ${path_clusters_wrapped_monopole_output} --path_clusters_wrapped_monopoless_output ${path_clusters_wrapped_monopoless_output} \
 --path_windings_abelian_output ${path_windings_abelian_output} --path_windings_monopole_output ${path_windings_monopole_output} --path_windings_monopoless_output ${path_windings_monopoless_output} \
 --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}