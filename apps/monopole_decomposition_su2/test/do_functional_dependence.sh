#!/bin/bash
path_conf="../../../tests/confs/su2/gluodynamics/32^3x8/beta2.779/CONF0001"
# path_conf="../../../tests/confs/su2/gluodynamics/64^4/beta2.9/CONF0002"
conf_format="qcdstag"
file_precision="double"
bytes_skip=0
path_inverse_laplacian="./result/inverse_laplacian_32x8"
N_dir_gevp=1
HYP_alpha1=1
HYP_alpha2=1
HYP_alpha3=0.5
APE_alpha=0.5
APE_steps=31
calculation_step_APE=10
calculation_APE_start=11
HYP_steps=0
x_size=32
y_size=32
z_size=32
t_size=8
copies_required=11
mag_steps=100
path_functional_output="./result/functional"
path_wilson_loops_abelian_output="./result/wilson_loops_abelian"
path_wilson_loops_monopole_output="./result/wilson_loops_monopole"
path_clusters_unwrapped_abelian_output="./result/clusters_unwrapped_abelian"
path_clusters_unwrapped_monopole_output="./result/clusters_unwrapped_monopole"
path_clusters_unwrapped_monopoless_output="./result/clusters_unwrapped_monopoless"
path_clusters_wrapped_abelian_output="./result/clusters_wrapped_abelian"
path_clusters_wrapped_monopole_output="./result/clusters_wrapped_monopole"
path_clusters_wrapped_monopoless_output="./result/clusters_wrapped_monopoless"
path_windings_abelian_output="./result/windings_abelian"
path_windings_monopole_output="./result/windings_monopole"
path_windings_monopoless_output="./result/windings_monopoless"

../functional_dependence_test --path_conf ${path_conf} --conf_format ${conf_format} --file_precision ${file_precision} --bytes_skip ${bytes_skip}\
 --path_inverse_laplacian ${path_inverse_laplacian} --N_dir_gevp ${N_dir_gevp} --HYP_alpha1 ${HYP_alpha1} --HYP_alpha2 ${HYP_alpha2} --HYP_alpha3 ${HYP_alpha3} \
 --copies_required ${copies_required} --path_functional_output ${path_functional_output} --mag_steps ${mag_steps} \
 --path_wilson_loops_abelian_output ${path_wilson_loops_abelian_output} --path_wilson_loops_monopole_output ${path_wilson_loops_monopole_output} \
 --path_clusters_unwrapped_abelian_output ${path_clusters_unwrapped_abelian_output} --path_clusters_unwrapped_monopole_output ${path_clusters_unwrapped_monopole_output} --path_clusters_unwrapped_monopoless_output ${path_clusters_unwrapped_monopoless_output} \
 --path_clusters_wrapped_abelian_output ${path_clusters_wrapped_abelian_output} --path_clusters_wrapped_monopole_output ${path_clusters_wrapped_monopole_output} --path_clusters_wrapped_monopoless_output ${path_clusters_wrapped_monopoless_output} \
 --path_windings_abelian_output ${path_windings_abelian_output} --path_windings_monopole_output ${path_windings_monopole_output} --path_windings_monopoless_output ${path_windings_monopoless_output} \
 --APE_alpha ${APE_alpha} --APE_steps ${APE_steps} --calculation_step_APE ${calculation_step_APE} --calculation_APE_start ${calculation_APE_start} --HYP_steps ${HYP_steps} \
 --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}