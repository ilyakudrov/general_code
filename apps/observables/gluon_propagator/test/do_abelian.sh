#!/bin/bash
conf_format="lexicographical"
file_precision="float"
# conf_format="double_qc2dstag"
conf_path="../../../../tests/confs/su2/gluodynamics/66^3x8/beta2.701/CONF_B2701L66T08"
# conf_path="../../../../tests/confs/su2/gluodynamics/32^3x8/beta2.779/CONF0001"
bytes_skip=76
# bytes_skip=0
beta=2.701
convert=1
L_spat=66
L_time=8
path_inverse_laplacian="../../../monopole_decomposition_su2/test/result/ALPHA66x8_d.LAT"
output_path_abelian=./result/gluon_propagator_abelian
output_path_monopole=./result/gluon_propagator_monopole
output_path_photon=./result/gluon_propagator_photon
parameters="--conf_format ${conf_format} --conf_path ${conf_path} --file_precision ${file_precision} --bytes_skip ${bytes_skip} --convert ${convert}\
    --path_inverse_laplacian ${path_inverse_laplacian} --output_path_abelian ${output_path_abelian} --output_path_monopole ${output_path_monopole} \
    --output_path_photon ${output_path_photon} --L_spat ${L_spat} --L_time ${L_time} --output_path ${output_path} --beta ${beta}"

../gluon_propagator_abelian_test ${parameters}