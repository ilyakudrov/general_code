#!/bin/bash
conf_format="lexicographical"
file_precision="double"
# conf_format="double_qc2dstag"
conf_path="../../../../tests/confs/su2/su2_suzuki/48^4/beta2.6/CON_MC_001.LAT"
# conf_path="../../../../tests/confs/su2/gluodynamics/32^3x8/beta2.779/CONF0001"
bytes_skip=4
# bytes_skip=0
matrix_type="su2"
beta=2.701
convert=0
L_spat=48
L_time=48
output_path=./result/gluon_propagator
parameters="--conf_format ${conf_format} --conf_path ${conf_path} --file_precision ${file_precision} --bytes_skip ${bytes_skip} --convert ${convert}\
    --L_spat ${L_spat} --L_time ${L_time} --output_path ${output_path} --beta ${beta}"

../gluon_propagator_${matrix_type}_test ${parameters}
