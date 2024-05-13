#!/bin/bash
conf_format="double"
conf_path="../../../tests/confs/monopole/su2/qc2dstag/40^4/mu0.00/conf_monopole_0201"
bytes_skip=0
# matrix_type="su2"
matrix_type="abelian"
L_spat=40
L_time=40
conf_path_output=./result/conf_monopole_qc2dstag0201
convert=0
parameters="-conf_format ${conf_format} -conf_path ${conf_path} -bytes_skip ${bytes_skip}\
    -convert ${convert} -L_spat ${L_spat} -L_time ${L_time} -conf_path_output ${conf_path_output}"

../transform_to_qc2dstag_${matrix_type}_test ${parameters}