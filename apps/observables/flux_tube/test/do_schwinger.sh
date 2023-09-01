#!/bin/bash
conf_format_wilson="double"
# conf_format_wilson="double_qc2dstag"
# conf_path_wilson="../../../../tests/confs/su2/smeared/qc2dstag/40^4/mu0.00/original/HYP1_alpha=1_1_0.5_APE_alpha=0.5/smeared_0201"
conf_path_wilson="../../../../tests/confs/su2/ml5/beta2.1_1conf.ml5"
# conf_path_wilson="../../../../tests/confs/su2/qc2dstag/40^4/mu0.00/CONF0201"
# conf_path_wilson="../../../../tests/confs/MA_gauge/su2/qc2dstag/40^4/mu0.00/conf_abelian_0201"
bytes_skip_wilson=0
matrix_type_wilson="su2"
convert_wilson=0
conf_format_plaket="double"
# conf_format_plaket="double_qc2dstag"
# conf_path_plaket="../../../../tests/confs/su2/qc2dstag/40^4/mu0.00/CONF0201"
# conf_path_plaket="../../../../tests/confs/su2/smeared/qc2dstag/40^4/mu0.00/original/HYP1_alpha=1_1_0.5_APE_alpha=0.5/smeared_0201"
conf_path_plaket="../../../../tests/confs/su2/ml5/beta2.1_1conf.ml5"
# conf_path_plaket="../../../../tests/confs/MA_gauge/su2/qc2dstag/40^4/mu0.00/conf_abelian_0201"
bytes_skip_plaket=0
matrix_type_plaket="su2"
convert_plaket=0
L_spat=20
L_time=20
output_path_electric_long_l=./result/flux_tube_schwinger_electric_long_l
output_path_electric_long_tr=./result/flux_tube_schwinger_electric_long_tr
output_path_magnetic_long_l=./result/flux_tube_schwinger_magnetic_long_l
output_path_magnetic_long_tr=./result/flux_tube_schwinger_magnetic_long_tr
output_path_electric_trans_l=./result/flux_tube_schwinger_electric_trans_l
output_path_electric_trans_tr=./result/flux_tube_schwinger_electric_trans_tr
output_path_magnetic_trans_l=./result/flux_tube_schwinger_magnetic_trans_l
output_path_magnetic_trans_tr=./result/flux_tube_schwinger_magnetic_trans_tr
T_min=4
T_max=4
R_min=5
R_max=5
d_ouside=10
d_max=10
parameters="-conf_format_plaket ${conf_format_plaket} -conf_path_plaket ${conf_path_plaket}\
    -bytes_skip_plaket ${bytes_skip_plaket} -convert_plaket ${convert_plaket}\
    -conf_format_wilson ${conf_format_wilson} -conf_path_wilson ${conf_path_wilson}\
    -bytes_skip_wilson ${bytes_skip_wilson} -convert_wilson ${convert_wilson}\
    -output_path_electric_long_l ${output_path_electric_long_l}\
    -output_path_electric_long_tr ${output_path_electric_long_tr}\
    -output_path_electric_trans_l ${output_path_electric_trans_l}\
    -output_path_electric_trans_tr ${output_path_electric_trans_tr}\
    -output_path_magnetic_long_l ${output_path_magnetic_long_l}\
    -output_path_magnetic_long_tr ${output_path_magnetic_long_tr}\
    -output_path_magnetic_trans_l ${output_path_magnetic_trans_l}\
    -output_path_magnetic_trans_tr ${output_path_magnetic_trans_tr}\
    -L_spat ${L_spat} -L_time ${L_time} -d_ouside ${d_ouside} -d_max ${d_max}\
    -T_min ${T_min} -T_max ${T_max} -R_min ${R_min} -R_max ${R_max}"

../flux_tube_schwinger_${matrix_type_plaket}_${matrix_type_wilson}_test ${parameters}