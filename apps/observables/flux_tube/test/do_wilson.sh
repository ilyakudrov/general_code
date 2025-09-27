#!/bin/bash
conf_format_wilson="lexicographical"
file_precision_wilson="double"
conf_path_wilson="../../../../tests/confs/su2/su2_suzuki/24^4/beta2.4/CON_MC_001.LAT"
bytes_skip_wilson=8
matrix_type_wilson="su2"
convert_wilson=0
conf_format_plaket="lexicographical"
file_precision_plaket="double"
conf_path_plaket="../../../../tests/confs/su2/su2_suzuki/24^4/beta2.4/CON_MC_001.LAT"
bytes_skip_plaket=8
matrix_type_plaket="su2"
convert_plaket=0
L_spat=24
L_time=24
output_path_electric_long=./result/flux_tube_electric_long
output_path_magnetic_long=./result/flux_tube_magnetic_long
output_path_electric_trans=./result/flux_tube_electric_trans
output_path_magnetic_trans=./result/flux_tube_magnetic_trans
T_min=4
T_max=12
R_min=4
R_max=12
parameters="--conf_format_plaket ${conf_format_plaket} --file_precision_wilson ${file_precision_wilson} --conf_path_plaket ${conf_path_plaket}\
    --bytes_skip_plaket ${bytes_skip_plaket} --convert_plaket ${convert_plaket}\
    --conf_format_wilson ${conf_format_wilson} --file_precision_plaket ${file_precision_plaket} --conf_path_wilson ${conf_path_wilson}\
    --bytes_skip_wilson ${bytes_skip_wilson} --convert_wilson ${convert_wilson}\
    --output_path_electric_long ${output_path_electric_long} --output_path_magnetic_long ${output_path_magnetic_long}\
    --output_path_electric_trans ${output_path_electric_trans} --output_path_magnetic_trans ${output_path_magnetic_trans}\
    --L_spat ${L_spat} --L_time ${L_time}\
    --T_min ${T_min} --T_max ${T_max} --R_min ${R_min} --R_max ${R_max}"

../flux_tube_wilson_${matrix_type_plaket}_${matrix_type_wilson}_test ${parameters}