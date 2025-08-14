#!/bin/bash
path_inverse_laplacian="./result/ALPHA24x24_d.LAT"
x_size=24
y_size=24
z_size=24
t_size=24

../inverse_laplacian_test --path_inverse_laplacian ${path_inverse_laplacian} --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}