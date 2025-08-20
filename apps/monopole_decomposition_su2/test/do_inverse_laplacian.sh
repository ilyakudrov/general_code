#!/bin/bash
path_inverse_laplacian="./result/inverse_laplacian_32x32"
x_size=32
y_size=32
z_size=32
t_size=32

../inverse_laplacian_test --path_inverse_laplacian ${path_inverse_laplacian} --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}
