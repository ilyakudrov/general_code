#!/bin/bash
path_inverse_laplacian="./result/inverse_laplacian_64x4"
x_size=64
y_size=64
z_size=64
t_size=4

../inverse_laplacian_test --path_inverse_laplacian ${path_inverse_laplacian} --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}
