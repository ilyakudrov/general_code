#!/bin/bash
path_inverse_laplacian="./result/inverse_laplacian_48x12"
x_size=48
y_size=48
z_size=48
t_size=12

../inverse_laplacian_test --path_inverse_laplacian ${path_inverse_laplacian} --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}
