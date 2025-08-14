#!/bin/bash
path_inverse_laplacian="./result/inverse_laplacian_40x40"
x_size=40
y_size=40
z_size=40
t_size=40

../inverse_laplacian_test --path_inverse_laplacian ${path_inverse_laplacian} --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}