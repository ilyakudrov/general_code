#!/bin/bash
path_inverse_laplacian="./result/ALPHA66x8_d.LAT"
x_size=66
y_size=66
z_size=66
t_size=8

../inverse_laplacian_test --path_inverse_laplacian ${path_inverse_laplacian} --x_size ${x_size} --y_size ${y_size} --z_size ${z_size} --t_size ${t_size}