#pragma once

#include "../../../cpu/include/link.h"
#include "../general/matrix.h"
#include "stdio.h"

namespace plaket_gpu {
__global__ void plaket(matrix_gpu::su3 *conf, double *traces, int *indices);
__device__ double plaket_site(matrix_gpu::su3 *conf, int *indices, int mu,
                              int nu);
__host__ void make_plaket_indices(int *indices, int lattice_size[4]);

__global__ void __launch_bounds__(128, 4)
    plaket_test1(matrix_gpu::su3 *conf, double *traces, int *indices);

__global__ void plaket(matrix_gpu::su3 *conf, double *traces, int *indices) {
  int lat_index = blockIdx.x * blockDim.x + threadIdx.x;
  int *ind_tmp = &indices[lat_index * 5];
  double trace = 0;

  // __shared__ su3 A1 = conf[0];

  trace += plaket_site(conf, ind_tmp, 0, 1);
  trace += plaket_site(conf, ind_tmp, 0, 2);
  trace += plaket_site(conf, ind_tmp, 0, 3);
  trace += plaket_site(conf, ind_tmp, 1, 2);
  trace += plaket_site(conf, ind_tmp, 1, 3);
  trace += plaket_site(conf, ind_tmp, 2, 3);
  traces[lat_index] = trace;
}

__device__ double plaket_site(matrix_gpu::su3 *conf, int *indices, int mu,
                              int nu) {
  matrix_gpu::su3 A = conf[indices[0] + mu] * conf[indices[mu + 1] + nu];
  A = A ^ conf[indices[nu + 1] + mu];
  return A.multiply_tr(conf[indices[0] + nu]);
}

__host__ void make_plaket_indices(int *indices, int lattice_size[4]) {
  link1 link(lattice_size[0], lattice_size[1], lattice_size[2],
             lattice_size[3]);
  int place = 0;
  for (int t = 0; t < lattice_size[0]; t++) {
    for (int z = 0; z < lattice_size[0]; z++) {
      for (int y = 0; y < lattice_size[0]; y++) {
        for (int x = 0; x < lattice_size[0]; x++) {
          link.go_update(x, y, z, t);
          indices[place] = link.place;
          place += 1;
          for (int mu = 0; mu < 4; mu++) {
            link.move(mu, 1);
            indices[place] = link.place;
            place += 1;
            link.move(mu, -1);
          }
        }
      }
    }
  }
}

__global__ void plaket_test1(matrix_gpu::su3 *__restrict__ conf, double *traces,
                             int *indices) {
  int lat_index = blockIdx.x * blockDim.x + threadIdx.x;
  int plaket_num = lat_index % 6;
  int mu = 0;
  int nu = 0;
  if (plaket_num == 0) {
    mu = 0;
    nu = 1;
  } else if (plaket_num == 1) {
    mu = 0;
    nu = 2;
  } else if (plaket_num == 2) {
    mu = 0;
    nu = 3;
  } else if (plaket_num == 3) {
    mu = 1;
    nu = 2;
  } else if (plaket_num == 4) {
    mu = 1;
    nu = 3;
  } else if (plaket_num == 5) {
    mu = 2;
    nu = 3;
  }

  traces[lat_index] = plaket_site(conf, &indices[(lat_index / 6) * 5], mu, nu);
}

} // namespace plaket_gpu