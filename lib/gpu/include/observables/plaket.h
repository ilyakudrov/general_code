#pragma once

#include "../../../cpu/include/link.h"
#include "../general/matrix.h"
#include "stdio.h"

__constant__ unsigned int steps[4];
__constant__ unsigned int arr_mu[6];
__constant__ unsigned int arr_nu[6];

namespace plaket_gpu {
__global__ void plaket(matrix_gpu::su3 *conf, double *traces,
                       unsigned int *indices);
__device__ double plaket_site(matrix_gpu::su3 *conf, unsigned int *indices,
                              unsigned short int mu, unsigned short int nu);
__host__ void make_plaket_indices(unsigned int *indices, int lattice_size[4]);

__global__ void __launch_bounds__(128, 4)
    plaket_test1(matrix_gpu::su3 *conf, double *traces, unsigned int *indices);

__device__ double plaket_site_test2(const matrix_gpu::su3 *conf,
                                    const unsigned int *indices);

__global__ void __launch_bounds__(128, 4)
    plaket_test2(const matrix_gpu::su3 *conf, double *traces,
                 const unsigned int *indices);

__host__ void make_plaket_indices_test2(unsigned int *indices,
                                        int lattice_size[4]);

__global__ void __launch_bounds__(128, 4)
    plaket_test3(const matrix_gpu::su3 *conf, double *traces);

__global__ void __launch_bounds__(128, 6)
    plaket_test4(matrix_gpu::su3 *conf, double *traces, unsigned int *indices);

__global__ void plaket(matrix_gpu::su3 *conf, double *traces,
                       unsigned int *indices) {
  int lat_index = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int *ind_tmp = &indices[lat_index * 5];
  double trace = 0;

  trace += plaket_site(conf, ind_tmp, 0, 1);
  trace += plaket_site(conf, ind_tmp, 0, 2);
  trace += plaket_site(conf, ind_tmp, 0, 3);
  trace += plaket_site(conf, ind_tmp, 1, 2);
  trace += plaket_site(conf, ind_tmp, 1, 3);
  trace += plaket_site(conf, ind_tmp, 2, 3);
  traces[lat_index] = trace;
}

__device__ double plaket_site(matrix_gpu::su3 *conf, unsigned int *indices,
                              unsigned short int mu, unsigned short int nu) {
  matrix_gpu::su3 A = conf[indices[0] + mu] * conf[indices[mu + 1] + nu];
  A = A ^ conf[indices[nu + 1] + mu];
  return A.multiply_tr(conf[indices[0] + nu]);
}

__host__ void make_plaket_indices(unsigned int *indices, int lattice_size[4]) {
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

__global__ void plaket_test1(matrix_gpu::su3 *conf, double *traces,
                             unsigned int *indices) {
  unsigned int lat_index = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int plaket_num = lat_index % 6;
  unsigned short int mu = 0;
  unsigned short int nu = 0;
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

__device__ double plaket_site_test2(const matrix_gpu::su3 *conf,
                                    const unsigned int *indices) {
  matrix_gpu::su3 A = conf[indices[0]] * conf[indices[1]];
  A = A ^ conf[indices[2]];
  return A.multiply_tr(conf[indices[3]]);
}

__global__ void plaket_test2(const matrix_gpu::su3 *conf, double *traces,
                             const unsigned int *indices) {
  int lat_index = blockIdx.x * blockDim.x + threadIdx.x;

  traces[lat_index] = plaket_site_test2(conf, &indices[lat_index * 4]);
}

__host__ void make_plaket_indices_test2(unsigned int *indices,
                                        int lattice_size[4]) {
  link1 link(lattice_size[0], lattice_size[1], lattice_size[2],
             lattice_size[3]);
  int place = 0;
  for (int t = 0; t < lattice_size[0]; t++) {
    for (int z = 0; z < lattice_size[0]; z++) {
      for (int y = 0; y < lattice_size[0]; y++) {
        for (int x = 0; x < lattice_size[0]; x++) {
          link.go_update(x, y, z, t);
          for (int mu = 0; mu < 3; mu++) {
            for (int nu = mu + 1; nu < 4; nu++) {
              indices[place] = link.place + mu;
              place += 1;
              link.move(mu, 1);
              indices[place] = link.place + nu;
              place += 1;
              link.move(mu, -1);
              link.move(nu, 1);
              indices[place] = link.place + mu;
              place += 1;
              link.move(nu, -1);
              indices[place] = link.place + nu;
              place += 1;
            }
          }
        }
      }
    }
  }
}

__global__ void plaket_test3(const matrix_gpu::su3 *conf, double *traces) {
  int lat_index = blockIdx.x * blockDim.x + threadIdx.x;
  int plaket_num = lat_index % 6;
  int lat_site = (lat_index / 6) * 4;
  unsigned int mu = arr_mu[plaket_num];
  unsigned int nu = arr_nu[plaket_num];
  unsigned int step_mu = steps[mu];
  unsigned int step_nu = steps[nu];

  matrix_gpu::su3 A = conf[lat_site + mu] * conf[lat_site + step_mu + nu];
  A = A ^ conf[lat_site + step_nu + mu];
  traces[lat_index] = A.multiply_tr(conf[lat_site + nu]);
}

__global__ void plaket_test4(matrix_gpu::su3 *conf, double *traces,
                             unsigned int *indices) {
  unsigned int lat_index =
      blockIdx.x * blockDim.x * blockDim.y + blockDim.x * threadIdx.y;
  unsigned short int mu = 0;
  unsigned short int nu = 0;
  if (threadIdx.x == 0) {
    mu = 0;
    nu = 1;
  } else if (threadIdx.x == 1) {
    mu = 0;
    nu = 2;
  } else if (threadIdx.x == 2) {
    mu = 0;
    nu = 3;
  } else if (threadIdx.x == 3) {
    mu = 1;
    nu = 2;
  } else if (threadIdx.x == 4) {
    mu = 1;
    nu = 3;
  } else if (threadIdx.x == 5) {
    mu = 2;
    nu = 3;
  }

  traces[lat_index + threadIdx.x] =
      plaket_site(conf, &indices[(lat_index / 6) * 5], mu, nu);
}

} // namespace plaket_gpu