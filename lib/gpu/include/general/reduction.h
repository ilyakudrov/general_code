#pragma once

// #ifndef THREADS_PER_BLOC
// #define THREADS_PER_BLOCK 256
// #endif

#include <__clang_cuda_device_functions.h>
namespace reduction {

// template <unsigned int blockSize>
// __device__ void warpReduce(volatile double *sdata, unsigned int tid) {
//   if (blockSize >= 64)
//     sdata[tid] += sdata[tid + 32];
//   if (blockSize >= 32)
//     sdata[tid] += sdata[tid + 16];
//   if (blockSize >= 16)
//     sdata[tid] += sdata[tid + 8];
//   if (blockSize >= 8)
//     sdata[tid] += sdata[tid + 4];
//   if (blockSize >= 4)
//     sdata[tid] += sdata[tid + 2];
//   if (blockSize >= 2)
//     sdata[tid] += sdata[tid + 1];
// }

// template <unsigned int blockSize>
// __global__ void reduction(double *g_idata, double *g_odata, unsigned int n) {
//   extern __shared__ double sdata[];
//   unsigned int tid = threadIdx.x;
//   unsigned int i = blockIdx.x * (blockSize * 2) + tid;
//   unsigned int gridSize = blockSize * 2 * gridDim.x;
//   sdata[tid] = 0;
//   while (i < n) {
//     sdata[tid] += g_idata[i] + g_idata[i + blockSize];
//     i += gridSize;
//   }
//   __syncthreads();
//   if (blockSize >= 512) {
//     if (tid < 256) {
//       sdata[tid] += sdata[tid + 256];
//     }
//     __syncthreads();
//   }
//   if (blockSize >= 256) {
//     if (tid < 128) {
//       sdata[tid] += sdata[tid + 128];
//     }
//     __syncthreads();
//   }
//   if (blockSize >= 128) {
//     if (tid < 64) {
//       sdata[tid] += sdata[tid + 64];
//     }
//     __syncthreads();
//   }
//   if (tid < 32)
//     warpReduce<THREADS_PER_BLOCK>(sdata, tid);
//   if (tid == 0)
//     g_odata[blockIdx.x] = sdata[0];
// }

__device__ void warpReduce(volatile double *sdata, int tid) {
  sdata[tid] += sdata[tid + 32];
  sdata[tid] += sdata[tid + 16];
  sdata[tid] += sdata[tid + 8];
  sdata[tid] += sdata[tid + 4];
  sdata[tid] += sdata[tid + 2];
  sdata[tid] += sdata[tid + 1];
}

template <unsigned int blockSize>
__global__ void reduction(double *g_idata, double *g_odata) {
  extern __shared__ double sdata[];
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * (blockDim.x * 2) + threadIdx.x;
  sdata[tid] = g_idata[i] + g_idata[i + blockDim.x];
  __syncthreads();
  for (unsigned int s = blockDim.x / 2; s > 32; s >>= 1) {
    if (tid < s)
      sdata[tid] += sdata[tid + s];
    __syncthreads();
  }
  if (tid < 32)
    warpReduce(sdata, tid);
  if (tid == 0)
    g_odata[blockIdx.x] = sdata[0];
}

// __global__ void reduction(double *g_idata, double *g_odata) {
//   extern __shared__ int sdata[];
//   unsigned int tid = threadIdx.x;
//   unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
//   sdata[tid] = g_idata[i];
//   __syncthreads();
//   for (unsigned int s = 1; s < blockDim.x; s *= 2) {
//     if (tid % (2 * s) == 0) {
//       sdata[tid] += sdata[tid + s];
//     }
//     __syncthreads();
//   }
//   if (tid == 0)
//     g_odata[blockIdx.x] = sdata[0];
// }

} // namespace reduction