#pragma once

#if defined(__CUDA_ARCH__) || defined(__CUDA__)
#define CUDA_HOST_DEVICE __device__ __host__
#define CUDA_DEVICE __device__
#define CUDA_HOST __host__
#else
#define CUDA_HOST_DEVICE
#define CUDA_DEVICE
#define CUDA_HOST
#endif