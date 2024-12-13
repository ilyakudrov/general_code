#include "../../../lib/cpu/include/data.h"
// #include "../../../lib/gpu/include/general/reduction.h"
#include "../../../lib/gpu/include/observables/plaket.h"

// #include "cuda_profiler_api.h"

#define gpuErrchk(ans)                                                         \
  { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line,
                      bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
            line);
    if (abort)
      exit(code);
  }
}

int x_size;
int y_size;
int z_size;
int t_size;

#define MATRIX_TYPE_CPU su3
#define MATRIX_TYPE_GPU matrix_gpu::su3

using namespace std;

#define NX 36
#define NY 36
#define NZ 36
#define NT 36
#define THREADS_PER_BLOCK 128
#define BLOCK_NUMBER NX *NY *NZ *NT / THREADS_PER_BLOCK

int main() {

  x_size = NX;
  y_size = NY;
  z_size = NZ;
  t_size = NT;

  int lattice_size = x_size * y_size * z_size * t_size;

  std::cout.precision(17);

  data<MATRIX_TYPE_CPU> conf;

  string conf_path = "../../confs/SU3_conf/gluodynamics/36^4/beta6.3/CONF0001";

  string conf_format = "double_qc2dstag";
  int bytes_skip = 0;
  bool convert = 0;

  get_data(conf, conf_path, conf_format, bytes_skip, convert);

  MATRIX_TYPE_GPU *dconf;
  gpuErrchk(cudaMalloc(&dconf, lattice_size * 4 * sizeof(MATRIX_TYPE_GPU)));
  gpuErrchk(cudaMemcpy(dconf, &conf.array[0],
                       lattice_size * 4 * sizeof(MATRIX_TYPE_GPU),
                       cudaMemcpyHostToDevice));

  unsigned int *lat_indices =
      (unsigned int *)malloc(lattice_size * 9 * sizeof(unsigned int));
  int lattice_sizes[4] = {NX, NY, NZ, NT};
  plaket_gpu::make_plaket_indices(lat_indices, lattice_sizes);

  unsigned int *dlat_indices;
  gpuErrchk(cudaMalloc(&dlat_indices, lattice_size * 9 * sizeof(unsigned int)));
  gpuErrchk(cudaMemcpy(dlat_indices, lat_indices,
                       lattice_size * 9 * sizeof(unsigned int),
                       cudaMemcpyHostToDevice));

  double *dtraces;
  gpuErrchk(cudaMalloc(&dtraces, lattice_size * sizeof(double)));

  dim3 dimBlock(4, THREADS_PER_BLOCK / 4, 1);

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // plaket
  cudaEventRecord(start);

  // cudaProfilerStart();

  plaket_gpu::plaket<<<BLOCK_NUMBER, THREADS_PER_BLOCK>>>(dconf, dtraces,
                                                          dlat_indices);

  // cudaProfilerStop();

  cudaEventRecord(stop);

  gpuErrchk(cudaEventSynchronize(stop));
  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);

  gpuErrchk(cudaPeekAtLastError());
  gpuErrchk(cudaDeviceSynchronize());

  cout << "kernal plaket time " << milliseconds / 1000 << endl;

  double *traces = (double *)malloc(lattice_size * sizeof(double));
  gpuErrchk(cudaMemcpy(traces, dtraces, lattice_size * sizeof(double),
                       cudaMemcpyDeviceToHost));

  double plaket = 0;
  for (int i = 0; i < lattice_size; i++) {
    plaket += traces[i];
  }

  cout << "plaket aver " << plaket / lattice_size / 6 << endl;

  gpuErrchk(cudaFree(dtraces));
  free(traces);

  // plaket_test1
  gpuErrchk(cudaMalloc(&dtraces, lattice_size * 6 * sizeof(double)));
  cudaEventRecord(start);

  plaket_gpu::plaket_test1<<<6 * BLOCK_NUMBER, THREADS_PER_BLOCK>>>(
      dconf, dtraces, dlat_indices);

  cudaEventRecord(stop);

  gpuErrchk(cudaEventSynchronize(stop));
  milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);

  gpuErrchk(cudaPeekAtLastError());
  gpuErrchk(cudaDeviceSynchronize());

  cout << "kernal plaket_test1 time " << milliseconds / 1000 << endl;

  traces = (double *)malloc(lattice_size * 6 * sizeof(double));
  gpuErrchk(cudaMemcpy(traces, dtraces, lattice_size * 6 * sizeof(double),
                       cudaMemcpyDeviceToHost));

  plaket = 0;
  for (int i = 0; i < lattice_size * 6; i++) {
    plaket += traces[i];
  }

  cout << "plaket_test1 aver " << plaket / lattice_size / 6 << endl;

  gpuErrchk(cudaFree(dtraces));
  free(traces);

  // plaket_test4

  gpuErrchk(cudaMalloc(&dtraces, lattice_size * 6 * sizeof(double)));
  cudaEventRecord(start);

  dim3 dimGrid(6, 16);
  cout << "number of blocks " << (NX * NY * NZ * NT * 6 + 6 * 16 - 1) / (6 * 16)
       << endl;

  plaket_gpu::plaket_test4<<<(NX * NY * NZ * NT * 6 + 6 * 16 - 1) / (6 * 16),
                             dimGrid>>>(dconf, dtraces, dlat_indices);

  cudaEventRecord(stop);

  gpuErrchk(cudaEventSynchronize(stop));
  milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);

  gpuErrchk(cudaPeekAtLastError());
  gpuErrchk(cudaDeviceSynchronize());

  cout << "kernal plaket_test4 time " << milliseconds / 1000 << endl;

  traces = (double *)malloc(lattice_size * 6 * sizeof(double));
  gpuErrchk(cudaMemcpy(traces, dtraces, lattice_size * 6 * sizeof(double),
                       cudaMemcpyDeviceToHost));

  plaket = 0;
  for (int i = 0; i < lattice_size * 6; i++) {
    plaket += traces[i];
  }

  cout << "plaket_test4 aver " << plaket / lattice_size / 6 << endl;

  gpuErrchk(cudaFree(dlat_indices));
  free(lat_indices);

  // plaket_test2

  lat_indices =
      (unsigned int *)malloc(lattice_size * 24 * sizeof(unsigned int));
  plaket_gpu::make_plaket_indices_test2(lat_indices, lattice_sizes);

  gpuErrchk(
      cudaMalloc(&dlat_indices, lattice_size * 24 * sizeof(unsigned int)));
  gpuErrchk(cudaMemcpy(dlat_indices, lat_indices,
                       lattice_size * 24 * sizeof(unsigned int),
                       cudaMemcpyHostToDevice));

  cudaEventRecord(start);

  plaket_gpu::plaket_test2<<<6 * BLOCK_NUMBER, THREADS_PER_BLOCK>>>(
      dconf, dtraces, dlat_indices);

  cudaEventRecord(stop);

  gpuErrchk(cudaEventSynchronize(stop));
  milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);

  gpuErrchk(cudaPeekAtLastError());
  gpuErrchk(cudaDeviceSynchronize());

  cout << "kernal plaket_test2 time " << milliseconds / 1000 << endl;

  gpuErrchk(cudaMemcpy(traces, dtraces, lattice_size * 6 * sizeof(double),
                       cudaMemcpyDeviceToHost));

  plaket = 0;
  for (int i = 0; i < lattice_size * 6; i++) {
    plaket += traces[i];
  }

  cout << "plaket_test2 aver " << plaket / lattice_size / 6 << endl;

  free(lat_indices);
  gpuErrchk(cudaFree(dlat_indices));
  gpuErrchk(cudaFree(dconf));

  // plaket_test3

  unsigned int hsteps[4] = {4, 4 * NX, 4 * NX * NY, 4 * NX * NY * NZ};
  unsigned int harr_mu[6] = {0, 0, 0, 1, 1, 2};
  unsigned int harr_nu[6] = {1, 2, 3, 2, 3, 3};

  gpuErrchk(cudaMemcpyToSymbol(steps, hsteps, 4 * sizeof(unsigned int)));
  gpuErrchk(cudaMemcpyToSymbol(arr_mu, harr_mu, 6 * sizeof(unsigned int)));
  gpuErrchk(cudaMemcpyToSymbol(arr_nu, harr_nu, 6 * sizeof(unsigned int)));

  gpuErrchk(cudaMalloc(&dconf, lattice_size * 4 * 2 * sizeof(MATRIX_TYPE_GPU)));
  conf.array.reserve(lattice_size * 4 * 2);
  for (int i = 0; i < lattice_size * 4; i++) {
    conf.array.push_back(conf.array[i]);
  }
  gpuErrchk(cudaMemcpy(dconf, &conf.array[0],
                       lattice_size * 4 * 2 * sizeof(MATRIX_TYPE_GPU),
                       cudaMemcpyHostToDevice));

  cudaEventRecord(start);

  plaket_gpu::plaket_test3<<<6 * BLOCK_NUMBER, THREADS_PER_BLOCK>>>(dconf,
                                                                    dtraces);

  cudaEventRecord(stop);

  gpuErrchk(cudaEventSynchronize(stop));
  milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);

  gpuErrchk(cudaPeekAtLastError());
  gpuErrchk(cudaDeviceSynchronize());

  cout << "kernal plaket_test3 time " << milliseconds / 1000 << endl;

  gpuErrchk(cudaMemcpy(traces, dtraces, lattice_size * 6 * sizeof(double),
                       cudaMemcpyDeviceToHost));

  plaket = 0;
  for (int i = 0; i < lattice_size * 6; i++) {
    plaket += traces[i];
  }

  cout << "plaket_test3 aver " << plaket / lattice_size / 6 << endl;

  gpuErrchk(cudaFree(dtraces));
  free(traces);

  gpuErrchk(cudaFree(dconf));

  // double *d_odata;
  // gpuErrchk(cudaMalloc(&d_odata, BLOCK_NUMBER * sizeof(double)));

  // reduction::reduction<THREADS_PER_BLOCK>
  //     <<<BLOCK_NUMBER / 2, THREADS_PER_BLOCK>>>(dtraces, d_odata);

  // gpuErrchk(cudaPeekAtLastError());
  // gpuErrchk(cudaDeviceSynchronize());

  // double *odata = (double *)malloc(BLOCK_NUMBER * sizeof(double));
  // gpuErrchk(cudaMemcpy((void *)odata, (const void *)d_odata, sizeof(double),
  //                      cudaMemcpyDeviceToHost));
}