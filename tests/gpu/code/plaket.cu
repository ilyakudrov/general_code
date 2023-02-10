#include "../../../lib/cpu/include/data.h"
#include "../../../lib/gpu/include/general/reduction.h"
#include "../../../lib/gpu/include/observables/plaket.h"

#include "cuda_profiler_api.h"

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

  x_size = 36;
  y_size = 36;
  z_size = 36;
  t_size = 36;

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

  int *lat_indices = (int *)malloc(lattice_size * 9 * sizeof(int));
  int lattice_sizes[4] = {NX, NY, NZ, NT};
  plaket_gpu::make_plaket_indices(lat_indices, lattice_sizes);

  int *dlat_indices;
  gpuErrchk(cudaMalloc(&dlat_indices, lattice_size * 9 * sizeof(int)));
  gpuErrchk(cudaMemcpy(dlat_indices, lat_indices,
                       lattice_size * 9 * sizeof(int), cudaMemcpyHostToDevice));

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

  cout << "kernal plaket time " << milliseconds / 1000 << endl;

  traces = (double *)malloc(lattice_size * 6 * sizeof(double));
  gpuErrchk(cudaMemcpy(traces, dtraces, lattice_size * 6 * sizeof(double),
                       cudaMemcpyDeviceToHost));

  plaket = 0;
  for (int i = 0; i < lattice_size * 6; i++) {
    plaket += traces[i];
  }

  cout << "plaket aver " << plaket / lattice_size / 6 << endl;

  gpuErrchk(cudaFree(dtraces));
  free(traces);

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