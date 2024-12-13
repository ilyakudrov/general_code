#include "../../../lib/cpu/include/Landau_U1.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/decomposition.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/monopoles.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <vector>

int x_size;
int y_size;
int z_size;
int t_size;

using namespace std;

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double search_time;

  x_size = 24;
  y_size = 24;
  z_size = 24;
  t_size = 24;

  double tolerance_maximal = 1e-5;
  double tolerance_average = 1e-7;
  int OR_steps = 4;

  std::cout.precision(17);

  link1 link(x_size, y_size, z_size, t_size);
  std::string path_conf =
      "../../confs/su2/mag/su2_suzuki/24^4/beta2.4/conf_0001";
  // std::string path_conf =
  //     "../../confs/su2/mag/qc2dstag/40^4/mu0.00/conf_abelian_0201";

  data<su2> conf_su2;
  conf_su2.read_double(path_conf, 0);

  string laplacian_path = "../../confs/su2/inverse_laplacian/ALPHA24x24_d.LAT";

  vector<double> inverse_laplacian = read_inverse_laplacian(laplacian_path);

  vector<complex_t> conf_complex = convert_to_complex(conf_su2.array);
  vector<complex_t> gauge_complex = generate_gauge_complex_unity();

  start_time = omp_get_wtime();

  make_maximization_final(gauge_complex, conf_complex, OR_steps,
                          tolerance_maximal, tolerance_average);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "Landau gauge time: " << search_time << std::endl;

  normalize_complex(gauge_complex);

  apply_gauge_Landau_complex(gauge_complex, conf_complex);

  vector<double> conf_angles_U1 = convert_complex_to_angles(conf_complex);

  vector<vector<int>> dirac_plakets =
      calculate_monopole_plaket_singular(conf_angles_U1);

  vector<double> monopole_angles;

  start_time = omp_get_wtime();

  monopole_angles = make_monopole_angles(dirac_plakets, inverse_laplacian);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "make_monopole_angles time: " << search_time << std::endl;

  double cos_sum = 0;
  for (int i = 0; i < monopole_angles.size(); i++) {
    cos_sum += cos(monopole_angles[i]);
  }
  cos_sum = cos_sum / monopole_angles.size();

  cout << "cos sum = " << cos_sum << endl;

  start_time = omp_get_wtime();

  monopole_angles =
      make_monopole_angles_parallel(dirac_plakets, inverse_laplacian);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "make_monopole_angles_parallel time: " << search_time
            << std::endl;

  cos_sum = 0;
  for (int i = 0; i < monopole_angles.size(); i++) {
    cos_sum += cos(monopole_angles[i]);
  }
  cos_sum = cos_sum / monopole_angles.size();

  cout << "cos sum = " << cos_sum << endl;
}