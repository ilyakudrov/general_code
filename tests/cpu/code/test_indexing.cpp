#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/indexing.h"
#include "../../../lib/cpu/include/link.h"
#include "../../../lib/cpu/include/matrix.h"
#include "../../../lib/cpu/include/polyakov_loops.h"
#include "../../../lib/cpu/include/smearing.h"

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <omp.h>

int x_size;
int y_size;
int z_size;
int t_size;
int size1;
int size2;

#define MATRIX_TYPE su3

using namespace std;

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double search_time;

  x_size = 24;
  y_size = 24;
  z_size = 24;
  t_size = 24;
  size1 = x_size * y_size;
  size2 = x_size * y_size * z_size;

  std::cout.precision(17);

  Data::data<MATRIX_TYPE> conf1;
  Data::data<MATRIX_TYPE> conf2;

  // string conf_path1 = "../../confs/su3/QCD/140MeV/nt4/conf.0501";
  // string conf_format1 = "ildg";
  string conf_path1 = "../../confs/su3/gluodynamics/24^4/beta6.0/CONF0001";
  string conf_format1 = "double_qc2dstag";
  // string conf_path1 =
  // "../../confs/MAG/su3/gluodynamics/40^4/beta6.4/steps_0/"
  //                     "copies=20/s1/conf_gaugefixed_0002_1";
  // string conf_format1 = "double";
  // string conf_path1 = "../../confs/su2/qc2dstag/40^4/mu0.00/CONF0201";
  // string conf_format1 = "double_qc2dstag";
  int bytes_skip = 0;
  bool convert = 0;

  int R_min = 1;
  int R_max = 4;
  int T_min = 1;
  int T_max = 4;

  map<tuple<int, int>, double> wilson_loops;

  get_data(conf1, conf_path1, conf_format1, bytes_skip, convert);

  start_time = omp_get_wtime();
  std::cout << polyakov_loop(conf1.array) << std::endl;
  std::cout << plaket(conf1.array) << std::endl;
  std::cout << plaket_space(conf1.array) << std::endl;
  std::cout << plaket_time(conf1.array) << std::endl;
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "observables time: " << search_time << std::endl;

  conf1.array.clear();
  conf1.array.shrink_to_fit();

  std::array<int, 4> lat_dim = {x_size, y_size, z_size, t_size};
  DataPatternLexicographical<4> data_pattern(lat_dim);
  Data::Data1<DataPatternLexicographical<4>, su3> data_indexed(
      (DataPatternLexicographical<4>(lat_dim)));
  FilePatternQCDSTAG<4> file_pattern;
  data_indexed.read_data(conf_path1, file_pattern, 0, "double");

  start_time = omp_get_wtime();
  std::cout << polyakov_loop(data_indexed.array) << std::endl;
  std::cout << plaket(data_indexed.array) << std::endl;
  std::cout << plaket_space(data_indexed.array) << std::endl;
  std::cout << plaket_time(data_indexed.array) << std::endl;
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "observables data_indexed time: " << search_time << std::endl;

  FilePatternLexicographical<4> file_pattern_lexicographical;
  data_indexed.write_data("../../confs/test/conf_test",
                          file_pattern_lexicographical);

  data_pattern = DataPatternLexicographical<4>(lat_dim);
  data_indexed.read_data("../../confs/test/conf_test",
                         file_pattern_lexicographical, 0, "double");

  start_time = omp_get_wtime();
  std::cout << polyakov_loop(data_indexed.array) << std::endl;
  std::cout << plaket(data_indexed.array) << std::endl;
  std::cout << plaket_space(data_indexed.array) << std::endl;
  std::cout << plaket_time(data_indexed.array) << std::endl;
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  std::cout << "observables file pattern lexicographical time: " << search_time
            << std::endl;
}