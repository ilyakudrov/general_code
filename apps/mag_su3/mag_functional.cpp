#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/mag.h"
#include "../../lib/cpu/include/matrix.h"

#include <fstream>
#include <iostream>

using namespace std;

// global variables for lattice size
int x_size;
int y_size;
int z_size;
int t_size;
int size1;
int size2;

int main(int argc, char **argv) {
  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  int x_size1;
  int y_size1;
  int z_size1;
  int t_size1;

  string path_conf;
  string conf_format;
  string file_precision;
  int bytes_skip = 0;
  string path_functional_output;

  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "--conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "--file_precision") {
      file_precision = argv[++i];
    } else if (string(argv[i]) == "--path_conf") {
      path_conf = argv[++i];
    } else if (string(argv[i]) == "--path_functional_output") {
      path_functional_output = argv[++i];
    } else if (string(argv[i]) == "--bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--x_size") {
      x_size1 = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--y_size") {
      y_size1 = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--z_size") {
      z_size1 = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--t_size") {
      t_size1 = stoi(string(argv[++i]));
    } else
      cout << "unknown parameter " << argv[i] << endl;
  }

  cout << "path_conf " << path_conf << endl;
  cout << "conf_format " << conf_format << endl;
  cout << "file_precision " << file_precision << endl;
  cout << "path_functional_output " << path_functional_output << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "x_size " << x_size1 << endl;
  cout << "y_size " << y_size1 << endl;
  cout << "z_size " << z_size1 << endl;
  cout << "t_size " << t_size1 << endl;

  bool convert = 0;
  Data::LatticeData<DataPatternLexicographical, su3> conf_su3(
      {x_size1, y_size1, z_size1, t_size1});
  Data::read_data_convert(conf_su3, path_conf, conf_format, bytes_skip,
                          file_precision, convert);

  cout.precision(17);

  std::ofstream output_functional(path_functional_output);
  output_functional.precision(17);
  output_functional << "functional" << endl;
  output_functional << mag_functional_su3(conf_su3.array) << endl;
  output_functional.close();
}
