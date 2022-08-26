#include "../../lib/cpu/include/abelian_projection_su3.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/matrix.h"
#include "../../lib/cpu/include/monopoles.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>

using namespace std;

// global variables for lattice size
int x_size;
int y_size;
int z_size;
int t_size;

int main(int argc, char **argv) {

  unsigned int start_time;
  unsigned int end_time;
  unsigned int search_time;

  string path_conf;
  string conf_format;
  int bytes_skip = 0;
  string path_output_clusters;
  string path_output_windings;
  string path_output_monopoles;

  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "-path_conf") {
      path_conf = argv[++i];
    } else if (string(argv[i]) == "-path_output_clusters") {
      path_output_clusters = argv[++i];
    } else if (string(argv[i]) == "-path_output_windings") {
      path_output_windings = argv[++i];
    } else if (string(argv[i]) == "-path_output_monopoles") {
      path_output_monopoles = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-x_size") {
      x_size = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-y_size") {
      y_size = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-z_size") {
      z_size = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-t_size") {
      t_size = stoi(string(argv[++i]));
    } else
      cout << "unknown parameter " << argv[i] << endl;
  }

  cout << "path_conf " << path_conf << endl;
  cout << "conf_format " << conf_format << endl;
  cout << "bytes_skip " << bytes_skip << endl;

  cout << "path_output_clusters " << path_output_clusters << endl;
  cout << "path_output_windings " << path_output_windings << endl;
  cout << "path_output_monopoles " << path_output_monopoles << endl;

  cout << "x_size " << x_size << endl;
  cout << "y_size " << y_size << endl;
  cout << "z_size " << z_size << endl;
  cout << "t_size " << t_size << endl;

  data<su3> conf_su3;

  // read configuration
  if (std::string(conf_format) == "float") {
    conf_su3.read_float(path_conf, bytes_skip);
  } else if (std::string(conf_format) == "double") {
    conf_su3.read_double(path_conf, bytes_skip);
  } else if (std::string(conf_format) == "double_qc2dstag") {
    conf_su3.read_double_qc2dstag(path_conf);
  } else if (std::string(conf_format) == "ildg") {
    conf_su3.read_ildg(path_conf);
  } else {
    cout << "wrong conf format: " << conf_format << endl;
    return 0;
  }

  cout.precision(17);

  vector<vector<double>> angles = make_angles_SU3(conf_su3.array);
  conf_su3.array.erase(conf_su3.array.begin(), conf_su3.array.end());

  std::ofstream output_stream_clusters(path_output_clusters);
  std::ofstream output_stream_windings(path_output_windings);
  std::ofstream output_stream_monopoles(path_output_monopoles);

  output_stream_clusters << "color,length,number" << endl;
  output_stream_windings << "color,length,number,type" << endl;
  output_stream_monopoles << "" << endl;

  for (int color = 0; color < angles.size(); color++) {

    vector<double> J = calculate_current(angles[color]);
    vector<loop *> LL = calculate_clusters(J);

    int length;

    map<int, int> lengths;
    map<int, int> space_windings;
    map<int, int> time_windings;
    vector<int> lengths_mu;
    int length_mu_test;
    vector<int> currents;

    for (int i = 0; i < LL.size(); i++) {
      length = cluster_length(LL[i]);
      lengths[length]++;
      lengths_mu = length_mu(LL[i]);

      for (int j = 0; j < 3; j++) {
        if (lengths_mu[j] != 0) {
          space_windings[abs(lengths_mu[j])]++;
        }
      }

      if (lengths_mu[3] != 0) {
        time_windings[abs(lengths_mu[3])]++;
      }
    }

    for (auto it = lengths.cbegin(); it != lengths.cend(); ++it) {
      output_stream_clusters << color + 1 << "," << it->first << ","
                             << it->second << endl;
    }

    for (auto it = space_windings.begin(); it != space_windings.end(); ++it) {
      output_stream_windings << color + 1 << "," << it->first << " "
                             << it->second << ",space" << endl;
    }

    for (auto it = time_windings.begin(); it != time_windings.end(); ++it) {
      output_stream_windings << color + 1 << "," << it->first << " "
                             << it->second << ",time" << endl;
    }
  }

  output_stream_clusters.close();
  output_stream_windings.close();
  output_stream_monopoles.close();
}
