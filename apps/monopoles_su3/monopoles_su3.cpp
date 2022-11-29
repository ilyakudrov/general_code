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
  string path_output_clusters_unwrapped;
  string path_output_clusters_wrapped;
  string path_output_windings;
  string path_output_monopoles;

  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "-path_conf") {
      path_conf = argv[++i];
    } else if (string(argv[i]) == "-path_output_clusters_unwrapped") {
      path_output_clusters_unwrapped = argv[++i];
    } else if (string(argv[i]) == "-path_output_clusters_wrapped") {
      path_output_clusters_wrapped = argv[++i];
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

  cout << "path_output_clusters_unwrapped " << path_output_clusters_unwrapped
       << endl;
  cout << "path_output_clusters_wrapped " << path_output_clusters_wrapped
       << endl;
  cout << "path_output_windings " << path_output_windings << endl;
  cout << "path_output_monopoles " << path_output_monopoles << endl;

  cout << "x_size " << x_size << endl;
  cout << "y_size " << y_size << endl;
  cout << "z_size " << z_size << endl;
  cout << "t_size " << t_size << endl;

  vector<int> lattice_sizes = {x_size, y_size, z_size, t_size};

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

  vector<vector<vector<double>>> monopole_plakets =
      make_monopole_plakets(angles);

  ofstream output_stream_clusters_unwrapped(path_output_clusters_unwrapped);
  ofstream output_stream_clusters_wrapped(path_output_clusters_wrapped);
  ofstream output_stream_windings(path_output_windings);
  ofstream output_stream_monopoles(path_output_monopoles);

  output_stream_clusters_unwrapped << "color,length,number" << endl;
  output_stream_clusters_wrapped << "color,length,number,direction" << endl;
  output_stream_windings << "color,winding_number,cluster_number,direction"
                         << endl;
  output_stream_monopoles << "color,asymmetry" << endl;

  for (int color = 0; color < angles.size(); color++) {

    vector<double> J =
        calculate_current_monopole_plakets(monopole_plakets[color]);
    vector<loop *> LL = calculate_clusters(J);

    int length;

    map<int, int> lengths_unwrapped;
    map<int, int> lengths_wrapped_time;
    map<int, int> lengths_wrapped_space;
    map<int, int> lengths_wrapped_both;
    map<int, int> space_windings;
    map<int, int> time_windings;
    vector<int> lengths_mu;
    int length_mu_test;
    vector<int> currents;

    int space_currents = 0;
    int time_currents = 0;

    for (int i = 0; i < LL.size(); i++) {
      length = cluster_length(LL[i]);
      lengths_mu = length_mu(LL[i]);

      currents = currents_directions(LL[i]);

      space_currents += currents[0];
      time_currents += currents[1];

      int wrappings_time_tmp = 0;
      int wrappings_space_tmp = 0;

      if (lengths_mu[0] == 0 && lengths_mu[1] == 0 && lengths_mu[2] == 0 &&
          lengths_mu[3] == 0) {
        lengths_unwrapped[length]++;
      }

      else {

        for (int mu = 0; mu < 3; mu++) {
          if (lengths_mu[mu] != 0) {
            wrappings_space_tmp += abs(lengths_mu[mu]) / lattice_sizes[mu];
          }
        }

        if (lengths_mu[3] != 0) {
          wrappings_time_tmp += abs(lengths_mu[3]) / lattice_sizes[3];
        }

        if (wrappings_space_tmp != 0) {
          space_windings[wrappings_space_tmp]++;
        }

        if (wrappings_time_tmp != 0) {
          time_windings[wrappings_time_tmp]++;
        }

        if (wrappings_space_tmp != 0 && wrappings_time_tmp != 0) {
          lengths_wrapped_both[length]++;
        } else if (wrappings_space_tmp != 0) {
          lengths_wrapped_space[length]++;
        } else if (wrappings_time_tmp != 0) {
          lengths_wrapped_time[length]++;
        }
      }
    }

    for (auto it = lengths_unwrapped.cbegin(); it != lengths_unwrapped.cend();
         ++it) {
      output_stream_clusters_unwrapped << color + 1 << "," << it->first << ","
                                       << it->second << endl;
    }

    for (auto it = lengths_wrapped_time.cbegin();
         it != lengths_wrapped_time.cend(); ++it) {
      output_stream_clusters_wrapped << color + 1 << "," << it->first << ","
                                     << it->second << ",time" << endl;
    }
    for (auto it = lengths_wrapped_space.cbegin();
         it != lengths_wrapped_space.cend(); ++it) {
      output_stream_clusters_wrapped << color + 1 << "," << it->first << ","
                                     << it->second << ",space" << endl;
    }
    for (auto it = lengths_wrapped_both.cbegin();
         it != lengths_wrapped_both.cend(); ++it) {
      output_stream_clusters_wrapped << color + 1 << "," << it->first << ","
                                     << it->second << ",both" << endl;
    }

    for (auto it = time_windings.begin(); it != time_windings.end(); ++it) {
      output_stream_windings << color + 1 << "," << it->first << ","
                             << it->second << ",time" << endl;
    }

    for (auto it = space_windings.begin(); it != space_windings.end(); ++it) {
      output_stream_windings << color + 1 << "," << it->first << ","
                             << it->second << ",space" << endl;
    }

    double asymmetry = (space_currents / 3. - time_currents) /
                       (space_currents / 3. + time_currents);

    output_stream_monopoles << color + 1 << "," << asymmetry << endl;
  }

  output_stream_clusters_unwrapped.close();
  output_stream_clusters_wrapped.close();
  output_stream_windings.close();
  output_stream_monopoles.close();
}
