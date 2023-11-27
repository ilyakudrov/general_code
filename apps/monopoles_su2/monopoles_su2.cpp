#include "../../lib/cpu/include/abelian_projection_su3.h"
#include "../../lib/cpu/include/basic_observables.h"
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
  bool convert = 0;

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
    } else if (string(argv[i]) == "-convert") {
      istringstream(string(argv[++i])) >> convert;
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
  cout << "convert " << convert << endl;

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

  data<abelian> conf;

  // read configuration
  get_data(conf, path_conf, conf_format, bytes_skip, convert);

  cout << plaket(conf.array) << endl;

  cout.precision(17);

  vector<double> angles = convert_to_angles(conf.array);

  conf.array.erase(conf.array.begin(), conf.array.end());

  std::ofstream output_stream_clusters_unwrapped(
      path_output_clusters_unwrapped);
  std::ofstream output_stream_clusters_wrapped(path_output_clusters_wrapped);
  std::ofstream output_stream_windings(path_output_windings);
  std::ofstream output_stream_monopoles(path_output_monopoles);

  output_stream_clusters_unwrapped << "length,number" << endl;
  output_stream_clusters_wrapped << "length,number,direction" << endl;
  output_stream_windings << "winding_number,cluster_number,direction" << endl;
  output_stream_monopoles << "asymmetry" << endl;

  vector<double> J = calculate_current(angles);

  vector<loop *> LL = calculate_clusters(J);

  int length;

  map<int, int> lengths_unwrapped;
  map<int, int> lengths_wrapped_time;
  map<int, int> lengths_wrapped_space;
  map<int, int> lengths_wrapped_both;
  map<int, int> space_windings;
  map<int, int> time_windings;
  vector<int> lengths_mu;
  vector<int> currents;

  int space_currents = 0;
  int time_currents = 0;

  vector<int> lattice_sizes = {x_size, y_size, z_size, t_size};

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
    output_stream_clusters_unwrapped << it->first << "," << it->second << endl;
  }

  for (auto it = lengths_wrapped_time.cbegin();
       it != lengths_wrapped_time.cend(); ++it) {
    output_stream_clusters_wrapped << it->first << "," << it->second << ",time"
                                   << endl;
  }
  for (auto it = lengths_wrapped_space.cbegin();
       it != lengths_wrapped_space.cend(); ++it) {
    output_stream_clusters_wrapped << it->first << "," << it->second << ",space"
                                   << endl;
  }
  for (auto it = lengths_wrapped_both.cbegin();
       it != lengths_wrapped_both.cend(); ++it) {
    output_stream_clusters_wrapped << it->first << "," << it->second << ",both"
                                   << endl;
  }

  for (auto it = time_windings.begin(); it != time_windings.end(); ++it) {
    output_stream_windings << it->first << "," << it->second << ",time" << endl;
  }

  for (auto it = space_windings.begin(); it != space_windings.end(); ++it) {
    output_stream_windings << it->first << "," << it->second << ",space"
                           << endl;
  }

  double asymmetry = (space_currents / 3. - time_currents) /
                     (space_currents / 3. + time_currents);

  output_stream_monopoles << asymmetry << endl;

  output_stream_clusters_unwrapped.close();
  output_stream_clusters_wrapped.close();
  output_stream_windings.close();
  output_stream_monopoles.close();
}
