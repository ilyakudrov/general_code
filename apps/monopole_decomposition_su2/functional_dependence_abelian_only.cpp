#include "../../lib/cpu/include/Landau_U1.h"
#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/decomposition.h"
#include "../../lib/cpu/include/mag.h"
#include "../../lib/cpu/include/matrix.h"
#include "../../lib/cpu/include/monopoles.h"
#include "../../lib/cpu/include/plaket.h"
#include "../../lib/cpu/include/smearing.h"
#include "../../lib/cpu/include/wilson_loops.h"

#include <filesystem>
#include <iostream>
#include <map>
#include <omp.h>

using namespace std;

// global variables for lattice size
int x_size;
int y_size;
int z_size;
int t_size;
int size1;
int size2;

void write_wilson_loops(map<tuple<int, int>, double> &wilson_tmp,
                        map<tuple<int, int, int, int>, double> &wilson_loops,
                        int smearing1, int smearing2) {
  for (auto it = wilson_tmp.begin(); it != wilson_tmp.end(); it++) {
    wilson_loops[tuple<int, int, int, int>(
        smearing1, smearing2, get<0>(it->first), get<1>(it->first))] +=
        it->second;
  }
}

int main(int argc, char **argv) {
  double omp_time;
  int x_size1;
  int y_size1;
  int z_size1;
  int t_size1;
  string path_conf;
  string conf_format;
  string file_precision;
  int bytes_skip = 0;
  string path_functional_output;
  string path_clusters_unwrapped_abelian_output;
  string path_clusters_wrapped_abelian_output;
  string path_windings_abelian_output;
  int copies_required;
  int mag_steps;
  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "--conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "--file_precision") {
      file_precision = argv[++i];
    } else if (string(argv[i]) == "--path_conf") {
      path_conf = argv[++i];
    } else if (string(argv[i]) == "--bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--copies_required") {
      copies_required = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--path_functional_output") {
      path_functional_output = argv[++i];
    } else if (string(argv[i]) == "--path_clusters_unwrapped_abelian_output") {
      path_clusters_unwrapped_abelian_output = argv[++i];
    } else if (string(argv[i]) == "--path_clusters_wrapped_abelian_output") {
      path_clusters_wrapped_abelian_output = argv[++i];
    } else if (string(argv[i]) == "--path_windings_abelian_output") {
      path_windings_abelian_output = argv[++i];
    } else if (string(argv[i]) == "--mag_steps") {
      mag_steps = stoi(string(argv[++i]));
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
  cout << "path_clusters_unwrapped_abelian_output "
       << path_clusters_unwrapped_abelian_output << endl;
  cout << "path_clusters_wrapped_abelian_output "
       << path_clusters_wrapped_abelian_output << endl;
  cout << "path_windings_abelian_output " << path_windings_abelian_output
       << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "copies_required " << copies_required << endl;
  cout << "mag_steps " << mag_steps << endl;
  cout << "x_size " << x_size1 << endl;
  cout << "y_size " << y_size1 << endl;
  cout << "z_size " << z_size1 << endl;
  cout << "t_size " << t_size1 << endl;

  // x_size = x_size1;
  // y_size = y_size1;
  // z_size = z_size1;
  // t_size = t_size1;

  string path_functional;
  string path_wilson_loops_monopole;
  string path_clusters_unwrapped_abelian;
  string path_clusters_wrapped_abelian;
  string path_windings_abelian;
  double MAG_functional;
  for (int copy = 0; copy < copies_required; copy++) {
    path_functional = path_functional_output + "_" + to_string(copy);
    path_clusters_unwrapped_abelian =
        path_clusters_unwrapped_abelian_output + "_" + to_string(copy);
    path_clusters_wrapped_abelian =
        path_clusters_wrapped_abelian_output + "_" + to_string(copy);
    path_windings_abelian =
        path_windings_abelian_output + "_" + to_string(copy);
    if (!(std::filesystem::exists(path_functional) &&
          std::filesystem::exists(path_wilson_loops_monopole) &&
          std::filesystem::exists(path_clusters_unwrapped_abelian) &&
          std::filesystem::exists(path_clusters_wrapped_abelian) &&
          std::filesystem::exists(path_windings_abelian))) {
      Data::LatticeData<DataPatternLexicographical, su2> conf_su2(
          {x_size1, y_size1, z_size1, t_size1});
      Data::read_data_convert(conf_su2, path_conf, conf_format, bytes_skip,
                              file_precision, 0);
      DataPatternLexicographical data_pattern(conf_su2.lat_dim);
      DataPatternLexicographical data_pattern_laplacian(
          {conf_su2.lat_dim[0] / 2 + 1, conf_su2.lat_dim[1] / 2 + 1,
           conf_su2.lat_dim[2] / 2 + 1, conf_su2.lat_dim[3] / 2 + 1});
      cout.precision(17);
      std::cout << "plaket before MAG: " << plaket(conf_su2) << std::endl;
      std::cout << "MAG functional before MAG: " << MAG_functional_su2(conf_su2)
                << std::endl;
      omp_time = omp_get_wtime();
      vector<spin> spins = generate_spins_uniform(data_pattern);
      if (mag_steps != 0) {
        make_simulated_annealing(conf_su2, spins, 2.5, 0.1, 2.4 / mag_steps, 6,
                                 50);
      }
      make_maximization_final(conf_su2, spins, 6, 1e-14, 1e-15);
      gauge_tranformation_spins(conf_su2, spins);
      std::cout << "MAG time: " << omp_get_wtime() - omp_time << std::endl;
      std::cout << "plaket after MAG: " << plaket(conf_su2) << std::endl;
      MAG_functional = MAG_functional_su2(conf_su2);
      std::cout << "MAG functional after MAG: " << MAG_functional << std::endl;

      Data::LatticeData<DataPatternLexicographical, abelian> conf_abelian(
          {x_size1, y_size1, z_size1, t_size1});
      conf_abelian.array = get_abelian(conf_su2);
      conf_su2.array.clear();
      conf_su2.array.shrink_to_fit();

      // monoples
      omp_time = omp_get_wtime();
      vector<double> J = calculate_current(conf_abelian);
      vector<loop_new *> LL = calculate_clusters(J, data_pattern);
      int length;
      vector<vector<int>> wrappings_abelian;
      vector<int> wrapped_lengths_abelian;
      map<int, int> lengths_unwrapped_abelian;
      map<int, int> lengths_wrapped_time;
      map<int, int> lengths_wrapped_space;
      map<int, int> lengths_wrapped_both;
      map<int, int> space_windings_abelian;
      map<int, int> time_windings_abelian;
      vector<int> lengths_mu;
      std::tuple<int, int> currents;
      int space_currents = 0;
      int time_currents = 0;
      for (int i = 0; i < LL.size(); i++) {
        length = cluster_length(LL[i]);
        lengths_mu = length_mu(LL[i]);
        currents = currents_directions(LL[i]);
        space_currents += std::get<0>(currents);
        time_currents += std::get<1>(currents);
        int wrappings_time_tmp = 0;
        int wrappings_space_tmp = 0;
        if (lengths_mu[0] == 0 && lengths_mu[1] == 0 && lengths_mu[2] == 0 &&
            lengths_mu[3] == 0) {
          lengths_unwrapped_abelian[length]++;
        } else {
          wrappings_abelian.push_back(
              vector<int>({lengths_mu[0] / conf_abelian.lat_dim[0],
                           lengths_mu[1] / conf_abelian.lat_dim[1],
                           lengths_mu[2] / conf_abelian.lat_dim[2],
                           lengths_mu[3] / conf_abelian.lat_dim[3]}));
          wrapped_lengths_abelian.push_back(length);
        }
      }
      // sort wrapped_lengths and wrappings accordingly
      std::vector<std::size_t> permutations(wrapped_lengths_abelian.size());
      std::iota(permutations.begin(), permutations.end(), 0);
      std::sort(permutations.begin(), permutations.end(),
                [&](std::size_t i, std::size_t j) {
                  return wrapped_lengths_abelian[i] >
                         wrapped_lengths_abelian[j];
                });
      std::vector<int> sorted_lengths(permutations.size());
      std::transform(permutations.begin(), permutations.end(),
                     sorted_lengths.begin(),
                     [&](std::size_t i) { return wrapped_lengths_abelian[i]; });
      wrapped_lengths_abelian = sorted_lengths;
      std::vector<std::vector<int>> sorted_wrappings(permutations.size(),
                                                     std::vector<int>());
      std::transform(permutations.begin(), permutations.end(),
                     sorted_wrappings.begin(),
                     [&](std::size_t i) { return wrappings_abelian[i]; });
      wrappings_abelian = sorted_wrappings;
      vector<int> positions_percolating_abelian;
      if (wrappings_abelian.size() > 1) {
        positions_percolating_abelian = group_percolating(wrappings_abelian);
      }
      std::cout << "monopoles time: " << omp_get_wtime() - omp_time
                << std::endl;

      std::ofstream output_stream_functional(path_functional);
      output_stream_functional << "functional" << endl;
      output_stream_functional << MAG_functional << endl;
      output_stream_functional.close();

      std::ofstream output_stream_clusters_unwrapped_abelian(
          path_clusters_unwrapped_abelian);
      output_stream_clusters_unwrapped_abelian << "length,number" << endl;
      for (auto it = lengths_unwrapped_abelian.cbegin();
           it != lengths_unwrapped_abelian.cend(); ++it) {
        output_stream_clusters_unwrapped_abelian << it->first << ","
                                                 << it->second << endl;
      }
      output_stream_clusters_unwrapped_abelian.close();
      std::ofstream output_stream_clusters_wrapped_abelian(
          path_clusters_wrapped_abelian);
      output_stream_clusters_wrapped_abelian
          << "length,x0_wrap,x1_wrap,x2_wrap,x3_wrap,percolating_group" << endl;
      for (int i = 0; i < wrapped_lengths_abelian.size(); i++) {
        if (std::find(positions_percolating_abelian.begin(),
                      positions_percolating_abelian.end(),
                      i) != positions_percolating_abelian.end()) {
          output_stream_clusters_wrapped_abelian
              << wrapped_lengths_abelian[i] << "," << wrappings_abelian[i][0]
              << "," << wrappings_abelian[i][1] << ","
              << wrappings_abelian[i][2] << "," << wrappings_abelian[i][3]
              << ",percolating" << endl;
        } else {
          output_stream_clusters_wrapped_abelian
              << wrapped_lengths_abelian[i] << "," << wrappings_abelian[i][0]
              << "," << wrappings_abelian[i][1] << ","
              << wrappings_abelian[i][2] << "," << wrappings_abelian[i][3]
              << ",non-percolating" << endl;
        }
      }
      output_stream_clusters_wrapped_abelian.close();
      std::ofstream output_stream_windings_abelian(path_windings_abelian);
      output_stream_windings_abelian
          << "winding_number,cluster_number,direction" << endl;
      for (auto it = time_windings_abelian.begin();
           it != time_windings_abelian.end(); ++it) {
        output_stream_windings_abelian << it->first << "," << it->second
                                       << ",time" << endl;
      }
      for (auto it = space_windings_abelian.begin();
           it != space_windings_abelian.end(); ++it) {
        output_stream_windings_abelian << it->first << "," << it->second
                                       << ",space" << endl;
      }
      output_stream_windings_abelian.close();
    }
  }
}