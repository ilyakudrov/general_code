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
  string path_wilson_loops_abelian_output;
  string path_wilson_loops_monopole_output;
  string path_clusters_unwrapped_abelian_output;
  string path_clusters_wrapped_abelian_output;
  string path_windings_abelian_output;
  string path_clusters_unwrapped_monopole_output;
  string path_clusters_wrapped_monopole_output;
  string path_windings_monopole_output;
  string path_clusters_unwrapped_monopoless_output;
  string path_clusters_wrapped_monopoless_output;
  string path_windings_monopoless_output;
  string path_inverse_laplacian;
  int N_dir_gevp;
  double HYP_alpha1, HYP_alpha2, HYP_alpha3;
  double APE_alpha;
  int HYP_steps, APE_steps;
  int calculation_APE_start, calculation_step_APE;
  int copies_required;
  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "--conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "--file_precision") {
      file_precision = argv[++i];
    } else if (string(argv[i]) == "--path_conf") {
      path_conf = argv[++i];
    } else if (string(argv[i]) == "--HYP_alpha1") {
      HYP_alpha1 = atof(argv[++i]);
    } else if (string(argv[i]) == "--HYP_alpha2") {
      HYP_alpha2 = atof(argv[++i]);
    } else if (string(argv[i]) == "--HYP_alpha3") {
      HYP_alpha3 = atof(argv[++i]);
    } else if (string(argv[i]) == "--APE_alpha") {
      APE_alpha = atof(argv[++i]);
    } else if (string(argv[i]) == "--APE_steps") {
      APE_steps = stoi(argv[++i]);
    } else if (string(argv[i]) == "--HYP_steps") {
      HYP_steps = stoi(argv[++i]);
    } else if (string(argv[i]) == "--calculation_step_APE") {
      calculation_step_APE = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--calculation_APE_start") {
      calculation_APE_start = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--N_dir_gevp") {
      N_dir_gevp = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--path_inverse_laplacian") {
      path_inverse_laplacian = argv[++i];
    } else if (string(argv[i]) == "--copies_required") {
      copies_required = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "--path_functional_output") {
      path_functional_output = argv[++i];
    } else if (string(argv[i]) == "--path_wilson_loops_abelian_output") {
      path_wilson_loops_abelian_output = argv[++i];
    } else if (string(argv[i]) == "--path_wilson_loops_monopole_output") {
      path_wilson_loops_monopole_output = argv[++i];
    } else if (string(argv[i]) == "--path_clusters_unwrapped_abelian_output") {
      path_clusters_unwrapped_abelian_output = argv[++i];
    } else if (string(argv[i]) == "--path_clusters_wrapped_abelian_output") {
      path_clusters_wrapped_abelian_output = argv[++i];
    } else if (string(argv[i]) == "--path_windings_abelian_output") {
      path_windings_abelian_output = argv[++i];
    } else if (string(argv[i]) == "--path_clusters_unwrapped_monopole_output") {
      path_clusters_unwrapped_monopole_output = argv[++i];
    } else if (string(argv[i]) == "--path_clusters_wrapped_monopole_output") {
      path_clusters_wrapped_monopole_output = argv[++i];
    } else if (string(argv[i]) == "--path_windings_monopole_output") {
      path_windings_monopole_output = argv[++i];
    } else if (string(argv[i]) ==
               "--path_clusters_unwrapped_monopoless_output") {
      path_clusters_unwrapped_monopoless_output = argv[++i];
    } else if (string(argv[i]) == "--path_clusters_wrapped_monopoless_output") {
      path_clusters_wrapped_monopoless_output = argv[++i];
    } else if (string(argv[i]) == "--path_windings_monopoless_output") {
      path_windings_monopoless_output = argv[++i];
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
  cout << "path_wilson_loops_abelian_output "
       << path_wilson_loops_abelian_output << endl;
  cout << "path_wilson_loops_monopole_output "
       << path_wilson_loops_monopole_output << endl;
  cout << "path_clusters_unwrapped_abelian_output "
       << path_clusters_unwrapped_abelian_output << endl;
  cout << "path_clusters_wrapped_abelian_output "
       << path_clusters_wrapped_abelian_output << endl;
  cout << "path_windings_abelian_output " << path_windings_abelian_output
       << endl;
  cout << "path_clusters_unwrapped_monopole_output "
       << path_clusters_unwrapped_monopole_output << endl;
  cout << "path_clusters_wrapped_monopole_output "
       << path_clusters_wrapped_monopole_output << endl;
  cout << "path_clusters_windings_output " << path_windings_monopole_output
       << endl;
  cout << "path_clusters_unwrapped_monopoless_output "
       << path_clusters_unwrapped_monopoless_output << endl;
  cout << "path_clusters_wrapped_monopoless_output "
       << path_clusters_wrapped_monopoless_output << endl;
  cout << "path_windings_monopoless_output " << path_windings_monopoless_output
       << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "path_inverse_laplacian " << path_inverse_laplacian << endl;
  cout << "N_dir_gevp " << N_dir_gevp << endl;
  cout << "HYP_alpha1 " << HYP_alpha1 << endl;
  cout << "HYP_alpha2 " << HYP_alpha2 << endl;
  cout << "HYP_alpha3 " << HYP_alpha3 << endl;
  cout << "APE_alpha " << APE_alpha << endl;
  cout << "APE_steps " << APE_steps << endl;
  cout << "HYP_steps " << HYP_steps << endl;
  cout << "copies_required " << copies_required << endl;
  cout << "x_size " << x_size1 << endl;
  cout << "y_size " << y_size1 << endl;
  cout << "z_size " << z_size1 << endl;
  cout << "t_size " << t_size1 << endl;

  // x_size = x_size1;
  // y_size = y_size1;
  // z_size = z_size1;
  // t_size = t_size1;

  string path_functional;
  string path_wilson_loops_abelian;
  string path_wilson_loops_monopole;
  string path_clusters_unwrapped_abelian;
  string path_clusters_wrapped_abelian;
  string path_windings_abelian;
  string path_clusters_unwrapped_monopole;
  string path_clusters_wrapped_monopole;
  string path_windings_monopole;
  string path_clusters_unwrapped_monopoless;
  string path_clusters_wrapped_monopoless;
  string path_windings_monopoless;
  double MAG_functional;
  for (int copy = 0; copy < copies_required; copy++) {
    path_functional = path_functional_output + "_" + to_string(copy);
    path_wilson_loops_abelian =
        path_wilson_loops_abelian_output + "_" + to_string(copy);
    path_wilson_loops_monopole =
        path_wilson_loops_monopole_output + "_" + to_string(copy);
    path_clusters_unwrapped_abelian =
        path_clusters_unwrapped_abelian_output + "_" + to_string(copy);
    path_clusters_wrapped_abelian =
        path_clusters_wrapped_abelian_output + "_" + to_string(copy);
    path_windings_abelian =
        path_windings_abelian_output + "_" + to_string(copy);
    path_clusters_unwrapped_monopole =
        path_clusters_unwrapped_monopole_output + "_" + to_string(copy);
    path_clusters_wrapped_monopole =
        path_clusters_wrapped_monopole_output + "_" + to_string(copy);
    path_windings_monopole =
        path_windings_monopole_output + "_" + to_string(copy);
    path_clusters_unwrapped_monopoless =
        path_clusters_unwrapped_monopoless_output + "_" + to_string(copy);
    path_clusters_wrapped_monopoless =
        path_clusters_wrapped_monopoless_output + "_" + to_string(copy);
    path_windings_monopoless =
        path_windings_monopoless_output + "_" + to_string(copy);
    if (!(std::filesystem::exists(path_functional) &&
          std::filesystem::exists(path_wilson_loops_abelian) &&
          std::filesystem::exists(path_wilson_loops_monopole) &&
          std::filesystem::exists(path_clusters_unwrapped_abelian) &&
          std::filesystem::exists(path_clusters_wrapped_abelian) &&
          std::filesystem::exists(path_windings_abelian) &&
          std::filesystem::exists(path_clusters_unwrapped_monopole) &&
          std::filesystem::exists(path_clusters_wrapped_monopole) &&
          std::filesystem::exists(path_windings_monopole) &&
          std::filesystem::exists(path_clusters_unwrapped_monopoless) &&
          std::filesystem::exists(path_clusters_wrapped_monopoless) &&
          std::filesystem::exists(path_windings_monopoless))) {
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
      make_simulated_annealing(conf_su2, spins, 2.5, 0.1, 0.1, 6, 20);
      make_maximization_final(conf_su2, spins, 6, 1e-14, 1e-15);
      gauge_tranformation_spins(conf_su2, spins);
      std::cout << "MAG time: " << omp_get_wtime() - omp_time << std::endl;
      std::cout << "plaket after MAG: " << plaket(conf_su2) << std::endl;
      MAG_functional = MAG_functional_su2(conf_su2);
      std::cout << "MAG functional after MAG: " << MAG_functional << std::endl;
      std::vector<std::complex<double>> conf_complex =
          convert_to_complex(conf_su2);
      std::vector<std::complex<double>> gauge_complex =
          generate_gauge_complex_uniform(data_pattern);
      std::cout << "Landau functional before Landau gauge: "
                << Landau_functional_conf_complex(conf_complex, gauge_complex,
                                                  data_pattern)
                << std::endl;

      omp_time = omp_get_wtime();
      make_simulated_annealing(conf_complex, gauge_complex, data_pattern, 10,
                               0.1, 0.1, 4, 20);
      make_maximization_final(conf_complex, gauge_complex, data_pattern, 4,
                              1e-14, 1e-15);
      std::cout << "Landau gauge time: " << omp_get_wtime() - omp_time
                << std::endl;
      std::cout << "Landau functional after Landau gauge: "
                << Landau_functional_conf_complex(conf_complex, gauge_complex,
                                                  data_pattern)
                << std::endl;
      apply_gauge_Landau(gauge_complex, conf_su2);
      std::cout << "Landau functional after apply_gauge_Landau: "
                << Landau_functional(conf_su2.array) << std::endl;
      std::cout << "plaket after Landau gauge: " << plaket(conf_su2)
                << std::endl;
      std::cout << "MAG functional after Landau gauge: "
                << MAG_functional_su2(conf_su2) << std::endl;

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
      for (auto it = lengths_unwrapped_abelian.cbegin();
           it != lengths_unwrapped_abelian.cend(); ++it) {
        std::cout << it->first << "," << it->second << endl;
      }
      for (int i = 0; i < wrapped_lengths_abelian.size(); i++) {
        if (std::find(positions_percolating_abelian.begin(),
                      positions_percolating_abelian.end(),
                      i) != positions_percolating_abelian.end()) {
          std::cout << wrapped_lengths_abelian[i] << ","
                    << wrappings_abelian[i][0] << "," << wrappings_abelian[i][1]
                    << "," << wrappings_abelian[i][2] << ","
                    << wrappings_abelian[i][3] << ",percolating" << endl;
        } else {
          std::cout << wrapped_lengths_abelian[i] << ","
                    << wrappings_abelian[i][0] << "," << wrappings_abelian[i][1]
                    << "," << wrappings_abelian[i][2] << ","
                    << wrappings_abelian[i][3] << ",non-percolating" << endl;
        }
      }
      for (auto it = time_windings_abelian.begin();
           it != time_windings_abelian.end(); ++it) {
        std::cout << it->first << "," << it->second << ",time" << endl;
      }
      for (auto it = space_windings_abelian.begin();
           it != space_windings_abelian.end(); ++it) {
        std::cout << it->first << "," << it->second << ",space" << endl;
      }
      double asymmetry_abelian = (space_currents / 3. - time_currents) /
                                 (space_currents / 3. + time_currents);
      std::cout << asymmetry_abelian << endl;
      std::cout << "monopoles time: " << omp_get_wtime() - omp_time
                << std::endl;

      // abelian wilson loop
      omp_time = omp_get_wtime();
      Data::LatticeData<DataPatternLexicographical, abelian> conf_abelian1(
          {x_size1, y_size1, z_size1, t_size1});
      Data::LatticeData<DataPatternLexicographical, abelian> conf_abelian2(
          {x_size1, y_size1, z_size1, t_size1});
      map<tuple<int, int>, double> wilson_tmp;
      map<tuple<int, int, int, int>, double> wilson_loops_abelian;
      for (int dir = 0; dir < N_dir_gevp; dir++) {
        conf_abelian1 = conf_abelian;
        if (dir > 0) {
          conf_abelian1.swap_directions(dir - 1, 3);
        }
        conf_abelian2 = conf_abelian1;
        if (HYP_steps > 0) {
          for (int HYP_step = 1; HYP_step <= HYP_steps; HYP_step++) {
            smearing_HYP(conf_abelian1, HYP_alpha1, HYP_alpha2, HYP_alpha3);
          }
        }
        // wilson loops at (APE_step1, APE_step2) APE_steps, APE_step1 <
        // APE_step2
        for (int APE_step1 = 1; APE_step1 <= APE_steps; APE_step1++) {
          smearing_APE(conf_abelian1, APE_alpha);
          if ((APE_step1 - calculation_APE_start) % calculation_step_APE == 0 &&
              APE_step1 >= calculation_APE_start) {
            wilson_tmp =
                wilson_loop(conf_abelian1, 1, x_size1 / 2, 1, t_size1 / 2);
            write_wilson_loops(wilson_tmp, wilson_loops_abelian, APE_step1,
                               APE_step1);
            conf_abelian2 = conf_abelian1;
            for (int APE_step2 = APE_step1 + 1; APE_step2 <= APE_steps;
                 APE_step2++) {
              smearing_APE(conf_abelian2, APE_alpha);
              if ((APE_step2 - APE_step1) % calculation_step_APE == 0) {
                wilson_tmp =
                    wilson_gevp_indexed(conf_abelian1, conf_abelian2, 1,
                                        x_size1 / 2, 1, t_size1 / 2);
                write_wilson_loops(wilson_tmp, wilson_loops_abelian, APE_step1,
                                   APE_step2);
              }
            }
          }
        }
      }
      for (auto it = wilson_loops_abelian.begin();
           it != wilson_loops_abelian.end(); it++) {
        it->second = it->second / N_dir_gevp;
      }
      conf_abelian1.array.clear();
      conf_abelian1.array.shrink_to_fit();
      conf_abelian2.array.clear();
      conf_abelian2.array.shrink_to_fit();
      std::cout << "abelian wilson loops time: " << omp_get_wtime() - omp_time
                << std::endl;
      std::cout << std::endl << "wilson_loops: " << std::endl;
      for (auto it = wilson_loops_abelian.begin();
           it != wilson_loops_abelian.end(); it++) {
        std::cout << get<0>(it->first) << "," << get<1>(it->first) << ","
                  << get<2>(it->first) << "," << get<3>(it->first) << ","
                  << it->second << endl;
      }

      // decomposition
      omp_time = omp_get_wtime();
      vector<double> inverse_laplacian = read_inverse_laplacian(
          path_inverse_laplacian, data_pattern_laplacian);
      vector<vector<int>> dirac_plakets =
          calculate_monopole_plaket_singular(conf_abelian);
      vector<double> monopole_angles =
          make_monopole_angles(dirac_plakets, inverse_laplacian, data_pattern,
                               data_pattern_laplacian);
      Data::LatticeData<DataPatternLexicographical, abelian> conf_monopole(
          {x_size1, y_size1, z_size1, t_size1});
      Data::LatticeData<DataPatternLexicographical, abelian> conf_photon(
          {x_size1, y_size1, z_size1, t_size1});
      for (int i = 0; i < monopole_angles.size(); i++) {
        conf_monopole.array[i] = abelian(1, monopole_angles[i]);
        conf_photon.array[i] =
            abelian(1, conf_abelian[i].phi + monopole_angles[i]);
      }
      monopole_angles.clear();
      monopole_angles.shrink_to_fit();
      std::cout << "decomposition time: " << omp_get_wtime() - omp_time
                << std::endl;

      // monoples
      omp_time = omp_get_wtime();
      J = calculate_current(conf_monopole);
      LL = calculate_clusters(J, data_pattern);
      vector<vector<int>> wrappings_monopole;
      vector<int> wrapped_lengths_monopole;
      map<int, int> lengths_unwrapped_monopole;
      lengths_wrapped_time = map<int, int>();
      lengths_wrapped_space = map<int, int>();
      lengths_wrapped_both = map<int, int>();
      map<int, int> space_windings_monopole;
      map<int, int> time_windings_monopole;
      lengths_mu = vector<int>();
      currents = std::tuple<int, int>();
      space_currents = 0;
      time_currents = 0;
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
          lengths_unwrapped_monopole[length]++;
        } else {
          wrappings_monopole.push_back(
              vector<int>({lengths_mu[0] / conf_monopole.lat_dim[0],
                           lengths_mu[1] / conf_monopole.lat_dim[1],
                           lengths_mu[2] / conf_monopole.lat_dim[2],
                           lengths_mu[3] / conf_monopole.lat_dim[3]}));
          wrapped_lengths_monopole.push_back(length);
        }
      }
      // sort wrapped_lengths and wrappings accordingly
      permutations = std::vector<std::size_t>(wrapped_lengths_monopole.size());
      std::iota(permutations.begin(), permutations.end(), 0);
      std::sort(permutations.begin(), permutations.end(),
                [&](std::size_t i, std::size_t j) {
                  return wrapped_lengths_monopole[i] >
                         wrapped_lengths_monopole[j];
                });
      sorted_lengths = std::vector<int>(permutations.size());
      std::transform(
          permutations.begin(), permutations.end(), sorted_lengths.begin(),
          [&](std::size_t i) { return wrapped_lengths_monopole[i]; });
      wrapped_lengths_monopole = sorted_lengths;
      sorted_wrappings = std::vector<std::vector<int>>(permutations.size(),
                                                       std::vector<int>());
      std::transform(permutations.begin(), permutations.end(),
                     sorted_wrappings.begin(),
                     [&](std::size_t i) { return wrappings_monopole[i]; });
      wrappings_monopole = sorted_wrappings;
      vector<int> positions_percolating_monopole;
      if (wrappings_monopole.size() > 1) {
        positions_percolating_monopole = group_percolating(wrappings_monopole);
      }
      for (auto it = lengths_unwrapped_monopole.cbegin();
           it != lengths_unwrapped_monopole.cend(); ++it) {
        std::cout << it->first << "," << it->second << endl;
      }
      for (int i = 0; i < wrapped_lengths_monopole.size(); i++) {
        if (std::find(positions_percolating_monopole.begin(),
                      positions_percolating_monopole.end(),
                      i) != positions_percolating_monopole.end()) {
          std::cout << wrapped_lengths_monopole[i] << ","
                    << wrappings_monopole[i][0] << ","
                    << wrappings_monopole[i][1] << ","
                    << wrappings_monopole[i][2] << ","
                    << wrappings_monopole[i][3] << ",percolating" << endl;
        } else {
          std::cout << wrapped_lengths_monopole[i] << ","
                    << wrappings_monopole[i][0] << ","
                    << wrappings_monopole[i][1] << ","
                    << wrappings_monopole[i][2] << ","
                    << wrappings_monopole[i][3] << ",non-percolating" << endl;
        }
      }
      for (auto it = time_windings_monopole.begin();
           it != time_windings_monopole.end(); ++it) {
        std::cout << it->first << "," << it->second << ",time" << endl;
      }
      for (auto it = space_windings_monopole.begin();
           it != space_windings_monopole.end(); ++it) {
        std::cout << it->first << "," << it->second << ",space" << endl;
      }
      double asymmetry_monopole = (space_currents / 3. - time_currents) /
                                  (space_currents / 3. + time_currents);
      std::cout << asymmetry_monopole << endl;
      std::cout << "monopoles from monopole angles time: "
                << omp_get_wtime() - omp_time << std::endl;

      omp_time = omp_get_wtime();
      J = calculate_current(conf_photon);
      LL = calculate_clusters(J, data_pattern);
      vector<vector<int>> wrappings_monopoless;
      vector<int> wrapped_lengths_monopoless;
      map<int, int> lengths_unwrapped_monopoless;
      lengths_wrapped_time = map<int, int>();
      lengths_wrapped_space = map<int, int>();
      lengths_wrapped_both = map<int, int>();
      map<int, int> space_windings_monopoless;
      map<int, int> time_windings_monopoless;
      lengths_mu = vector<int>();
      currents = std::tuple<int, int>();
      space_currents = 0;
      time_currents = 0;
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
          lengths_unwrapped_monopoless[length]++;
        } else {
          wrappings_monopoless.push_back(
              vector<int>({lengths_mu[0] / conf_photon.lat_dim[0],
                           lengths_mu[1] / conf_photon.lat_dim[1],
                           lengths_mu[2] / conf_photon.lat_dim[2],
                           lengths_mu[3] / conf_photon.lat_dim[3]}));
          wrapped_lengths_monopoless.push_back(length);
        }
      }
      // sort wrapped_lengths and wrappings accordingly
      permutations =
          std::vector<std::size_t>(wrapped_lengths_monopoless.size());
      std::iota(permutations.begin(), permutations.end(), 0);
      std::sort(permutations.begin(), permutations.end(),
                [&](std::size_t i, std::size_t j) {
                  return wrapped_lengths_monopoless[i] >
                         wrapped_lengths_monopoless[j];
                });
      sorted_lengths = std::vector<int>(permutations.size());
      std::transform(
          permutations.begin(), permutations.end(), sorted_lengths.begin(),
          [&](std::size_t i) { return wrapped_lengths_monopoless[i]; });
      wrapped_lengths_monopoless = sorted_lengths;
      sorted_wrappings = std::vector<std::vector<int>>(permutations.size(),
                                                       std::vector<int>());
      std::transform(permutations.begin(), permutations.end(),
                     sorted_wrappings.begin(),
                     [&](std::size_t i) { return wrappings_monopoless[i]; });
      wrappings_monopoless = sorted_wrappings;
      vector<int> positions_percolating_monopoless;
      if (wrappings_monopoless.size() > 1) {
        positions_percolating_monopoless =
            group_percolating(wrappings_monopoless);
      }
      for (auto it = lengths_unwrapped_monopoless.cbegin();
           it != lengths_unwrapped_monopoless.cend(); ++it) {
        std::cout << it->first << "," << it->second << endl;
      }
      for (int i = 0; i < wrapped_lengths_monopoless.size(); i++) {
        if (std::find(positions_percolating_monopoless.begin(),
                      positions_percolating_monopoless.end(),
                      i) != positions_percolating_monopoless.end()) {
          std::cout << wrapped_lengths_monopoless[i] << ","
                    << wrappings_monopoless[i][0] << ","
                    << wrappings_monopoless[i][1] << ","
                    << wrappings_monopoless[i][2] << ","
                    << wrappings_monopoless[i][3] << ",percolating" << endl;
        } else {
          std::cout << wrapped_lengths_monopoless[i] << ","
                    << wrappings_monopoless[i][0] << ","
                    << wrappings_monopoless[i][1] << ","
                    << wrappings_monopoless[i][2] << ","
                    << wrappings_monopoless[i][3] << ",non-percolating" << endl;
        }
      }
      for (auto it = time_windings_monopoless.begin();
           it != time_windings_monopoless.end(); ++it) {
        std::cout << it->first << "," << it->second << ",time" << endl;
      }
      for (auto it = space_windings_monopoless.begin();
           it != space_windings_monopoless.end(); ++it) {
        std::cout << it->first << "," << it->second << ",space" << endl;
      }
      double asymmetry_monopoless = (space_currents / 3. - time_currents) /
                                    (space_currents / 3. + time_currents);
      std::cout << asymmetry_monopoless << endl;
      std::cout << "monopoles from photon angles time: "
                << omp_get_wtime() - omp_time << std::endl;

      // monopole wilson loop
      omp_time = omp_get_wtime();
      Data::LatticeData<DataPatternLexicographical, abelian> conf_monopole1(
          {x_size1, y_size1, z_size1, t_size1});
      Data::LatticeData<DataPatternLexicographical, abelian> conf_monopole2(
          {x_size1, y_size1, z_size1, t_size1});
      wilson_tmp = map<tuple<int, int>, double>();
      map<tuple<int, int, int, int>, double> wilson_loops_monopole;
      for (int dir = 0; dir < N_dir_gevp; dir++) {
        conf_monopole1 = conf_monopole;
        if (dir > 0) {
          conf_monopole1.swap_directions(dir - 1, 3);
        }
        conf_monopole2 = conf_monopole1;
        if (HYP_steps > 0) {
          for (int HYP_step = 1; HYP_step <= HYP_steps; HYP_step++) {
            smearing_HYP(conf_monopole1, HYP_alpha1, HYP_alpha2, HYP_alpha3);
          }
        }
        // wilson loops at (APE_step1, APE_step2) APE_steps, APE_step1 <
        // APE_step2
        for (int APE_step1 = 1; APE_step1 <= APE_steps; APE_step1++) {
          smearing_APE(conf_monopole1, APE_alpha);
          if ((APE_step1 - calculation_APE_start) % calculation_step_APE == 0 &&
              APE_step1 >= calculation_APE_start) {
            wilson_tmp =
                wilson_loop(conf_monopole1, 1, x_size1 / 2, 1, t_size1 / 2);
            write_wilson_loops(wilson_tmp, wilson_loops_monopole, APE_step1,
                               APE_step1);
            conf_monopole2 = conf_monopole1;
            for (int APE_step2 = APE_step1 + 1; APE_step2 <= APE_steps;
                 APE_step2++) {
              smearing_APE(conf_monopole2, APE_alpha);
              if ((APE_step2 - APE_step1) % calculation_step_APE == 0) {
                wilson_tmp =
                    wilson_gevp_indexed(conf_monopole1, conf_monopole2, 1,
                                        x_size1 / 2, 1, t_size1 / 2);
                write_wilson_loops(wilson_tmp, wilson_loops_monopole, APE_step1,
                                   APE_step2);
              }
            }
          }
        }
      }
      for (auto it = wilson_loops_monopole.begin();
           it != wilson_loops_monopole.end(); it++) {
        it->second = it->second / N_dir_gevp;
      }
      conf_monopole1.array.clear();
      conf_monopole1.array.shrink_to_fit();
      conf_monopole2.array.clear();
      conf_monopole2.array.shrink_to_fit();
      std::cout << "monopole wilson loops time: " << omp_get_wtime() - omp_time
                << std::endl;
      std::cout << std::endl << "wilson_loops: " << std::endl;
      for (auto it = wilson_loops_monopole.begin();
           it != wilson_loops_monopole.end(); it++) {
        std::cout << get<0>(it->first) << "," << get<1>(it->first) << ","
                  << get<2>(it->first) << "," << get<3>(it->first) << ","
                  << it->second << endl;
      }

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

      std::ofstream output_stream_clusters_unwrapped_monopole(
          path_clusters_unwrapped_monopole);
      output_stream_clusters_unwrapped_monopole << "length,number" << endl;
      for (auto it = lengths_unwrapped_monopole.cbegin();
           it != lengths_unwrapped_monopole.cend(); ++it) {
        output_stream_clusters_unwrapped_monopole << it->first << ","
                                                  << it->second << endl;
      }
      output_stream_clusters_unwrapped_monopole.close();
      std::ofstream output_stream_clusters_wrapped_monopole(
          path_clusters_wrapped_monopole);
      output_stream_clusters_wrapped_monopole
          << "length,x0_wrap,x1_wrap,x2_wrap,x3_wrap,percolating_group" << endl;
      for (int i = 0; i < wrapped_lengths_monopole.size(); i++) {
        if (std::find(positions_percolating_monopole.begin(),
                      positions_percolating_monopole.end(),
                      i) != positions_percolating_monopole.end()) {
          output_stream_clusters_wrapped_monopole
              << wrapped_lengths_monopole[i] << "," << wrappings_monopole[i][0]
              << "," << wrappings_monopole[i][1] << ","
              << wrappings_monopole[i][2] << "," << wrappings_monopole[i][3]
              << ",percolating" << endl;
        } else {
          output_stream_clusters_wrapped_monopole
              << wrapped_lengths_monopole[i] << "," << wrappings_monopole[i][0]
              << "," << wrappings_monopole[i][1] << ","
              << wrappings_monopole[i][2] << "," << wrappings_monopole[i][3]
              << ",non-percolating" << endl;
        }
      }
      output_stream_clusters_wrapped_monopole.close();
      std::ofstream output_stream_windings_monopole(path_windings_monopole);
      output_stream_windings_monopole
          << "winding_number,cluster_number,direction" << endl;
      for (auto it = time_windings_monopole.begin();
           it != time_windings_monopole.end(); ++it) {
        output_stream_windings_monopole << it->first << "," << it->second
                                        << ",time" << endl;
      }
      for (auto it = space_windings_monopole.begin();
           it != space_windings_monopole.end(); ++it) {
        output_stream_windings_monopole << it->first << "," << it->second
                                        << ",space" << endl;
      }
      output_stream_windings_monopole.close();

      std::ofstream output_stream_clusters_unwrapped_monopoless(
          path_clusters_unwrapped_monopoless);
      output_stream_clusters_unwrapped_monopoless << "length,number" << endl;
      for (auto it = lengths_unwrapped_monopoless.cbegin();
           it != lengths_unwrapped_monopoless.cend(); ++it) {
        output_stream_clusters_unwrapped_monopoless << it->first << ","
                                                    << it->second << endl;
      }
      output_stream_clusters_unwrapped_monopoless.close();
      std::ofstream output_stream_clusters_wrapped_monopoless(
          path_clusters_wrapped_monopoless);
      output_stream_clusters_wrapped_monopoless
          << "length,x0_wrap,x1_wrap,x2_wrap,x3_wrap,percolating_group" << endl;
      for (int i = 0; i < wrapped_lengths_monopoless.size(); i++) {
        if (std::find(positions_percolating_monopoless.begin(),
                      positions_percolating_monopoless.end(),
                      i) != positions_percolating_monopoless.end()) {
          output_stream_clusters_wrapped_monopoless
              << wrapped_lengths_monopoless[i] << ","
              << wrappings_monopoless[i][0] << "," << wrappings_monopoless[i][1]
              << "," << wrappings_monopoless[i][2] << ","
              << wrappings_monopoless[i][3] << ",percolating" << endl;
        } else {
          output_stream_clusters_wrapped_monopoless
              << wrapped_lengths_monopoless[i] << ","
              << wrappings_monopoless[i][0] << "," << wrappings_monopoless[i][1]
              << "," << wrappings_monopoless[i][2] << ","
              << wrappings_monopoless[i][3] << ",non-percolating" << endl;
        }
      }
      output_stream_clusters_wrapped_monopoless.close();
      std::ofstream output_stream_windings_monopoless(path_windings_monopoless);
      output_stream_windings_monopoless
          << "winding_number,cluster_number,direction" << endl;
      for (auto it = time_windings_monopoless.begin();
           it != time_windings_monopoless.end(); ++it) {
        output_stream_windings_monopoless << it->first << "," << it->second
                                          << ",time" << endl;
      }
      for (auto it = space_windings_monopoless.begin();
           it != space_windings_monopoless.end(); ++it) {
        output_stream_windings_monopoless << it->first << "," << it->second
                                          << ",space" << endl;
      }
      output_stream_windings_monopoless.close();

      std::ofstream output_stream_wilson_loops_abelian(
          path_wilson_loops_abelian);
      output_stream_wilson_loops_abelian
          << "smearing_step1,smearing_step2,time_size,space_size,wilson_loop"
          << endl;
      for (auto it = wilson_loops_abelian.begin();
           it != wilson_loops_abelian.end(); it++) {
        output_stream_wilson_loops_abelian
            << get<0>(it->first) << "," << get<1>(it->first) << ","
            << get<2>(it->first) << "," << get<3>(it->first) << ","
            << it->second << endl;
      }
      output_stream_wilson_loops_abelian.close();

      std::ofstream output_stream_wilson_loops_monopole(
          path_wilson_loops_monopole);
      output_stream_wilson_loops_monopole
          << "smearing_step1,smearing_step2,time_size,space_size,wilson_loop"
          << endl;
      for (auto it = wilson_loops_monopole.begin();
           it != wilson_loops_monopole.end(); it++) {
        output_stream_wilson_loops_monopole
            << get<0>(it->first) << "," << get<1>(it->first) << ","
            << get<2>(it->first) << "," << get<3>(it->first) << ","
            << it->second << endl;
      }
      output_stream_wilson_loops_monopole.close();
    }
  }
}