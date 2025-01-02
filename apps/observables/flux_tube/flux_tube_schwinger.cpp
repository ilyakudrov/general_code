#include "../../../lib/cpu/include/basic_observables.h"
#include "../../../lib/cpu/include/data.h"
#include "../../../lib/cpu/include/flux_tube.h"
#include "../../../lib/cpu/include/matrix.h"

#include <algorithm>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <vector>

#ifndef MATRIX_PLAKET
#define MATRIX_PLAKET su2
#endif
#ifndef MATRIX_WILSON
#define MATRIX_WILSON su2
#endif

using namespace std;

int x_size;
int y_size;
int z_size;
int t_size;

int main(int argc, char *argv[]) {
  double start_time;
  double end_time;
  double search_time;

  std::string conf_format_plaket;
  std::string conf_format_wilson;
  std::string conf_path_plaket;
  std::string conf_path_wilson;
  std::string output_path_electric_long_l;
  std::string output_path_electric_long_tr;
  std::string output_path_magnetic_long_l;
  std::string output_path_magnetic_long_tr;
  std::string output_path_electric_trans_l;
  std::string output_path_electric_trans_tr;
  std::string output_path_magnetic_trans_l;
  std::string output_path_magnetic_trans_tr;
  int L_spat, L_time;
  int x_trans = 0;
  int bytes_skip_plaket = 0;
  int bytes_skip_wilson = 0;
  int T_min, T_max;
  int R_min, R_max;
  int d_ouside = 10;
  int d_max = 10;
  bool convert_plaket = 0;
  bool convert_wilson = 0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-conf_format_plaket") == 0) {
      conf_format_plaket = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip_plaket") {
      bytes_skip_plaket = stoi(string(argv[++i]));
    } else if (std::string(argv[i]) == "-conf_format_wilson") {
      conf_format_wilson = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip_wilson") {
      bytes_skip_wilson = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-convert_wilson") {
      istringstream(string(argv[++i])) >> convert_wilson;
    } else if (string(argv[i]) == "-convert_plaket") {
      istringstream(string(argv[++i])) >> convert_plaket;
    } else if (std::string(argv[i]) == "-conf_path_plaket") {
      conf_path_plaket = argv[++i];
    } else if (std::string(argv[i]) == "-conf_path_wilson") {
      conf_path_wilson = argv[++i];
    } else if (std::string(argv[i]) == "-output_path_electric_long_l") {
      output_path_electric_long_l = argv[++i];
    } else if (std::string(argv[i]) == "-output_path_electric_long_tr") {
      output_path_electric_long_tr = argv[++i];
    } else if (std::string(argv[i]) == "-output_path_magnetic_long_l") {
      output_path_magnetic_long_l = argv[++i];
    } else if (std::string(argv[i]) == "-output_path_magnetic_long_tr") {
      output_path_magnetic_long_tr = argv[++i];
    } else if (std::string(argv[i]) == "-output_path_electric_trans_l") {
      output_path_electric_trans_l = argv[++i];
    } else if (std::string(argv[i]) == "-output_path_electric_trans_tr") {
      output_path_electric_trans_tr = argv[++i];
    } else if (std::string(argv[i]) == "-output_path_magnetic_trans_l") {
      output_path_magnetic_trans_l = argv[++i];
    } else if (std::string(argv[i]) == "-output_path_magnetic_trans_tr") {
      output_path_magnetic_trans_tr = argv[++i];
    } else if (std::string(argv[i]) == "-L_spat") {
      L_spat = stoi(std::string(argv[++i]));
    } else if (std::string(argv[i]) == "-L_time") {
      L_time = stoi(std::string(argv[++i]));
    } else if (std::string(argv[i]) == "-x_trans") {
      x_trans = stoi(std::string(argv[++i]));
    } else if (std::string(argv[i]) == "-T_min") {
      T_min = stoi(std::string(argv[++i]));
    } else if (std::string(argv[i]) == "-T_max") {
      T_max = stoi(std::string(argv[++i]));
    } else if (std::string(argv[i]) == "-R_min") {
      R_min = stoi(std::string(argv[++i]));
    } else if (std::string(argv[i]) == "-R_max") {
      R_max = stoi(std::string(argv[++i]));
    } else if (std::string(argv[i]) == "-d_ouside") {
      d_ouside = stoi(std::string(argv[++i]));
    } else if (std::string(argv[i]) == "-d_max") {
      d_max = stoi(std::string(argv[++i]));
    }
  }

  std::cout << "conf_format_plaket " << conf_format_plaket << std::endl;
  cout << "bytes_skip_plaket " << bytes_skip_plaket << endl;
  cout << "convert_plaket " << convert_plaket << endl;
  std::cout << "conf_format_wilson " << conf_format_wilson << std::endl;
  cout << "bytes_skip_wilson " << bytes_skip_wilson << endl;
  cout << "convert_wilson " << convert_wilson << endl;
  std::cout << "conf_path_plaket " << conf_path_plaket << std::endl;
  std::cout << "conf_path_wilson " << conf_path_wilson << std::endl;
  std::cout << "output_path_electric_long_l " << output_path_electric_long_l
            << std::endl;
  std::cout << "output_path_electric_long_tr " << output_path_electric_long_tr
            << std::endl;
  std::cout << "output_path_magnetic_long_l " << output_path_magnetic_long_l
            << std::endl;
  std::cout << "output_path_magnetic_long_tr " << output_path_magnetic_long_tr
            << std::endl;
  std::cout << "output_path_electric_trans_l " << output_path_electric_trans_l
            << std::endl;
  std::cout << "output_path_electric_trans_tr " << output_path_electric_trans_tr
            << std::endl;
  std::cout << "output_path_magnetic_trans_l " << output_path_magnetic_trans_l
            << std::endl;
  std::cout << "output_path_magnetic_trans_tr " << output_path_magnetic_trans_tr
            << std::endl;
  std::cout << "L_spat " << L_spat << std::endl;
  std::cout << "L_time " << L_time << std::endl;
  std::cout << "R_min " << R_min << std::endl;
  std::cout << "R_max " << R_max << std::endl;
  std::cout << "T_min " << T_min << std::endl;
  std::cout << "T_max " << T_max << std::endl;
  std::cout << "d_ouside " << d_ouside << std::endl;
  std::cout << "d_max " << d_ouside << std::endl;
  std::cout << "x_trans " << x_trans << std::endl;

  x_size = L_spat;
  y_size = L_spat;
  z_size = L_spat;
  t_size = L_time;

  Data::data<MATRIX_PLAKET> conf_plaket;
  Data::data<MATRIX_WILSON> conf_wilson;

  get_data(conf_plaket, conf_path_plaket, conf_format_plaket, bytes_skip_plaket,
           convert_plaket);
  // vector<float> conf_full = read_full_ml5(conf_path_plaket, 1);
  // conf_plaket.read_float_ml5(conf_full, 0);
  get_data(conf_wilson, conf_path_wilson, conf_format_wilson, bytes_skip_wilson,
           convert_wilson);
  // conf_full = read_full_ml5(conf_path_wilson, 1);
  // conf_wilson.read_float_ml5(conf_full, 0);

  double plaket_time_average = plaket_time(conf_plaket.array);
  double plaket_space_average = plaket_space(conf_plaket.array);

  std::cout << "plaket_time " << plaket_time_average << " smeared_plaket_time "
            << plaket_time(conf_wilson.array) << std::endl;
  std::cout << "plaket_space " << plaket_space_average
            << " smeared_plaket_space " << plaket_space(conf_wilson.array)
            << std::endl;

  std::ofstream stream_electric_long_l;
  std::ofstream stream_electric_long_tr;
  std::ofstream stream_magnetic_long_l;
  std::ofstream stream_magnetic_long_tr;
  std::ofstream stream_electric_trans_l;
  std::ofstream stream_electric_trans_tr;
  std::ofstream stream_magnetic_trans_l;
  std::ofstream stream_magnetic_trans_tr;

  stream_electric_long_l.precision(17);
  stream_electric_long_tr.precision(17);
  stream_magnetic_long_l.precision(17);
  stream_magnetic_long_tr.precision(17);
  stream_electric_trans_l.precision(17);
  stream_electric_trans_tr.precision(17);
  stream_magnetic_trans_l.precision(17);
  stream_magnetic_trans_tr.precision(17);

  stream_electric_long_l.open(output_path_electric_long_l);
  stream_electric_long_tr.open(output_path_electric_long_tr);
  stream_magnetic_long_l.open(output_path_magnetic_long_l);
  stream_magnetic_long_tr.open(output_path_magnetic_long_tr);
  stream_electric_trans_l.open(output_path_electric_trans_l);
  stream_electric_trans_tr.open(output_path_electric_trans_tr);
  stream_magnetic_trans_l.open(output_path_magnetic_trans_l);
  stream_magnetic_trans_tr.open(output_path_magnetic_trans_tr);

  stream_electric_long_l
      << "T,R,d,correlator_schwinger,correlator_wilson,wilson_loop" << endl;
  stream_electric_long_tr
      << "T,R,d,correlator_schwinger,correlator_wilson,wilson_loop" << endl;
  stream_magnetic_long_l
      << "T,R,d,correlator_schwinger,correlator_wilson,wilson_loop" << endl;
  stream_magnetic_long_tr
      << "T,R,d,correlator_schwinger,correlator_wilson,wilson_loop" << endl;
  stream_electric_trans_l
      << "T,R,d,correlator_schwinger,correlator_wilson,wilson_loop" << endl;
  stream_electric_trans_tr
      << "T,R,d,correlator_schwinger,correlator_wilson,wilson_loop" << endl;
  stream_magnetic_trans_l
      << "T,R,d,correlator_schwinger,correlator_wilson,wilson_loop" << endl;
  stream_magnetic_trans_tr
      << "T,R,d,correlator_schwinger,correlator_wilson,wilson_loop" << endl;

  int schwinger_line_max = max(d_max, d_ouside);
  schwinger_line_max = max(schwinger_line_max, (R_max + 1) / 2);

  std::vector<std::vector<MATRIX_WILSON>> schwinger_lines_short(
      schwinger_line_max, std::vector<MATRIX_WILSON>());

  for (int d = 0; d < schwinger_line_max; d++) {
    schwinger_lines_short[d] =
        calculate_schwinger_lines_short(conf_wilson.array, d + 1);
  }

  std::map<std::tuple<int, int, int>, double>
      flux_tube_schwinger_electric_long_l;
  std::map<std::tuple<int, int, int>, double>
      flux_tube_schwinger_electric_long_tr;
  std::map<std::tuple<int, int, int>, double>
      flux_tube_schwinger_electric_trans_l;
  std::map<std::tuple<int, int, int>, double>
      flux_tube_schwinger_electric_trans_tr;
  std::map<std::tuple<int, int, int>, double>
      flux_tube_schwinger_magnetic_long_l;
  std::map<std::tuple<int, int, int>, double>
      flux_tube_schwinger_magnetic_long_tr;
  std::map<std::tuple<int, int, int>, double>
      flux_tube_schwinger_magnetic_trans_l;
  std::map<std::tuple<int, int, int>, double>
      flux_tube_schwinger_magnetic_trans_tr;

  // start_time = omp_get_wtime();

  // std::map<std::tuple<int, int, int>, double>
  //     flux_tube_schwinger_electric_long_l =
  //         flux_schwinger_electric_longitudinal_l(conf_wilson.array,
  //                                                schwinger_lines_short,
  //                                                T_min, T_max, R_min, R_max,
  //                                                d_ouside);

  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // cout << "flux_schwinger_electric_long_l: " << search_time << endl;

  // start_time = omp_get_wtime();

  // std::map<std::tuple<int, int, int>, double>
  //     flux_tube_schwinger_electric_long_tr =
  //         flux_schwinger_electric_longitudinal_tr(
  //             conf_wilson.array, schwinger_lines_short, T_min, T_max, R_min,
  //             R_max, d_ouside);

  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // cout << "flux_schwinger_electric_long_tr: " << search_time << endl;

  // start_time = omp_get_wtime();

  // std::map<std::tuple<int, int, int>, double>
  //     flux_tube_schwinger_electric_trans_l =
  //         flux_schwinger_electric_transversal_l(conf_wilson.array,
  //                                               schwinger_lines_short, T_min,
  //                                               T_max, R_min, R_max, d_max);

  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // cout << "flux_schwinger_electric_trans_l: " << search_time << endl;

  // start_time = omp_get_wtime();

  // std::map<std::tuple<int, int, int>, double>
  //     flux_tube_schwinger_electric_trans_tr =
  //         flux_schwinger_electric_transversal_tr(conf_wilson.array,
  //                                                schwinger_lines_short,
  //                                                T_min, T_max, R_min, R_max,
  //                                                d_max);

  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // cout << "flux_schwinger_electric_trans_tr: " << search_time << endl;

  // start_time = omp_get_wtime();

  // std::map<std::tuple<int, int, int>, double>
  //     flux_tube_schwinger_magnetic_long_l =
  //         flux_schwinger_magnetic_longitudinal_l(conf_wilson.array,
  //                                                schwinger_lines_short,
  //                                                T_min, T_max, R_min, R_max,
  //                                                d_ouside);

  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // cout << "flux_schwinger_magnetic_long_l: " << search_time << endl;

  // start_time = omp_get_wtime();

  // std::map<std::tuple<int, int, int>, double>
  //     flux_tube_schwinger_magnetic_long_tr =
  //         flux_schwinger_magnetic_longitudinal_tr(
  //             conf_wilson.array, schwinger_lines_short, T_min, T_max, R_min,
  //             R_max, d_ouside);

  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // cout << "flux_schwinger_magnetic_long_tr: " << search_time << endl;

  // start_time = omp_get_wtime();

  // std::map<std::tuple<int, int, int>, double>
  //     flux_tube_schwinger_magnetic_trans_l =
  //         flux_schwinger_magnetic_transversal_l(conf_wilson.array,
  //                                               schwinger_lines_short, T_min,
  //                                               T_max, R_min, R_max, d_max);

  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // cout << "flux_schwinger_magnetic_trans_l: " << search_time << endl;

  // start_time = omp_get_wtime();

  // std::map<std::tuple<int, int, int>, double>
  //     flux_tube_schwinger_magnetic_trans_tr =
  //         flux_schwinger_magnetic_transversal_tr(conf_wilson.array,
  //                                                schwinger_lines_short,
  //                                                T_min, T_max, R_min, R_max,
  //                                                d_max);

  // end_time = omp_get_wtime();
  // search_time = end_time - start_time;
  // cout << "flux_schwinger_magnetic_trans_tr: " << search_time << endl;

  start_time = omp_get_wtime();

  flux_schwinger_all(
      conf_wilson.array, schwinger_lines_short,
      flux_tube_schwinger_electric_long_l, flux_tube_schwinger_electric_long_tr,
      flux_tube_schwinger_electric_trans_l,
      flux_tube_schwinger_electric_trans_tr,
      flux_tube_schwinger_magnetic_long_l, flux_tube_schwinger_magnetic_long_tr,
      flux_tube_schwinger_magnetic_trans_l,
      flux_tube_schwinger_magnetic_trans_tr, T_min, T_max, R_min, R_max,
      d_ouside, d_max);

  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux_schwinger all: " << search_time << endl;
  vector<vector<MATRIX_WILSON>> separated_wilson =
      separate_wilson(conf_wilson.array);

  map<tuple<int, int>, double> wilson_loops =
      wilson_parallel(separated_wilson, R_min, R_max, T_min, T_max);

  map<tuple<int, int, int>, double> flux_tube_wilson_electric_long_l;
  map<tuple<int, int, int>, double> flux_tube_wilson_electric_long_tr;
  map<tuple<int, int, int>, double> flux_tube_wilson_electric_trans_l;
  map<tuple<int, int, int>, double> flux_tube_wilson_electric_trans_tr;
  map<tuple<int, int, int>, double> flux_tube_wilson_magnetic_long_l;
  map<tuple<int, int, int>, double> flux_tube_wilson_magnetic_long_tr;
  map<tuple<int, int, int>, double> flux_tube_wilson_magnetic_trans_l;
  map<tuple<int, int, int>, double> flux_tube_wilson_magnetic_trans_tr;

  vector<double> plaket_tr = calculate_plaket_time_trace_l(conf_plaket.array);

  start_time = omp_get_wtime();
  flux_tube_wilson_electric_long_l =
      calculate_wilson_plaket_correlator_electric_longitudinal(
          plaket_tr, conf_wilson.array, T_min, T_max, R_min, R_max, d_ouside);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube wilson longitudinal_l electric time: " << search_time
       << endl;

  start_time = omp_get_wtime();
  flux_tube_wilson_electric_trans_l =
      calculate_wilson_plaket_correlator_electric_transversal(
          plaket_tr, conf_wilson.array, T_min, T_max, R_min, R_max, d_max);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube wilson transversal_l electric time: " << search_time
       << endl;

  plaket_tr = calculate_plaket_time_trace_tr(conf_plaket.array);

  start_time = omp_get_wtime();
  flux_tube_wilson_electric_long_tr =
      calculate_wilson_plaket_correlator_electric_longitudinal(
          plaket_tr, conf_wilson.array, T_min, T_max, R_min, R_max, d_ouside);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube wilson longitudinal_tr electric time: " << search_time
       << endl;

  start_time = omp_get_wtime();
  flux_tube_wilson_electric_trans_tr =
      calculate_wilson_plaket_correlator_electric_transversal(
          plaket_tr, conf_wilson.array, T_min, T_max, R_min, R_max, d_ouside);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube wilson transversal_tr electric time: " << search_time
       << endl;

  plaket_tr = calculate_plaket_space_trace_l(conf_plaket.array);

  start_time = omp_get_wtime();
  flux_tube_wilson_magnetic_long_l =
      calculate_wilson_plaket_correlator_electric_longitudinal(
          plaket_tr, conf_wilson.array, T_min, T_max, R_min, R_max, d_ouside);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube wilson longitudinal_l magnetic time: " << search_time
       << endl;

  start_time = omp_get_wtime();
  flux_tube_wilson_magnetic_trans_l =
      calculate_wilson_plaket_correlator_electric_transversal(
          plaket_tr, conf_wilson.array, T_min, T_max, R_min, R_max, d_max);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube wilson transversal_l magnetic time: " << search_time
       << endl;

  plaket_tr = calculate_plaket_space_trace_tr(conf_plaket.array);

  start_time = omp_get_wtime();
  flux_tube_wilson_magnetic_long_tr =
      calculate_wilson_plaket_correlator_electric_longitudinal(
          plaket_tr, conf_wilson.array, T_min, T_max, R_min, R_max, d_ouside);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube wilson longitudinal_tr magnetic time: " << search_time
       << endl;

  start_time = omp_get_wtime();
  flux_tube_wilson_magnetic_trans_tr =
      calculate_wilson_plaket_correlator_electric_transversal(
          plaket_tr, conf_wilson.array, T_min, T_max, R_min, R_max, d_ouside);
  end_time = omp_get_wtime();
  search_time = end_time - start_time;
  cout << "flux tube wilson transversal_tr magnetic time: " << search_time
       << endl;

  for (auto it = flux_tube_schwinger_electric_long_l.begin();
       it != flux_tube_schwinger_electric_long_l.end(); it++) {
    stream_electric_long_l
        << get<0>(it->first) << "," << get<1>(it->first) << ","
        << get<2>(it->first) << "," << it->second << ","
        << flux_tube_wilson_electric_long_l[tuple<int, int, int>(
               get<0>(it->first), get<1>(it->first), get<2>(it->first))]
        << ","
        << wilson_loops[tuple<int, int>(get<0>(it->first), get<1>(it->first))]
        << endl;
  }
  for (auto it = flux_tube_schwinger_electric_long_tr.begin();
       it != flux_tube_schwinger_electric_long_tr.end(); it++) {
    stream_electric_long_tr
        << get<0>(it->first) << "," << get<1>(it->first) << ","
        << get<2>(it->first) << "," << it->second << ","
        << flux_tube_wilson_electric_long_tr[tuple<int, int, int>(
               get<0>(it->first), get<1>(it->first), get<2>(it->first))]
        << ","
        << wilson_loops[tuple<int, int>(get<0>(it->first), get<1>(it->first))]
        << endl;
  }
  for (auto it = flux_tube_schwinger_electric_trans_l.begin();
       it != flux_tube_schwinger_electric_trans_l.end(); it++) {
    stream_electric_trans_l
        << get<0>(it->first) << "," << get<1>(it->first) << ","
        << get<2>(it->first) << "," << it->second << ","
        << flux_tube_wilson_electric_trans_l[tuple<int, int, int>(
               get<0>(it->first), get<1>(it->first), get<2>(it->first))]
        << ","
        << wilson_loops[tuple<int, int>(get<0>(it->first), get<1>(it->first))]
        << endl;
  }
  for (auto it = flux_tube_schwinger_electric_trans_tr.begin();
       it != flux_tube_schwinger_electric_trans_tr.end(); it++) {
    stream_electric_trans_tr
        << get<0>(it->first) << "," << get<1>(it->first) << ","
        << get<2>(it->first) << "," << it->second << ","
        << flux_tube_wilson_electric_trans_tr[tuple<int, int, int>(
               get<0>(it->first), get<1>(it->first), get<2>(it->first))]
        << ","
        << wilson_loops[tuple<int, int>(get<0>(it->first), get<1>(it->first))]
        << endl;
  }
  for (auto it = flux_tube_schwinger_magnetic_long_l.begin();
       it != flux_tube_schwinger_magnetic_long_l.end(); it++) {
    stream_magnetic_long_l
        << get<0>(it->first) << "," << get<1>(it->first) << ","
        << get<2>(it->first) << "," << it->second << ","
        << flux_tube_wilson_magnetic_long_l[tuple<int, int, int>(
               get<0>(it->first), get<1>(it->first), get<2>(it->first))]
        << ","
        << wilson_loops[tuple<int, int>(get<0>(it->first), get<1>(it->first))]
        << endl;
  }
  for (auto it = flux_tube_schwinger_magnetic_long_tr.begin();
       it != flux_tube_schwinger_magnetic_long_tr.end(); it++) {
    stream_magnetic_long_tr
        << get<0>(it->first) << "," << get<1>(it->first) << ","
        << get<2>(it->first) << "," << it->second << ","
        << flux_tube_wilson_magnetic_long_tr[tuple<int, int, int>(
               get<0>(it->first), get<1>(it->first), get<2>(it->first))]
        << ","
        << wilson_loops[tuple<int, int>(get<0>(it->first), get<1>(it->first))]
        << endl;
  }
  for (auto it = flux_tube_schwinger_magnetic_trans_l.begin();
       it != flux_tube_schwinger_magnetic_trans_l.end(); it++) {
    stream_magnetic_trans_l
        << get<0>(it->first) << "," << get<1>(it->first) << ","
        << get<2>(it->first) << "," << it->second << ","
        << flux_tube_wilson_magnetic_trans_l[tuple<int, int, int>(
               get<0>(it->first), get<1>(it->first), get<2>(it->first))]
        << ","
        << wilson_loops[tuple<int, int>(get<0>(it->first), get<1>(it->first))]
        << endl;
  }
  for (auto it = flux_tube_schwinger_magnetic_trans_tr.begin();
       it != flux_tube_schwinger_magnetic_trans_tr.end(); it++) {
    stream_magnetic_trans_tr
        << get<0>(it->first) << "," << get<1>(it->first) << ","
        << get<2>(it->first) << "," << it->second << ","
        << flux_tube_wilson_magnetic_trans_tr[tuple<int, int, int>(
               get<0>(it->first), get<1>(it->first), get<2>(it->first))]
        << ","
        << wilson_loops[tuple<int, int>(get<0>(it->first), get<1>(it->first))]
        << endl;
  }
}
