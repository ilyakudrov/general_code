#pragma once

#include "basic_observables.h"
#include "data.h"
#include "matrix.h"
#include "smearing.h"

#include <numeric>

// preserves data_pattern
// clockwise direction
template <class DataPattern, class MatrixType>
MatrixType
wilson_loop_schwinger(const Data::LatticeData<DataPattern, MatrixType> &conf,
                      DataPattern &data_pattern, int r, int t, int mu) {
  MatrixType A;
  for (int i = 0; i < t / 2; i++) {
    A = A * conf[data_pattern.get_index_link(3)];
    data_pattern.move_forward(1, 3);
  }
  for (int i = 0; i < r; i++) {
    A = A * conf[data_pattern.get_index_link(mu)];
    data_pattern.move_forward(1, mu);
  }
  for (int i = 0; i < t; i++) {
    data_pattern.move_backward(1, 3);
    A = A ^ conf[data_pattern.get_index_link(3)];
  }
  for (int i = 0; i < r; i++) {
    data_pattern.move_backward(1, mu);
    A = A ^ conf[data_pattern.get_index_link(mu)];
  }
  for (int i = 0; i < t / 2; i++) {
    A = A * conf[data_pattern.get_index_link(3)];
    data_pattern.move_forward(1, 3);
  }
  return A;
}

// preserves data_pattern
// counter-clockwise direction
template <class DataPattern, class MatrixType>
MatrixType wilson_loop_schwinger_opposite(
    const Data::LatticeData<DataPattern, MatrixType> &conf,
    DataPattern &data_pattern, int r, int t, int mu) {
  MatrixType A;
  for (int i = 0; i < t / 2; i++) {
    data_pattern.move_backward(1, 3);
    A = A ^ conf[data_pattern.get_index_link(3)];
  }
  for (int i = 0; i < r; i++) {
    data_pattern.move_backward(1, mu);
    A = A ^ conf[data_pattern.get_index_link(mu)];
  }
  for (int i = 0; i < t; i++) {
    A = A * conf[data_pattern.get_index_link(3)];
    data_pattern.move_forward(1, 3);
  }
  for (int i = 0; i < r; i++) {
    A = A * conf[data_pattern.get_index_link(mu)];
    data_pattern.move_forward(1, mu);
  }
  for (int i = 0; i < t / 2; i++) {
    data_pattern.move_backward(1, 3);
    A = A ^ conf[data_pattern.get_index_link(3)];
  }
  return A;
}

// preserves data_pattern
template <class DataPattern, class MatrixType>
inline MatrixType
wilson_line(const Data::LatticeData<DataPattern, MatrixType> &conf,
            DataPattern &data_pattern, int mu, int length) {
  MatrixType A;
  for (int i = 0; i < length; i++) {
    A = A * conf[data_pattern.get_index_link(mu)];
    data_pattern.move_forward(1, mu);
  }
  data_pattern.move_backward(length, mu);
  return A;
}

// calculates offaxis spatial line of wilson loop according to pattern
// pattern defines directions of links along the line
// it consists of positive and negative integers which correspond to spatial
// directions x - 1, y - 2, z - 3
template <class DataPattern, class MatrixType>
MatrixType
wilson_line_offaxis(const Data::LatticeData<DataPattern, MatrixType> &conf,
                    DataPattern &data_pattern,
                    const std::vector<int> &pattern) {
  MatrixType A;
  // iterate through pattern
  for (int i = 0; i < pattern.size(); i++) {
    // if positive direction
    if (pattern[i] > 0) {
      A = A * conf[data_pattern.get_index_link(abs(pattern[i]) - 1)];
      data_pattern.move_forward(1, abs(pattern[i]) - 1);
    }
    // if negative direction
    else if (pattern[i] < 0) {
      data_pattern.move_backward(1, abs(pattern[i]) - 1);
      A = A ^ conf[data_pattern.get_index_link(abs(pattern[i]) - 1)];
    } else {
      std::cout << "wilson_line_offaxis pattern direction error " << pattern[i]
                << std::endl;
    }
  }
  return A;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType>
wilson_lines(const Data::LatticeData<DataPattern, MatrixType> &conf, int mu,
             int length) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> wilson_lines(data_pattern.get_lattice_size());
  int place;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(place)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          place = data_pattern.get_index_site();
          wilson_lines[place] = wilson_line(conf, data_pattern, mu, length);
        }
      }
    }
  }
  return wilson_lines;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType>
wilson_line_increase(const Data::LatticeData<DataPattern, MatrixType> &conf,
                     const std::vector<MatrixType> &lines, int mu, int length) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> lines_new(data_pattern.get_lattice_size());
  int place;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(place)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          place = data_pattern.get_index_site();
          data_pattern.move_forward(length, mu);
          lines_new[place] =
              lines[place] * conf[data_pattern.get_index_link(mu)];
        }
      }
    }
  }
  return lines_new;
}

// calculate space lines in particular direction on a lattice
template <class DataPattern, class MatrixType>
std::vector<MatrixType>
wilson_lines_offaxis(const Data::LatticeData<DataPattern, MatrixType> &conf,
                     const std::vector<int> pattern) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> lines(data_pattern.get_lattice_size());
  int place;
#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(place)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          place = data_pattern.get_index_site();
          lines[place] = wilson_line_offaxis(conf, data_pattern, pattern);
        }
      }
    }
  }

  return lines;
}

// increase offaxis line by the length of the pattern
// takes direction of previous line
template <class DataPattern, class MatrixType>
std::vector<MatrixType> wilson_lines_offaxis_increase(
    const Data::LatticeData<DataPattern, MatrixType> &conf,
    const std::vector<MatrixType> &lines1, const std::vector<int> pattern,
    const std::vector<int> direction) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> lines(data_pattern.get_lattice_size());
  int place;

#pragma omp parallel for collapse(4) firstprivate(data_pattern) private(place)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          place = data_pattern.get_index_site();
          for (int i = 0; i < 3; i++) {
            data_pattern.move(direction[i], i);
          }
          lines[place] =
              lines1[place] * wilson_line_offaxis(conf, data_pattern, pattern);
        }
      }
    }
  }
  return lines;
}

// calculates wilson loop in particular direction on configuration
template <class DataPattern, class MatrixType>
double calculate_wilson_loop_offaxis(const std::vector<MatrixType> &time_lines,
                                     DataPattern &data_pattern, int time,
                                     const std::vector<MatrixType> &space_lines,
                                     const std::vector<int> &direction) {
  MatrixType A;
  double result = 0;

#pragma omp parallel for collapse(4) private(A)                                \
    firstprivate(direction, data_pattern) reduction(+ : result)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          A = time_lines[data_pattern.get_index_site()];
          data_pattern.move_forward(time, 3);
          A = A * space_lines[data_pattern.get_index_site()];
          data_pattern.move_backward(time, 3);
          for (int i = 0; i < 3; i++) {
            data_pattern.move(direction[i], i);
          }
          A = A * time_lines[data_pattern.get_index_site()].conj();
          for (int i = 0; i < 3; i++) {
            data_pattern.move(-direction[i], i);
          }
          A = A * space_lines[data_pattern.get_index_site()].conj();
          result += A.tr_real();
        }
      }
    }
  }
  return result / (data_pattern.get_lattice_size());
}

template <class DataPattern, class MatrixType>
double calculate_wilson_loop_offaxis_adjoint(
    const std::vector<MatrixType> &time_lines, DataPattern &data_pattern,
    int time, const std::vector<MatrixType> &space_lines,
    const std::vector<int> &direction) {
  MatrixType A;
  double result = 0;

#pragma omp parallel for collapse(4) private(A)                                \
    firstprivate(direction, data_pattern) reduction(+ : result)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          A = time_lines[data_pattern.get_index_site()];
          data_pattern.move_forward(time, 3);
          A = A * space_lines[data_pattern.get_index_site()];
          data_pattern.move_backward(time, 3);
          for (int i = 0; i < 3; i++) {
            data_pattern.move(direction[i], i);
          }
          A = A * time_lines[data_pattern.get_index_site()].conj();
          for (int i = 0; i < 3; i++) {
            data_pattern.move(-direction[i], i);
          }
          A = A * space_lines[data_pattern.get_index_site()].conj();
          result += A.tr_adjoint();
        }
      }
    }
  }
  return result / (data_pattern.get_lattice_size());
}

// off-axis wilson loop
// directions are for space lines of wilson loops
template <class DataPattern, class MatrixType>
std::vector<wilson_result>
wilson_offaxis(const Data::LatticeData<DataPattern, MatrixType> &conf,
               const std::vector<std::vector<int>> directions, double r_min,
               double r_max, int time_min, int time_max) {
  // std::vector of resulting wilson loops and their sizes
  std::vector<wilson_result> wilson;
  std::vector<std::vector<MatrixType>> time_lines(time_max - time_min + 1);
  std::vector<MatrixType> space_lines;
  DataPattern data_pattern(conf.lat_dim);

  // calculate time lines
  for (int i = time_min; i <= time_max; i++) {
    time_lines[i - time_min] = wilson_lines(conf, 3, i);
  }

  MatrixType A;
  double length_initial;
  std::vector<int> pattern;
  // this direction will be permutated
  std::vector<int> direction_permutated(3);

  // patterns after permutation
  std::vector<int> pattern_permutated;
  // pattern after permutation and reflection
  std::vector<int> pattern_reflected;

  // iterate through all spatial directions
  for (const auto &direction : directions) {

    // length of the initial line pattern
    length_initial =
        sqrt(direction[0] * direction[0] + direction[1] * direction[1] +
             direction[2] * direction[2]);

    // calculate pattern for this direction
    pattern = make_offaxis_pattern(direction);

    // wilson loops for one direction
    std::vector<std::vector<std::vector<double>>> wilson_tmp(
        round(r_max / length_initial - 0.5) -
            round(r_min / length_initial + 0.5) + 1,
        std::vector<std::vector<double>>(time_max - time_min + 1));

    // generate all reflections of the pattern
    std::vector<std::vector<int>> reflections;
    reflections = generate_reflections(direction);

    // generate all permutations of the pattern
    std::vector<std::vector<int>> permutations;
    permutations = generate_permutations(direction);

    // iterate through all permutations of coordinates
    for (const auto &permutation : permutations) {

      // pattern after permutation
      pattern_permutated = permutate_pattern(pattern, permutation);

      // direction after permutation
      std::vector<int> direction_permutated(3);

      // permutate direction
      for (int k = 0; k < 3; k++) {
        direction_permutated[permutation[k] - 1] = direction[k];
      }

      // iterate through all reflections of coordinates
      for (const auto &reflection : reflections) {

        // direction after permutation and reflection
        std::vector<int> direction_reflected(3);

        // reflect direction and prolong
        for (int k = 0; k < 3; k++) {
          direction_reflected[k] = direction_permutated[k] * reflection[k];
        }

        // reflect pattern
        pattern_reflected = reflect_pattern(pattern_permutated, reflection);

        // iterate through all lengths in range [r_min, r_max]
        // + 0.5 to round up, -0.5 - to round down
        for (int length_multiplier = round(r_min / length_initial + 0.5);
             length_multiplier <= round(r_max / length_initial - 0.5);
             length_multiplier++) {

          // pattern after permutation, reflection and increase
          std::vector<int> pattern_prolonged;

          // prolong pattern fo fit the length
          for (int j = 0; j < length_multiplier; j++) {
            for (int k = 0; k < pattern.size(); k++) {
              pattern_prolonged.push_back(pattern_reflected[k]);
            }
          }

          // direction after permutation, reflection and increase
          std::vector<int> direction_prolonged(3);

          for (int i = 0; i < 3; i++) {
            direction_prolonged[i] =
                direction_reflected[i] * (length_multiplier - 1);
          }

          space_lines = wilson_lines_offaxis(conf, pattern_prolonged);

          for (int i = 0; i < 3; i++) {
            direction_prolonged[i] += direction_reflected[i];
          }

          // claculate wilson loop in this direction for different
          // time sizes
          for (int time = time_min; time <= time_max; time++) {
            wilson_tmp[length_multiplier - round(r_min / length_initial + 0.5)]
                      [time - time_min]
                          .push_back(calculate_wilson_loop_offaxis(
                              time_lines[time - time_min], data_pattern, time,
                              space_lines, direction_prolonged));
          }
        }
      }
    }
    // push back the result
    wilson_result result;
    for (int length_multiplier = round(r_min / length_initial + 0.5);
         length_multiplier <= round(r_max / length_initial - 0.5);
         length_multiplier++) {

      for (int time = time_min; time <= time_max; time++) {

        result.statistics_size =
            wilson_tmp[length_multiplier - round(r_min / length_initial + 0.5)]
                      [time - time_min]
                          .size();
        result.wilson_loop =
            accumulate(
                wilson_tmp[length_multiplier -
                           round(r_min / length_initial + 0.5)][time - time_min]
                    .cbegin(),
                wilson_tmp[length_multiplier -
                           round(r_min / length_initial + 0.5)][time - time_min]
                    .cend(),
                0.0) /
            wilson_tmp[length_multiplier - round(r_min / length_initial + 0.5)]
                      [time - time_min]
                          .size();
        ;
        result.time_size = time;
        result.space_size = length_initial * length_multiplier;

        wilson.push_back(result);
      }
    }
  }
  return wilson;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, double>, double>
wilson_offaxis_result(const Data::LatticeData<DataPattern, MatrixType> &conf,
                      double r_min, double r_max, int time_min, int time_max) {
  std::vector<std::vector<int>> directions = generate_directions(4);
  std::vector<wilson_result> wilson_offaxis_result = wilson_offaxis(
      conf, directions, r_min - 0.1, r_max + 0.1, time_min, time_max);
  wilson_offaxis_reduce(wilson_offaxis_result);
  std::map<std::tuple<int, double>, double> result;
  for (int i = 0; i < wilson_offaxis_result.size(); i++) {
    result[std::tuple<int, double>(wilson_offaxis_result[i].time_size,
                                   wilson_offaxis_result[i].space_size)] =
        wilson_offaxis_result[i].wilson_loop;
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::vector<wilson_result>
wilson_offaxis_adjoint(const Data::LatticeData<DataPattern, MatrixType> &conf,
                       const std::vector<std::vector<int>> directions,
                       double r_min, double r_max, int time_min, int time_max) {
  // std::vector of resulting wilson loops and their sizes
  std::vector<wilson_result> wilson;
  std::vector<std::vector<MatrixType>> time_lines(time_max - time_min + 1);
  std::vector<MatrixType> space_lines;
  DataPattern data_pattern(conf.lat_dim);

  // calculate time lines
  for (int i = time_min; i <= time_max; i++) {
    time_lines[i - time_min] = wilson_lines(conf, 3, i);
  }

  MatrixType A;
  double length_initial;
  std::vector<int> pattern;
  // this direction will be permutated
  std::vector<int> direction_permutated(3);

  // patterns after permutation
  std::vector<int> pattern_permutated;
  // pattern after permutation and reflection
  std::vector<int> pattern_reflected;

  int test_count = 0;

  // iterate through all spatial directions
  for (const auto &direction : directions) {

    // length of the initial line pattern
    length_initial =
        sqrt(direction[0] * direction[0] + direction[1] * direction[1] +
             direction[2] * direction[2]);

    // calculate pattern for this direction
    pattern = make_offaxis_pattern(direction);

    // wilson loops for one direction
    std::vector<std::vector<std::vector<double>>> wilson_tmp(
        round(r_max / length_initial - 0.5) -
            round(r_min / length_initial + 0.5) + 1,
        std::vector<std::vector<double>>(time_max - time_min + 1));

    // generate all reflections of the pattern
    std::vector<std::vector<int>> reflections;
    reflections = generate_reflections(direction);

    // generate all permutations of the pattern
    std::vector<std::vector<int>> permutations;
    permutations = generate_permutations(direction);

    // iterate through all permutations of coordinates
    for (const auto &permutation : permutations) {

      // pattern after permutation
      pattern_permutated = permutate_pattern(pattern, permutation);

      // direction after permutation
      std::vector<int> direction_permutated(3);

      // permutate direction
      for (int k = 0; k < 3; k++) {
        direction_permutated[permutation[k] - 1] = direction[k];
      }

      // iterate through all reflections of coordinates
      for (const auto &reflection : reflections) {

        // direction after permutation and reflection
        std::vector<int> direction_reflected(3);

        // reflect direction and prolong
        for (int k = 0; k < 3; k++) {
          direction_reflected[k] = direction_permutated[k] * reflection[k];
        }

        // reflect pattern
        pattern_reflected = reflect_pattern(pattern_permutated, reflection);

        // iterate through all lengths in range [r_min, r_max]
        // + 0.5 to round up, -0.5 - to round down
        for (int length_multiplier = round(r_min / length_initial + 0.5);
             length_multiplier <= round(r_max / length_initial - 0.5);
             length_multiplier++) {

          // pattern after permutation, reflection and increase
          std::vector<int> pattern_prolonged;

          // prolong pattern fo fit the length
          for (int j = 0; j < length_multiplier; j++) {
            for (int k = 0; k < pattern.size(); k++) {
              pattern_prolonged.push_back(pattern_reflected[k]);
            }
          }

          // direction after permutation, reflection and increase
          std::vector<int> direction_prolonged(3);

          for (int i = 0; i < 3; i++) {
            direction_prolonged[i] =
                direction_reflected[i] * (length_multiplier - 1);
          }

          space_lines = wilson_lines_offaxis(conf, pattern_prolonged);

          for (int i = 0; i < 3; i++) {
            direction_prolonged[i] += direction_reflected[i];
          }

          // claculate wilson loop in this direction for different
          // time sizes
          for (int time = time_min; time <= time_max; time++) {
            test_count++;
            wilson_tmp[length_multiplier - round(r_min / length_initial + 0.5)]
                      [time - time_min]
                          .push_back(calculate_wilson_loop_offaxis_adjoint(
                              time_lines[time - time_min], data_pattern, time,
                              space_lines, direction_prolonged));
          }
        }
      }
    }
    std::cout << "test_count " << test_count << std::endl;
    // push back the result
    wilson_result result;
    for (int length_multiplier = round(r_min / length_initial + 0.5);
         length_multiplier <= round(r_max / length_initial - 0.5);
         length_multiplier++) {

      for (int time = time_min; time <= time_max; time++) {

        result.statistics_size =
            wilson_tmp[length_multiplier - round(r_min / length_initial + 0.5)]
                      [time - time_min]
                          .size();
        result.wilson_loop =
            accumulate(
                wilson_tmp[length_multiplier -
                           round(r_min / length_initial + 0.5)][time - time_min]
                    .cbegin(),
                wilson_tmp[length_multiplier -
                           round(r_min / length_initial + 0.5)][time - time_min]
                    .cend(),
                0.0) /
            wilson_tmp[length_multiplier - round(r_min / length_initial + 0.5)]
                      [time - time_min]
                          .size();
        ;
        result.time_size = time;
        result.space_size = length_initial * length_multiplier;

        wilson.push_back(result);
      }
    }
  }
  return wilson;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, double>, double> wilson_offaxis_adjoint_result(
    const Data::LatticeData<DataPattern, MatrixType> &conf, double r_min,
    double r_max, int time_min, int time_max) {
  std::vector<std::vector<int>> directions = generate_directions(4);
  std::vector<wilson_result> wilson_offaxis_result = wilson_offaxis_adjoint(
      conf, directions, r_min - 0.1, r_max + 0.1, time_min, time_max);
  wilson_offaxis_reduce(wilson_offaxis_result);
  std::map<std::tuple<int, double>, double> result;
  for (int i = 0; i < wilson_offaxis_result.size(); i++) {
    result[std::tuple<int, double>(wilson_offaxis_result[i].time_size,
                                   wilson_offaxis_result[i].space_size)] =
        wilson_offaxis_result[i].wilson_loop;
  }
  return result;
}

template <class DataPattern, class MatrixType>
void wilson_lines_prolong(
    const Data::LatticeData<DataPattern, MatrixType> &conf,
    std::vector<MatrixType> &wilson_lines, int length, int mu) {
  DataPattern data_pattern(conf.lat_dim);
  int index;
#pragma omp parallel for collapse(4) private(index)                            \
    firstprivate(data_pattern, mu)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site();
          data_pattern.move_forward(length, mu);
          wilson_lines[index] =
              wilson_lines[index] * conf[data_pattern.get_index_link(mu)];
        }
      }
    }
  }
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType> wilson_lines_get_length_one(
    const Data::LatticeData<DataPattern, MatrixType> &conf, int mu) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> wilson_lines(data_pattern.get_lattice_size());
#pragma omp parallel for collapse(4) firstprivate(mu, data_pattern)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          wilson_lines[data_pattern.get_index_site()] =
              conf[data_pattern.get_index_link(mu)];
        }
      }
    }
  }
  return wilson_lines;
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType>
wilson_lines_indexed(const Data::LatticeData<DataPattern, MatrixType> &conf,
                     int length, int mu) {
  std::vector<MatrixType> wilson_lines = wilson_lines_get_length_one(conf, mu);
  for (int i = 1; i < length; i++) {
    wilson_lines_prolong(conf, wilson_lines, i, mu);
  }
  return wilson_lines;
}

template <class DataPattern, class MatrixType>
double wilson_plane_indexed_single_rxt(
    const std::vector<MatrixType> &wilson_lines_mu,
    const std::vector<std::vector<MatrixType>> &wilson_lines_nu,
    DataPattern &data_pattern, int mu, int length_mu, int length_nu) {
  MatrixType wilson_loop;
  double result = 0;
  int index1;
  int index2;
#pragma omp parallel for collapse(3) private(wilson_loop, index1, index2)      \
    firstprivate(data_pattern, mu, length_mu, length_nu) reduction(+ : result)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index1 = data_pattern.get_index_site();
          data_pattern.move_forward(length_mu, mu);
          index2 = data_pattern.get_index_site();
          data_pattern.move_backward(length_mu, mu);
          for (int nu = 0; nu < 3; nu++) {
            wilson_loop = wilson_lines_mu[index1] * wilson_lines_nu[nu][index2];
            data_pattern.move_forward(length_nu, nu);
            wilson_loop =
                wilson_loop ^ wilson_lines_mu[data_pattern.get_index_site()];
            data_pattern.move_backward(length_nu, nu);
            result += wilson_loop.multiply_conj_tr(
                wilson_lines_nu[nu][data_pattern.get_index_site()]);
          }
        }
      }
    }
  }
  return result / (data_pattern.get_lattice_size() * 3);
}

template <class DataPattern>
double wilson_plane_indexed_single_rxt(
    const std::vector<su3> &wilson_lines_mu,
    const std::vector<std::vector<su3>> &wilson_lines_nu,
    DataPattern &data_pattern, int mu, int length_mu, int length_nu) {
  double result = 0;
  int index1;
  int index2;
  int index3;
#pragma omp parallel for collapse(3) private(index1, index2, index3)           \
    firstprivate(data_pattern, mu, length_mu, length_nu) reduction(+ : result) \
    schedule(dynamic)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index1 = data_pattern.get_index_site();
          data_pattern.move_forward(length_mu, mu);
          index2 = data_pattern.get_index_site();
          data_pattern.move_backward(length_mu, mu);
          for (int nu = 0; nu < 3; nu++) {
            data_pattern.move_forward(length_nu, nu);
            index3 = data_pattern.get_index_site();
            data_pattern.move_backward(length_nu, nu);
            result +=
                (wilson_lines_mu[index1] * wilson_lines_nu[nu][index2] *
                 wilson_lines_mu[index3].adjoint() *
                 wilson_lines_nu[nu][data_pattern.get_index_site()].adjoint())
                    .tr_real();
          }
        }
      }
    }
  }
  return result / (data_pattern.get_lattice_size() * 3);
}

template <class DataPattern, class MatrixType>
double wilson_adjoint_plane_indexed_single_rxt(
    const std::vector<MatrixType> &wilson_lines_mu,
    const std::vector<std::vector<MatrixType>> &wilson_lines_nu,
    DataPattern &data_pattern, int mu, int length_mu, int length_nu) {
  MatrixType wilson_loop;
  double result = 0;
  int index1;
  int index2;
#pragma omp parallel for collapse(4) private(wilson_loop, index1, index2)      \
    firstprivate(data_pattern, mu, length_mu, length_nu) reduction(+ : result)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index1 = data_pattern.get_index_site();
          data_pattern.move_forward(length_mu, mu);
          index2 = data_pattern.get_index_site();
          data_pattern.move_backward(length_mu, mu);
          for (int nu = 0; nu < 3; nu++) {
            wilson_loop = wilson_lines_mu[index1] * wilson_lines_nu[nu][index2];
            data_pattern.move_forward(length_nu, nu);
            wilson_loop =
                wilson_loop ^ wilson_lines_mu[data_pattern.get_index_site()];
            data_pattern.move_backward(length_nu, nu);
            result += wilson_loop.multiply_conj_tr_adjoint(
                wilson_lines_nu[nu][data_pattern.get_index_site()]);
          }
        }
      }
    }
  }
  return result / (data_pattern.get_lattice_size() * 3);
}

template <class DataPattern>
double wilson_adjoint_plane_indexed_single_rxt(
    const std::vector<su3> &wilson_lines_mu,
    const std::vector<std::vector<su3>> &wilson_lines_nu,
    DataPattern &data_pattern, int mu, int length_mu, int length_nu) {
  double result = 0;
  int index1;
  int index2;
  int index3;
#pragma omp parallel for collapse(3) private(index1, index2, index3)           \
    firstprivate(data_pattern, mu, length_mu, length_nu) reduction(+ : result) \
    schedule(dynamic)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index1 = data_pattern.get_index_site();
          data_pattern.move_forward(length_mu, mu);
          index2 = data_pattern.get_index_site();
          data_pattern.move_backward(length_mu, mu);
          for (int nu = 0; nu < 3; nu++) {
            data_pattern.move_forward(length_nu, nu);
            index3 = data_pattern.get_index_site();
            data_pattern.move_backward(length_nu, nu);
            result +=
                (wilson_lines_mu[index1] * wilson_lines_nu[nu][index2] *
                 wilson_lines_mu[index3].adjoint() *
                 wilson_lines_nu[nu][data_pattern.get_index_site()].adjoint())
                    .tr_adjoint();
          }
        }
      }
    }
  }
  return result / (data_pattern.get_lattice_size() * 3);
}

template <class DataPattern, class MatrixType>
double wilson_plane_gevp_indexed_single_rxt(
    const std::vector<MatrixType> &wilson_lines_mu,
    const std::vector<std::vector<MatrixType>> &wilson_lines_nu1,
    const std::vector<std::vector<MatrixType>> &wilson_lines_nu2,
    DataPattern &data_pattern, int mu, int length_mu, int length_nu) {
  MatrixType wilson_loop1;
  MatrixType wilson_loop2;
  double result = 0;
  int index1;
  int index2;
#pragma omp parallel for collapse(4) private(wilson_loop1, wilson_loop2,       \
                                                 index1, index2)               \
    firstprivate(data_pattern, mu, length_mu, length_nu) reduction(+ : result)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index1 = data_pattern.get_index_site();
          data_pattern.move_forward(length_mu, mu);
          index2 = data_pattern.get_index_site();
          data_pattern.move_backward(length_mu, mu);
          for (int nu = 0; nu < 3; nu++) {
            wilson_loop1 =
                wilson_lines_mu[index1] * wilson_lines_nu1[nu][index2];
            wilson_loop2 =
                wilson_lines_mu[index1] * wilson_lines_nu2[nu][index2];
            data_pattern.move_forward(length_nu, nu);
            wilson_loop1 =
                wilson_loop1 ^ wilson_lines_mu[data_pattern.get_index_site()];
            wilson_loop2 =
                wilson_loop2 ^ wilson_lines_mu[data_pattern.get_index_site()];
            data_pattern.move_backward(length_nu, nu);
            result += wilson_loop1.multiply_conj_tr(
                wilson_lines_nu2[nu][data_pattern.get_index_site()]);
            result += wilson_loop2.multiply_conj_tr(
                wilson_lines_nu1[nu][data_pattern.get_index_site()]);
          }
        }
      }
    }
  }
  return result / (data_pattern.get_lattice_size() * 6);
}

template <class DataPattern>
double wilson_plane_gevp_indexed_single_rxt(
    const std::vector<su3> &wilson_lines_mu,
    const std::vector<std::vector<su3>> &wilson_lines_nu1,
    const std::vector<std::vector<su3>> &wilson_lines_nu2,
    DataPattern &data_pattern, int mu, int length_mu, int length_nu) {
  double result = 0;
  int index1;
  int index2;
  int index3;
#pragma omp parallel for collapse(3) private(index1, index2, index3)           \
    firstprivate(data_pattern, mu, length_mu, length_nu) reduction(+ : result) \
    schedule(dynamic)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index1 = data_pattern.get_index_site();
          data_pattern.move_forward(length_mu, mu);
          index2 = data_pattern.get_index_site();
          data_pattern.move_backward(length_mu, mu);
          for (int nu = 0; nu < 3; nu++) {
            data_pattern.move_forward(length_nu, nu);
            index3 = data_pattern.get_index_site();
            data_pattern.move_backward(length_nu, nu);
            result +=
                (wilson_lines_mu[index1] * wilson_lines_nu1[nu][index2] *
                 wilson_lines_mu[index3].adjoint() *
                 wilson_lines_nu2[nu][data_pattern.get_index_site()].adjoint())
                    .tr_real();
            result +=
                (wilson_lines_mu[index1] * wilson_lines_nu2[nu][index2] *
                 wilson_lines_mu[index3].adjoint() *
                 wilson_lines_nu1[nu][data_pattern.get_index_site()].adjoint())
                    .tr_real();
          }
        }
      }
    }
  }
  return result / (data_pattern.get_lattice_size() * 6);
}

template <class DataPattern, class MatrixType>
double wilson_adjoint_plane_gevp_indexed_single_rxt(
    const std::vector<MatrixType> &wilson_lines_mu,
    const std::vector<std::vector<MatrixType>> &wilson_lines_nu1,
    const std::vector<std::vector<MatrixType>> &wilson_lines_nu2,
    DataPattern &data_pattern, int mu, int length_mu, int length_nu) {
  MatrixType wilson_loop1;
  MatrixType wilson_loop2;
  double result = 0;
  int index1;
  int index2;
#pragma omp parallel for collapse(4) private(wilson_loop1, wilson_loop2,       \
                                                 index1, index2)               \
    firstprivate(data_pattern, mu, length_mu, length_nu) reduction(+ : result)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index1 = data_pattern.get_index_site();
          data_pattern.move_forward(length_mu, mu);
          index2 = data_pattern.get_index_site();
          data_pattern.move_backward(length_mu, mu);
          for (int nu = 0; nu < 3; nu++) {
            wilson_loop1 =
                wilson_lines_mu[index1] * wilson_lines_nu1[nu][index2];
            wilson_loop2 =
                wilson_lines_mu[index1] * wilson_lines_nu2[nu][index2];
            data_pattern.move_forward(length_nu, nu);
            wilson_loop1 =
                wilson_loop1 ^ wilson_lines_mu[data_pattern.get_index_site()];
            wilson_loop2 =
                wilson_loop2 ^ wilson_lines_mu[data_pattern.get_index_site()];
            data_pattern.move_backward(length_nu, nu);
            result += wilson_loop1.multiply_conj_tr_adjoint(
                wilson_lines_nu2[nu][data_pattern.get_index_site()]);
            result += wilson_loop2.multiply_conj_tr_adjoint(
                wilson_lines_nu1[nu][data_pattern.get_index_site()]);
          }
        }
      }
    }
  }
  return result / (data_pattern.get_lattice_size() * 6);
}

template <class DataPattern>
double wilson_adjoint_plane_gevp_indexed_single_rxt(
    const std::vector<su3> &wilson_lines_mu,
    const std::vector<std::vector<su3>> &wilson_lines_nu1,
    const std::vector<std::vector<su3>> &wilson_lines_nu2,
    DataPattern &data_pattern, int mu, int length_mu, int length_nu) {
  double result = 0;
  int index1;
  int index2;
  int index3;
#pragma omp parallel for collapse(2) private(index1, index2, index3)           \
    firstprivate(data_pattern, mu, length_mu, length_nu) reduction(+ : result) \
    schedule(dynamic)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index1 = data_pattern.get_index_site();
          data_pattern.move_forward(length_mu, mu);
          index2 = data_pattern.get_index_site();
          data_pattern.move_backward(length_mu, mu);
          for (int nu = 0; nu < 3; nu++) {
            data_pattern.move_forward(length_nu, nu);
            index3 = data_pattern.get_index_site();
            data_pattern.move_backward(length_nu, nu);
            result +=
                (wilson_lines_mu[index1] * wilson_lines_nu1[nu][index2] *
                 wilson_lines_mu[index3].adjoint() *
                 wilson_lines_nu2[nu][data_pattern.get_index_site()].adjoint())
                    .tr_adjoint();
            result +=
                (wilson_lines_mu[index1] * wilson_lines_nu2[nu][index2] *
                 wilson_lines_mu[index3].adjoint() *
                 wilson_lines_nu1[nu][data_pattern.get_index_site()].adjoint())
                    .tr_adjoint();
          }
        }
      }
    }
  }
  return result / (data_pattern.get_lattice_size() * 6);
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int>, double>
wilson_loop(const Data::LatticeData<DataPattern, MatrixType> &conf, int r_min,
            int r_max, int time_min, int time_max) {
  std::map<std::tuple<int, int>, double> wilson_loops;
  std::vector<std::vector<MatrixType>> time_lines(time_max - time_min + 1);
  std::vector<std::vector<MatrixType>> space_lines(3);
  DataPattern data_pattern(conf.lat_dim);
  time_lines[0] = wilson_lines_indexed(conf, time_min, 3);
  for (int t = time_min + 1; t <= time_max; t++) {
    time_lines[t - time_min] = time_lines[t - time_min - 1];
    wilson_lines_prolong(conf, time_lines[t - time_min], t - 1, 3);
  }
  for (int mu = 0; mu < 3; mu++) {
    space_lines[mu] = wilson_lines_indexed(conf, r_min, mu);
  }
  for (int r = r_min; r <= r_max; r++) {
    for (int t = time_min; t <= time_max; t++) {
      wilson_loops[std::tuple<int, int>(t, r)] +=
          wilson_plane_indexed_single_rxt(time_lines[t - time_min], space_lines,
                                          data_pattern, 3, t, r);
    }
    for (int mu = 0; mu < 3; mu++) {
      wilson_lines_prolong(conf, space_lines[mu], r, mu);
    }
  }
  return wilson_loops;
}

// wilson loop goes in positive mu, nu directions starting from the corner
template <class DataPattern, class MatrixType>
std::vector<double>
calculate_wilson_loop_plane_tr(const std::vector<MatrixType> &wilson_lines_mu,
                               const std::vector<MatrixType> &wilson_lines_nu,
                               DataPattern &data_pattern, int mu, int nu,
                               int length_mu, int length_nu) {
  MatrixType wilson_loop;
  std::vector<double> result(data_pattern.get_lattice_size());
  int index1;
#pragma omp parallel for collapse(4) private(wilson_loop, index1)              \
    firstprivate(data_pattern, mu, nu, length_mu, length_nu)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index1 = data_pattern.get_index_site();
          data_pattern.move_forward(length_mu, mu);
          wilson_loop = wilson_lines_mu[index1] *
                        wilson_lines_nu[data_pattern.get_index_site()];
          data_pattern.move_backward(length_mu, mu);
          data_pattern.move_forward(length_nu, nu);
          wilson_loop =
              wilson_loop ^ wilson_lines_mu[data_pattern.get_index_site()];
          data_pattern.move_backward(length_nu, nu);
          result[index1] = wilson_loop.multiply_conj_tr(
              wilson_lines_nu[data_pattern.get_index_site()]);
        }
      }
    }
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::vector<double> calculate_wilson_loop_time_tr(
    const std::vector<MatrixType> &time_lines,
    const std::array<std::vector<MatrixType>, 3> &space_lines,
    DataPattern &data_pattern, int t, int r) {
  MatrixType wilson_loop;
  std::vector<double> result(data_pattern.get_lattice_size() * 3);
  int index;
#pragma omp parallel for collapse(4) private(wilson_loop, index)               \
    firstprivate(data_pattern, t, r)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site();
          for (int mu = 0; mu < 3; mu++) {
            data_pattern.move_forward(r, mu);
            wilson_loop = space_lines[mu][index] *
                          time_lines[data_pattern.get_index_site()];
            data_pattern.move_backward(r, mu);
            data_pattern.move_forward(t, 3);
            wilson_loop =
                wilson_loop ^ space_lines[mu][data_pattern.get_index_site()];
            data_pattern.move_backward(t, 3);
            result[index * 3 + mu] = wilson_loop.multiply_conj_tr(
                time_lines[data_pattern.get_index_site()]);
          }
        }
      }
    }
  }
  return result;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int>, double>
wilson_loop_adjoint(const Data::LatticeData<DataPattern, MatrixType> &conf,
                    int r_min, int r_max, int time_min, int time_max) {
  std::map<std::tuple<int, int>, double> wilson_loops;
  std::vector<std::vector<MatrixType>> time_lines(time_max - time_min + 1);
  std::vector<std::vector<MatrixType>> space_lines(3);
  DataPattern data_pattern(conf.lat_dim);
  time_lines[0] = wilson_lines_indexed(conf, time_min, 3);
  for (int t = time_min + 1; t <= time_max; t++) {
    time_lines[t - time_min] = time_lines[t - time_min - 1];
    wilson_lines_prolong(conf, time_lines[t - time_min], t - 1, 3);
  }
  for (int mu = 0; mu < 3; mu++) {
    space_lines[mu] = wilson_lines_indexed(conf, r_min, mu);
  }
  for (int r = r_min; r <= r_max; r++) {
    for (int t = time_min; t <= time_max; t++) {
      wilson_loops[std::tuple<int, int>(t, r)] +=
          wilson_adjoint_plane_indexed_single_rxt(
              time_lines[t - time_min], space_lines, data_pattern, 3, t, r);
    }
    for (int mu = 0; mu < 3; mu++) {
      wilson_lines_prolong(conf, space_lines[mu], r, mu);
    }
  }
  return wilson_loops;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int>, double>
wilson_gevp_indexed(const Data::LatticeData<DataPattern, MatrixType> &conf1,
                    const Data::LatticeData<DataPattern, MatrixType> &conf2,
                    int r_min, int r_max, int time_min, int time_max) {
  std::map<std::tuple<int, int>, double> wilson_loops;
  std::vector<std::vector<MatrixType>> time_lines(time_max - time_min + 1);
  std::vector<std::vector<MatrixType>> space_lines1(3), space_lines2(3);
  DataPattern data_pattern(conf1.lat_dim);
  time_lines[0] = wilson_lines_indexed(conf1, time_min, 3);
  for (int t = time_min + 1; t <= time_max; t++) {
    time_lines[t - time_min] = time_lines[t - time_min - 1];
    wilson_lines_prolong(conf1, time_lines[t - time_min], t - 1, 3);
  }
  for (int mu = 0; mu < 3; mu++) {
    space_lines1[mu] = wilson_lines_indexed(conf1, r_min, mu);
    space_lines2[mu] = wilson_lines_indexed(conf2, r_min, mu);
  }
  for (int r = r_min; r <= r_max; r++) {
    for (int t = time_min; t <= time_max; t++) {
      wilson_loops[std::tuple<int, int>(t, r)] +=
          wilson_plane_gevp_indexed_single_rxt(time_lines[t - time_min],
                                               space_lines1, space_lines2,
                                               data_pattern, 3, t, r);
    }
    for (int mu = 0; mu < 3; mu++) {
      wilson_lines_prolong(conf1, space_lines1[mu], r, mu);
      wilson_lines_prolong(conf2, space_lines2[mu], r, mu);
    }
  }
  return wilson_loops;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int>, double> wilson_gevp_adjoint_indexed(
    const Data::LatticeData<DataPattern, MatrixType> &conf1,
    const Data::LatticeData<DataPattern, MatrixType> &conf2, int r_min,
    int r_max, int time_min, int time_max) {
  std::map<std::tuple<int, int>, double> wilson_loops;
  std::vector<std::vector<MatrixType>> time_lines(time_max - time_min + 1);
  std::vector<std::vector<MatrixType>> space_lines1(3), space_lines2(3);
  DataPattern data_pattern(conf1.lat_dim);
  time_lines[0] = wilson_lines_indexed(conf1, time_min, 3);
  for (int t = time_min + 1; t <= time_max; t++) {
    time_lines[t - time_min] = time_lines[t - time_min - 1];
    wilson_lines_prolong(conf1, time_lines[t - time_min], t - 1, 3);
  }
  for (int mu = 0; mu < 3; mu++) {
    space_lines1[mu] = wilson_lines_indexed(conf1, r_min, mu);
    space_lines2[mu] = wilson_lines_indexed(conf2, r_min, mu);
  }
  for (int r = r_min; r <= r_max; r++) {
    for (int t = time_min; t <= time_max; t++) {
      wilson_loops[std::tuple<int, int>(t, r)] +=
          wilson_adjoint_plane_gevp_indexed_single_rxt(
              time_lines[t - time_min], space_lines1, space_lines2,
              data_pattern, 3, t, r);
    }
    for (int mu = 0; mu < 3; mu++) {
      wilson_lines_prolong(conf1, space_lines1[mu], r, mu);
      wilson_lines_prolong(conf2, space_lines2[mu], r, mu);
    }
  }
  return wilson_loops;
}

template <class DataPattern, class MatrixType>
void wilson_lines_single_direction_prolong(
    const std::vector<MatrixType> &conf, DataPattern &data_pattern,
    std::vector<MatrixType> &wilson_lines, int length, int mu) {
  int index;
#pragma omp parallel for collapse(4) private(index)                            \
    firstprivate(data_pattern, mu)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site();
          data_pattern.move_forward(length, mu);
          wilson_lines[index] =
              wilson_lines[index] * conf[data_pattern.get_index_site()];
        }
      }
    }
  }
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType>
wilson_lines_single_direction_indexed(const std::vector<MatrixType> &conf,
                                      DataPattern &data_pattern, int length,
                                      int mu) {
  std::vector<MatrixType> wilson_lines = conf;
  for (int i = 1; i < length; i++) {
    wilson_lines_single_direction_prolong(conf, data_pattern, wilson_lines,
                                          length, mu);
  }
  return wilson_lines;
}

template <class DataPattern, class MatrixType>
double
wilson_spatial_plane_indexed(const std::vector<MatrixType> &wilson_lines_mu,
                             const std::vector<MatrixType> &wilson_lines_nu,
                             DataPattern &data_pattern, int mu, int nu,
                             int length_mu, int length_nu) {
  MatrixType wilson_loop;
  double result = 0;
  int index;
#pragma omp parallel for collapse(4) private(wilson_loop, index)               \
    firstprivate(data_pattern, mu, nu, length_mu, length_nu)                   \
    reduction(+ : result)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index = data_pattern.get_index_site();
          data_pattern.move_forward(length_mu, mu);
          wilson_loop = wilson_lines_mu[index] *
                        wilson_lines_nu[data_pattern.get_index_site()];
          data_pattern.move_backward(length_mu, mu);
          data_pattern.move_forward(length_nu, nu);
          wilson_loop =
              wilson_loop ^ wilson_lines_mu[data_pattern.get_index_site()];
          result += wilson_loop.multiply_conj_tr(wilson_lines_nu[index]);
        }
      }
    }
  }
  return result / data_pattern.get_lattice_size();
}

template <class DataPattern, class MatrixType>
void wilson_spatial_3d_step_indexed(
    std::map<std::tuple<int, int, int>, double> &wilson_loops,
    const std::vector<MatrixType> &conf_nu,
    const std::vector<MatrixType> &conf_eta, DataPattern &data_pattern, int mu,
    int nu, int eta, std::vector<std::vector<MatrixType>> &quark_lines,
    int r_min, int r_max, int time_min, int time_max, int smearing) {
  std::vector<MatrixType> space_lines_nu, space_lines_eta;
  space_lines_nu = wilson_lines_single_direction_indexed(conf_nu, data_pattern,
                                                         time_min, nu);
  space_lines_eta = wilson_lines_single_direction_indexed(
      conf_eta, data_pattern, time_min, eta);
  for (int r = r_min; r <= r_max; r++) {
    for (int t = time_min; t <= time_max; t++) {
      wilson_loops[{smearing, t, r}] += wilson_spatial_plane_indexed(
          quark_lines[t - time_min], space_lines_nu, data_pattern, mu, nu, t,
          r);
      wilson_loops[{smearing, t, r}] += wilson_spatial_plane_indexed(
          quark_lines[t - time_min], space_lines_eta, data_pattern, mu, eta, t,
          r);
    }
    wilson_lines_single_direction_prolong(conf_nu, data_pattern, space_lines_nu,
                                          r, nu);
    wilson_lines_single_direction_prolong(conf_eta, data_pattern,
                                          space_lines_eta, r, eta);
  }
}

template <class DataPattern, class MatrixType>
std::vector<MatrixType> get_data_single_direction(
    const Data::LatticeData<DataPattern, MatrixType> &conf, int mu) {
  DataPattern data_pattern(conf.lat_dim);
  std::vector<MatrixType> data(data_pattern.get_lattice_size());
#pragma omp parallel for collapse(4) firstprivate(data_pattern, mu)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          data[data_pattern.get_index_site()] =
              conf[data_pattern.get_index_link(mu)];
        }
      }
    }
  }
  return data;
}

template <class DataPattern, class MatrixType>
std::map<std::tuple<int, int, int>, double> wilson_spatial_3d_indexed(
    const Data::LatticeData<DataPattern, MatrixType> &conf, int r_min,
    int r_max, int time_min, int time_max, double alpha, int smearing_start,
    int smearing_end, int smearing_step) {
  double start_time;
  double end_time;
  DataPattern data_pattern(conf.lat_dim);
  std::map<std::tuple<int, int, int>, double> wilson_loops;
  std::vector<MatrixType> space_lines;
  std::vector<std::vector<MatrixType>> quark_lines(time_max - time_min + 1);
  std::vector<MatrixType> smeared1, smeared2;

  for (int nu = 0; nu < 2; nu++) {
    for (int eta = nu + 1; eta < 3; eta++) {
      smeared1 = get_data_single_direction(conf, nu);
      smeared2 = get_data_single_direction(conf, eta);
      for (int mu = 0; mu < 3; mu++) {
        if (mu != nu && mu != eta) {
          quark_lines[0] = wilson_lines_indexed(conf, time_min, mu);
          for (int t = time_min + 1; t <= time_max; t++) {
            quark_lines[t - time_min] = quark_lines[t - time_min - 1];
            wilson_lines_prolong(conf, quark_lines[t - time_min], t - 1, mu);
          }
          wilson_spatial_3d_step_indexed(wilson_loops, smeared1, smeared2,
                                         data_pattern, mu, nu, eta, quark_lines,
                                         r_min, r_max, time_min, time_max, 0);
          for (int smearing = 1; smearing <= smearing_end; smearing++) {
            smearing_APE_2d(smeared1, smeared2, data_pattern, nu, eta, alpha);
            if ((smearing - smearing_start) % smearing_step == 0 &&
                smearing_step >= smearing_start) {
              wilson_spatial_3d_step_indexed(
                  wilson_loops, smeared1, smeared2, data_pattern, mu, nu, eta,
                  quark_lines, r_min, r_max, time_min, time_max, smearing);
            }
          }
        }
      }
    }
  }
  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    it->second = it->second / 6;
  }
  return wilson_loops;
}