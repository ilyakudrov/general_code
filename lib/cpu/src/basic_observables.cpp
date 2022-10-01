#define DATA_SIZE 4 * x_size *y_size *z_size *t_size
#define PLACE3_DIR                                                             \
  (t) * 3 * x_size *y_size *z_size + (z)*3 * x_size *y_size + (y)*3 * x_size + \
      (x)*3 + dir
#define PLACE3_LINK_NODIR                                                      \
  (link.coordinate[3]) * 3 * x_size *y_size *z_size +                          \
      (link.coordinate[2]) * 3 * x_size *y_size +                              \
      (link.coordinate[1]) * 3 * x_size + (link.coordinate[0]) * 3
#define PLACE1_NODIR                                                           \
  (t) * x_size *y_size *z_size + (z)*x_size *y_size + (y)*x_size + (x)
#define PLACE_PLAKET_TIME                                                      \
  (link.coordinate[3]) * 3 * x_size *y_size *z_size +                          \
      (link.coordinate[2]) * 3 * x_size *y_size +                              \
      (link.coordinate[1]) * 3 * x_size + (link.coordinate[0]) * 3 +           \
      link.direction
#define PLACE_PLAKET_SPACE                                                     \
  (link.coordinate[3]) * 3 * x_size *y_size *z_size +                          \
      (link.coordinate[2]) * 3 * x_size *y_size +                              \
      (link.coordinate[1]) * 3 * x_size + (link.coordinate[0]) * 3
#define PLACE1_LINK_NODIR                                                      \
  (link.coordinate[3]) * x_size *y_size *z_size +                              \
      (link.coordinate[2]) * x_size *y_size + (link.coordinate[1]) * x_size +  \
      (link.coordinate[0])

#define SPACE_ITER_START                                                       \
  for (int t = 0; t < t_size; t++) {                                           \
    for (int z = 0; z < z_size; z++) {                                         \
      for (int y = 0; y < y_size; y++) {                                       \
        for (int x = 0; x < x_size; x++) {                                     \
          link.go(x, y, z, t);                                                 \
          link.update(0);                                                      \
          link.update(1);                                                      \
          link.update(2);                                                      \
          link.update(3);

#define SPACE_ITER_END                                                         \
  }                                                                            \
  }                                                                            \
  }                                                                            \
  }

#define ITER_START_ZYX                                                         \
  for (int z = 0; z < z_size; z++) {                                           \
    for (int y = 0; y < y_size; y++) {                                         \
      for (int x = 0; x < x_size; x++) {                                       \
        link.go(x, y, z, t);                                                   \
        link.update(0);                                                        \
        link.update(1);                                                        \
        link.update(2);

#define ITER_END_3                                                             \
  }                                                                            \
  }                                                                            \
  }

#define SPACE_ITER_START_3D                                                    \
  for (int z = 0; z < z_size; z++) {                                           \
    for (int y = 0; y < y_size; y++) {                                         \
      for (int x = 0; x < x_size; x++) {                                       \
        link.go(x, y, z, t);                                                   \
        link.update(0);                                                        \
        link.update(1);                                                        \
        link.update(2);                                                        \
        link.update(3);

#define SPACE_ITER_END_3D                                                      \
  }                                                                            \
  }                                                                            \
  }

#include "../include/basic_observables.h"
#include "../include/flux_tube.h"
#include "../include/link.h"

#include <numeric>
#include <omp.h>

template <class T> double plaket_time(const std::vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  double plaket = 0;
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    SPACE_ITER_START;
    plaket += link.plaket_mu(array, 3).tr();
    SPACE_ITER_END;
  }
  return plaket / (x_size * y_size * z_size * t_size * 3);
}

template <class T> double plaket_space(const std::vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  double plaket = 0;
  SPACE_ITER_START;
  for (int mu = 0; mu < 3; mu++) {
    for (int nu = mu + 1; nu < 3; nu++) {
      link.move_dir(nu);
      plaket += link.plaket_mu(array, mu).tr();
    }
  }
  SPACE_ITER_END;
  return plaket / (x_size * y_size * z_size * t_size * 3);
}

template <class T> double plaket(const std::vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  double plaket = 0;
  SPACE_ITER_START;
  for (int mu = 0; mu < 4; mu++) {
    for (int nu = mu + 1; nu < 4; nu++) {
      link.move_dir(nu);
      plaket += link.plaket_mu(array, mu).tr();
    }
  }
  SPACE_ITER_END;
  return plaket / (x_size * y_size * z_size * t_size * 6);
}

// fast wilson_loop
template <class T>
std::vector<double> wilson(const std::vector<T> &array, int r_min, int r_max,
                           int time_min, int time_max) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<double> wilson((time_max - time_min + 1) * (r_max - r_min + 1));
  std::vector<std::vector<T>> time_lines(time_max - time_min + 1);
  std::vector<T> space_lines;
  for (int i = time_min; i <= time_max; i++) {
    time_lines[i - time_min] = wilson_lines(array, 3, i);
  }
  T A;
  for (int dir = 0; dir < 3; dir++) {
    for (int r = r_min; r <= r_max; r++) {
      // if (r == r_min)
      //   space_lines = wilson_lines(array, dir, r);
      // else
      //   space_lines = wilson_line_increase(array, space_lines, dir, r - 1);
      space_lines = wilson_lines(array, dir, r);
      for (int time = time_min; time <= time_max; time++) {

        SPACE_ITER_START

        A = time_lines[time - time_min][link.place / 4];
        link.move(3, time);
        A = A * space_lines[link.place / 4];
        link.move(3, -time);
        link.move(dir, r);
        A = A * time_lines[time - time_min][link.place / 4].conj();
        link.move(dir, -r);
        A = A * space_lines[link.place / 4].conj();

        wilson[(r - r_min) + (time - time_min) * (r_max - r_min + 1)] += A.tr();

        SPACE_ITER_END
      }
    }
  }
  for (int i = 0; i < (time_max - time_min + 1) * (r_max - r_min + 1); i++) {
    wilson[i] = wilson[i] / ((double)DATA_SIZE / 4 * 3);
  }
  return wilson;
}

template <class T>
std::vector<T> wilson_lines_single(const std::vector<T> &array, int length) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> vec(DATA_SIZE / 4);
  int place;

  for (int t = 0; t < t_size; t++) {

    SPACE_ITER_START_3D

    place = link.place;
    place = place / 4;
    vec[place] = link.wilson_line_single(array, length);

    SPACE_ITER_END_3D
  }
  return vec;
}

template <class T>
double wilson_loop_single_size(std::vector<T> lines1, std::vector<T> lines2,
                               int mu, int nu, int r1, int r2) {
  link1 link(x_size, y_size, z_size, t_size);

  T A;
  double wilson = 0;

  for (int t = 0; t < t_size; t++) {

    SPACE_ITER_START_3D

    A = lines1[link.place / 4];
    link.move(mu, r1);
    A = A * lines2[link.place / 4];
    link.move(mu, -r1);
    link.move(nu, r2);
    A = A * lines1[link.place / 4].conj();
    link.move(nu, -r2);
    A = A * lines2[link.place / 4].conj();

    wilson += A.tr();

    SPACE_ITER_END_3D
  }

  wilson = wilson / (x_size * y_size * z_size * t_size);

  return wilson;
}

// spatial wilson_loops
template <class T>
std::map<std::tuple<int, int>, double>
wilson_spatial(const std::vector<T> &array,
               std::map<std::tuple<int, int>, std::vector<T>> smeared,
               int time_min, int time_max, int r_min, int r_max) {
  link1 link(x_size, y_size, z_size, t_size);
  std::map<std::tuple<int, int>, double> wilson;
  std::vector<T> space_lines;
  std::unordered_map<int, std::vector<T>> time_lines;

  for (int i = 0; i < 3; i++) {
    for (int t = time_min; t <= time_max; t++) {
      time_lines[t] = wilson_lines(array, i, t);
    }

    for (int j = 0; j < 3; j++) {
      if (i != j) {
        for (int r = r_min; r <= r_max; r++) {
          space_lines =
              wilson_lines_single(smeared[std::tuple<int, int>{i, j}], r);

          for (int t = time_min; t <= time_max; t++) {
            wilson[std::tuple<int, int>{t, r}] +=
                wilson_loop_single_size(time_lines[t], space_lines, i, j, t, r);
          }
        }
      }
    }
  }

  for (int t = time_min; t <= time_max; t++) {
    for (int r = r_min; r <= r_max; r++) {
      wilson[std::tuple<int, int>{t, r}] =
          wilson[std::tuple<int, int>{t, r}] / 6;
    }
  }
  return wilson;
}

std::vector<double> read_abelian_fortran(std::string path_abelian) {
  std::vector<double> conf_abelian;
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  conf_abelian.reserve(data_size1);
  std::ifstream stream(path_abelian);
  std::vector<double> v(data_size1);
  stream.ignore(4);
  if (!stream.read((char *)&v[0], (data_size1) * sizeof(double)))
    std::cout << "read_abelian_fortran_test error: " << path_abelian
              << std::endl;
  for (int i = 0; i < data_size1; i++) {
    conf_abelian.push_back(v[i]);
  }
  stream.close();
  return conf_abelian;
}

std::vector<double> read_abelian_fortran_float(std::string path_abelian) {
  std::vector<double> conf_abelian;
  int data_size1 = 4 * x_size * y_size * z_size * t_size;
  conf_abelian.reserve(data_size1);
  std::ifstream stream(path_abelian);
  std::vector<float> v(data_size1);
  stream.ignore(4);
  if (!stream.read((char *)&v[0], (data_size1) * sizeof(float)))
    std::cout << "read_abelian_fortran_test error: " << path_abelian
              << std::endl;
  for (int i = 0; i < data_size1; i++) {
    conf_abelian.push_back((double)v[i]);
  }
  stream.close();
  return conf_abelian;
}

double wilson_abelian(const std::vector<double> &array, int r, int time) {
  link1 link(x_size, y_size, z_size, t_size);

  double wilson = 0;

  SPACE_ITER_START

  for (int mu = 0; mu < 3; mu++) {
    link.move_dir(mu);
    wilson += cos(link.wilson_loop_abelian(array, r, time));
  }

  SPACE_ITER_END

  wilson = wilson / (3 * x_size * y_size * z_size * t_size);
  return wilson;
}

template <class T>
std::vector<T> wilson_lines(const std::vector<T> &array, int mu, int length) {
  link1 link(x_size, y_size, z_size, t_size);
  link.move_dir(mu);
  std::vector<T> vec(DATA_SIZE / 4);
  int place;
  SPACE_ITER_START;
  place = link.place;
  place = place / 4;
  vec[place] = link.wilson_line(array, length);
  SPACE_ITER_END;
  return vec;
}

template <class T>
std::vector<T> wilson_line_increase(const std::vector<T> &array,
                                    const std::vector<T> &lines, int mu,
                                    int length) {
  link1 link(x_size, y_size, z_size, t_size);
  link.move_dir(mu);
  std::vector<T> lines_new(DATA_SIZE / 4);
  SPACE_ITER_START;
  link.move(mu, length);
  lines_new[PLACE1_NODIR] = lines[PLACE1_NODIR] * array[link.place + mu];
  SPACE_ITER_END;
  return lines_new;
}

// average over equal spacial sizes
void wilson_offaxis_reduce(std::vector<wilson_result> &wilson_offaxis_result) {
  for (int i = 0; i < wilson_offaxis_result.size(); i++) {
    for (int j = 0; j < wilson_offaxis_result.size(); j++) {
      if (wilson_offaxis_result[i].space_size ==
              wilson_offaxis_result[j].space_size &&
          wilson_offaxis_result[i].time_size ==
              wilson_offaxis_result[j].time_size &&
          i != j) {
        wilson_offaxis_result[i].wilson_loop =
            (wilson_offaxis_result[i].wilson_loop *
                 wilson_offaxis_result[i].statistics_size +
             wilson_offaxis_result[j].wilson_loop *
                 wilson_offaxis_result[j].statistics_size) /
            (wilson_offaxis_result[i].statistics_size +
             wilson_offaxis_result[j].statistics_size);
        wilson_offaxis_result.erase(wilson_offaxis_result.begin() + j);
      }
    }
  }
}

// off-axis wilson loop
// directions are for space lines of wilson loops
template <class T>
std::vector<wilson_result>
wilson_offaxis(const std::vector<T> &array,
               const std::vector<std::vector<int>> directions, double r_min,
               double r_max, int time_min, int time_max) {
  // std::vector of resulting wilson loops and their sizes
  std::vector<wilson_result> wilson;

  std::vector<std::vector<T>> time_lines(time_max - time_min + 1);
  std::vector<T> space_lines;

  // calculate time lines
  for (int i = time_min; i <= time_max; i++) {
    time_lines[i - time_min] = wilson_lines(array, 3, i);
  }

  T A;
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

          space_lines = wilson_lines_offaxis(array, pattern_prolonged);

          for (int i = 0; i < 3; i++) {
            direction_prolonged[i] += direction_reflected[i];
          }

          // claculate wilson loop in this direction for different
          // time sizes
          for (int time = time_min; time <= time_max; time++) {

            wilson_tmp[length_multiplier - round(r_min / length_initial + 0.5)]
                      [time - time_min]
                          .push_back(calculate_wilson_loop_offaxis(
                              time_lines[time - time_min], time, space_lines,
                              direction_prolonged));
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

template <class T>
std::vector<wilson_result>
wilson_offaxis_adjoint(const std::vector<T> &array,
                       const std::vector<std::vector<int>> directions,
                       double r_min, double r_max, int time_min, int time_max) {
  // std::vector of resulting wilson loops and their sizes
  std::vector<wilson_result> wilson;

  std::vector<std::vector<T>> time_lines(time_max - time_min + 1);
  std::vector<T> space_lines;

  // calculate time lines
  for (int i = time_min; i <= time_max; i++) {
    time_lines[i - time_min] = wilson_lines(array, 3, i);
  }

  T A;
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

          space_lines = wilson_lines_offaxis(array, pattern_prolonged);

          for (int i = 0; i < 3; i++) {
            direction_prolonged[i] += direction_reflected[i];
          }

          // claculate wilson loop in this direction for different
          // time sizes
          for (int time = time_min; time <= time_max; time++) {

            wilson_tmp[length_multiplier - round(r_min / length_initial + 0.5)]
                      [time - time_min]
                          .push_back(calculate_wilson_loop_offaxis_adjoint(
                              time_lines[time - time_min], time, space_lines,
                              direction_prolonged));
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

// generate all possible permutations of the  direction std::vector
std::vector<std::vector<int>>
generate_permutations(const std::vector<int> &direction) {
  std::vector<std::vector<int>> permutations;

  // check if there are two zeros in direction
  bool two_zeros = false;
  int count = 0;
  for (int i = 0; i < 3; i++) {
    if (direction[i] == 0)
      count++;
  }
  if (count == 2)
    two_zeros = true;

  // if 1 or less zeros
  if (!two_zeros) {
    for (int p1 = 1; p1 <= 3; p1++) {
      for (int p2 = 1; p2 <= 3; p2++) {
        for (int p3 = 1; p3 <= 3; p3++) {
          if (p1 != p2 && p2 != p3 && p3 != p1) {
            permutations.push_back({p1, p2, p3});
          }
        }
      }
    }
  }

  // if there are two zeros (not offaxis) - there are less permutations
  else {
    permutations.push_back({1, 2, 3});
    permutations.push_back({2, 1, 3});
    permutations.push_back({3, 1, 2});
  }

  return permutations;
}

// generate all possible reflections of the direction std::vector
std::vector<std::vector<int>>
generate_reflections(const std::vector<int> &direction) {
  std::vector<std::vector<int>> reflections;

  std::vector<int> vec1{1, -1};

  // check if there are two zeros in direction
  int count = 0;
  for (int i = 0; i < 3; i++) {
    if (direction[i] == 0)
      count++;
  }

  if (count == 0) {
    for (int r1 = 0; r1 < 2; r1++) {
      for (int r2 = 0; r2 < 2; r2++) {
        reflections.push_back({1, vec1[r1], vec1[r2]});
      }
    }
  }

  if (count == 1) {
    std::vector<int> reflection(3);
    for (int r1 = 0; r1 < 2; r1++) {
      int count1 = 0;
      while (direction[count1] == 0)
        count1++;

      // generate reflection
      for (int k = 0; k < 3; k++) {
        if (k == count1) {
          reflection[k] = vec1[r1];
        } else {
          reflection[k] = 1;
        }
      }
      reflections.push_back(reflection);
    }

  } else {
    reflections.push_back({1, 1, 1});
  }

  return reflections;
}

// calculates wilson loop in particular direction on configuration
template <class T>
double calculate_wilson_loop_offaxis(const std::vector<T> &time_lines, int time,
                                     const std::vector<T> &space_lines,
                                     const std::vector<int> &direction) {
  T A;
  link1 link(x_size, y_size, z_size, t_size);
  double result = 0;

  SPACE_ITER_START;

  A = time_lines[link.place / 4];

  link.move(3, time);

  A = A * space_lines[link.place / 4];

  link.move(3, -time);
  for (int i = 0; i < 3; i++) {
    link.move(i, direction[i]);
  }

  A = A * time_lines[link.place / 4].conj();

  for (int i = 0; i < 3; i++) {
    link.move(i, -direction[i]);
  }

  A = A * space_lines[link.place / 4].conj();

  result += A.tr();

  SPACE_ITER_END

  return result / ((double)DATA_SIZE / 4);
}

template <class T>
double calculate_wilson_loop_offaxis_adjoint(
    const std::vector<T> &time_lines, int time,
    const std::vector<T> &space_lines, const std::vector<int> &direction) {
  T A;
  link1 link(x_size, y_size, z_size, t_size);
  double result = 0;
  double trace;

  SPACE_ITER_START;

  A = time_lines[link.place / 4];

  link.move(3, time);

  A = A * space_lines[link.place / 4];

  link.move(3, -time);
  for (int i = 0; i < 3; i++) {
    link.move(i, direction[i]);
  }

  A = A * time_lines[link.place / 4].conj();

  for (int i = 0; i < 3; i++) {
    link.move(i, -direction[i]);
  }

  A = A * space_lines[link.place / 4].conj();

  trace = A.tr();
  result += trace * trace;

  SPACE_ITER_END

  return 4 * result / (x_size * y_size * z_size * t_size) - 1;
}

// calculate space lines in particular direction on a lattice
template <class T>
std::vector<T> wilson_lines_offaxis(const std::vector<T> &array,
                                    const std::vector<int> pattern) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> lines(DATA_SIZE / 4);
  int place;

  SPACE_ITER_START;

  place = link.place;
  place = place / 4;
  lines[place] = link.wilson_line_offaxis(array, pattern);

  SPACE_ITER_END;

  return lines;
}

// increase offaxis line by the length of the pattern
// takes direction of previous line
template <class T>
std::vector<T> wilson_lines_offaxis_increase(const std::vector<T> &array,
                                             const std::vector<T> &lines1,
                                             const std::vector<int> pattern,
                                             const std::vector<int> direction) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> lines(DATA_SIZE / 4);
  int place;

  SPACE_ITER_START;

  place = link.place / 4;

  for (int i = 0; i < 3; i++) {
    link.move(i, direction[i]);
  }

  lines[place] = lines1[place] * link.wilson_line_offaxis(array, pattern);

  SPACE_ITER_END;

  return lines;
}

// generate proper directions for calculating offaxis wilson loops
// slopes can not be the same
// length max is maximum length of the pattern
std::vector<std::vector<int>> generate_directions(int length_max) {
  std::vector<std::vector<int>> directions;
  std::vector<int> new_direction;
  bool new_dir;

  // iterate through all directions
  for (int a1 = 1; a1 <= length_max; a1++) {
    for (int a2 = 0; a2 <= length_max; a2++) {
      for (int a3 = 0; a3 <= length_max; a3++) {
        new_direction = {a1, a2, a3};
        if (a1 >= a2 && a2 >= a3) {
          new_dir = true;
          for (int i = 0; i < directions.size(); i++) {
            if (is_direction_same(directions[i], new_direction))
              new_dir = false;
          }
          if (new_dir) {
            directions.push_back(new_direction);
          }
        }
      }
    }
  }

  return directions;
}

// checks if two directions are the same
bool is_direction_same(const std::vector<int> &direction1,
                       const std::vector<int> &direction2) {
  return (direction1[0] * direction2[1] == direction1[1] * direction2[0] &&
          direction1[0] * direction2[2] == direction1[2] * direction2[0]);
}

// permutates directions in the pattern to change direction of the line
// replaces direction i by parmutation[i - 1] conserving sing of the
// direction
std::vector<int> permutate_pattern(const std::vector<int> &pattern,
                                   const std::vector<int> permutation) {
  std::vector<int> permutated_pattern(pattern.size());
  for (int i = 0; i < pattern.size(); i++) {
    permutated_pattern[i] =
        permutation[abs(pattern[i]) - 1] * (pattern[i] / abs(pattern[i]));
  }
  return permutated_pattern;
}

// reflects pattern according to reflection
// reflection consists of 1 and -1, which correspond to positive and
// negative directions
std::vector<int> reflect_pattern(const std::vector<int> &pattern,
                                 const std::vector<int> reflection) {
  std::vector<int> reflected_pattern(pattern.size());
  for (int i = 0; i < pattern.size(); i++) {
    reflected_pattern[i] = pattern[i] * reflection[pattern[i] - 1];
  }
  return reflected_pattern;
}

// create a pattern of offaxis line for wilson loop
// direction consists of projections of the line onto spatial directions
// pattern consists of positive integers which correspond to spatial
// directions 1, 2, 3
// projection of the line on direction 1 is the largest
// it utilizes Bresenham's line algorithm
std::vector<int> make_offaxis_pattern(const std::vector<int> &line_direction) {
  // std::vector for pattern
  std::vector<int> pattern;
  int slope1 = line_direction[1], slope2 = line_direction[2];
  int err1 = -slope1, err2 = -slope2;

  for (int i = 0; i < line_direction[0]; i++) {
    err1 += slope1 * 2;
    err2 += slope2 * 2;
    if (err1 >= line_direction[0]) {
      pattern.push_back(2);
      err1 -= 2 * line_direction[0];
    }
    if (err2 >= line_direction[0]) {
      pattern.push_back(3);
      err2 -= 2 * line_direction[0];
    }
    pattern.push_back(1);
  }
  return pattern;
}

template <class T> double polyakov_loop(const std::vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  double polyakov_loop;
  link.move_dir(3);
  SPACE_ITER_START;
  polyakov_loop += link.polyakov_loop(array).tr();
  SPACE_ITER_END;
  return polyakov_loop / (x_size * y_size * z_size * t_size);
}

template <class T>
std::vector<double> calculate_polyakov_loops_tr(const std::vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<double> polyakov_loops(x_size * y_size * z_size * t_size);
  link.move_dir(3);
  SPACE_ITER_START;

  polyakov_loops[link.place / 4] = link.polyakov_loop(array).tr();

  SPACE_ITER_END;

  return polyakov_loops;
}

template <class T>
std::map<int, double> polyakov_loop_correlator(const std::vector<T> &conf,
                                               int D_min, int D_max) {
  std::vector<double> polyakov_loops = calculate_polyakov_loops_tr(conf);
  std::map<int, double> polyakov_loop_correlator;
  link1 link(x_size, y_size, z_size, t_size);

  double polyakov_tmp;
  SPACE_ITER_START

  polyakov_tmp = polyakov_loops[link.place / 4];

  for (int mu = 0; mu < 3; mu++) {
    link.move(mu, D_min - 1);
    for (int D = D_min; D <= D_max; D++) {
      link.move(mu, 1);
      polyakov_loop_correlator[D] +=
          polyakov_tmp * polyakov_loops[link.place / 4];
    }
    link.move(mu, -D_max);
  }

  SPACE_ITER_END;

  for (auto it = polyakov_loop_correlator.begin();
       it != polyakov_loop_correlator.end(); ++it) {
    it->second = it->second / (x_size * y_size * z_size * t_size * 3);
  }

  return polyakov_loop_correlator;
}

template <class T>
std::vector<T> calculate_polyakov_loops(const std::vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> polyakov_loops(x_size * y_size * z_size * t_size);
  link.move_dir(3);
  SPACE_ITER_START;

  polyakov_loops[link.place / 4] = link.polyakov_loop(array);

  SPACE_ITER_END;

  return polyakov_loops;
}

template <class T>
std::map<int, double>
polyakov_loop_correlator_singlet(const std::vector<T> &conf, int D_min,
                                 int D_max) {
  std::vector<T> polyakov_loops = calculate_polyakov_loops(conf);
  std::map<int, double> polyakov_loop_correlator;
  link1 link(x_size, y_size, z_size, t_size);

  T polyakov_tmp;
  SPACE_ITER_START

  polyakov_tmp = polyakov_loops[link.place / 4];

  // for (int mu = 0; mu < 3; mu++) {
  //   link.move(mu, D_min);
  //   for (int D = D_min; D <= D_max; D++) {
  //     polyakov_loop_correlator[D] +=
  //         // polyakov_tmp.multiply_tr(polyakov_loops[link.place / 4]);
  //         (polyakov_tmp * polyakov_loops[link.place / 4].conj()).tr();
  //     link.move(mu, 1);
  //   }
  //   link.move(mu, -D_max - 1);
  // }

  polyakov_tmp = polyakov_loops[link.place / 4];

  // for (int mu = 0; mu < 3; mu++) {
  link.move(0, 1);
  polyakov_loop_correlator[1] +=
      (polyakov_tmp ^ polyakov_loops[link.place / 4]).tr();
  link.move(0, -1);
  // }

  SPACE_ITER_END;

  for (auto it = polyakov_loop_correlator.begin();
       it != polyakov_loop_correlator.end(); ++it) {
    it->second = it->second / (x_size * y_size * z_size * t_size /* * 3*/);
  }

  return polyakov_loop_correlator;
}

template <class T>
double plaket_plane(std::vector<T> &conf_mu, std::vector<T> &conf_nu,
                    int size_mu1, int size_mu2, int size_nu1, int size_nu2) {
  int data_size1 = x_size * y_size * z_size * t_size;

  std::vector<double> wilson_loops_tr(x_size * y_size * z_size * t_size);

  T loops;
  double result = 0;

#pragma omp parallel for collapse(3) private(loops) reduction(+ : result)
  for (int k = 0; k < data_size1; k += size_nu2) {
    for (int i = 0; i < size_nu2; i += size_mu2) {
      for (int j = 0; j < size_mu2; j++) {
        if (j < size_mu2 - size_mu1)
          loops = conf_mu[i + k + j] * conf_nu[i + k + j + size_mu1];
        else
          loops = conf_mu[i + k + j] * conf_nu[i + k + j - size_mu2 + size_mu1];
        if (i + j < size_nu2 - size_nu1)
          loops = loops ^ conf_mu[i + k + j + size_nu1];
        else
          loops = loops ^ conf_mu[i + k + j - size_nu2 + size_nu1];

        result += loops.multiply_tr(conf_nu[i + k + j]);
      }
    }
  }

  return result / data_size1;
}

template <class T>
double plaket_time_parallel(std::vector<std::vector<T>> conf) {

  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  double plaket_time = 0;

  for (int mu = 0; mu < 3; mu++) {

    plaket_time += plaket_plane(conf[mu], conf[3], steps[mu], steps[mu + 1],
                                steps[3], steps[4]);
  }

  return plaket_time / 3;
}

template <class T>
double plaket_space_parallel(std::vector<std::vector<T>> conf) {

  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  double plaket_space = 0;

  for (int mu = 0; mu < 3; mu++) {
    for (int nu = mu + 1; nu < 3; nu++) {

      plaket_space += plaket_plane(conf[mu], conf[nu], steps[mu], steps[mu + 1],
                                   steps[nu], steps[nu + 1]);
    }
  }

  return plaket_space / 3;
}

template <class T> double plaket_parallel(std::vector<std::vector<T>> conf) {

  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  double plaket = 0;

  for (int mu = 0; mu < 3; mu++) {
    for (int nu = mu + 1; nu < 4; nu++) {

      plaket += plaket_plane(conf[mu], conf[nu], steps[mu], steps[mu + 1],
                             steps[nu], steps[nu + 1]);
    }
  }

  return plaket / 6;
}

template <class T>
std::vector<std::vector<T>> separate_wilson(std::vector<T> &conf) {
  int data_size = x_size * y_size * z_size * t_size;
  std::vector<std::vector<T>> result(4, std::vector<T>(data_size));

  link1 link(x_size, y_size, z_size, t_size);

  for (int mu = 0; mu < 4; ++mu) {

    SPACE_ITER_START

    result[mu][link.place / 4] = conf[link.place + mu];

    SPACE_ITER_END
  }

  return result;
}

template <class T>
std::vector<T> merge_wilson(std::vector<std::vector<T>> &conf_separated) {
  int data_size = 4 * x_size * y_size * z_size * t_size;

  link1 link(x_size, y_size, z_size, t_size);

  std::vector<T> conf_merged(data_size);

  for (int mu = 0; mu < 4; mu++) {

    SPACE_ITER_START

    conf_merged[link.place + mu] = conf_separated[mu][link.place / 4];

    SPACE_ITER_END
  }

  return conf_merged;
}

template <class T>
std::vector<T> wilson_lines(std::vector<T> separated, int length, int size1,
                            int size2) {
  int data_size = x_size * y_size * z_size * t_size;

  std::vector<T> wilson_lines(data_size);
  T A;

  // #pragma omp parallel for collapse(2) private(A)
  for (int i = 0; i < data_size; i += size2) {
    for (int k = i; k < i + size1; k++) {
      A = separated[k];
      for (int j = k + size1; j < k + length * size1; j += size1) {
        A = A * separated[j];
      }
      wilson_lines[k] = A;
      for (int j = k + size1; j < k + size2 - (length - 1) * size1;
           j += size1) {
        A = separated[j - size1] % A;
        A = A * separated[j + (length - 1) * size1];
        wilson_lines[j] = A;
      }

      for (int j = k + size2 - (length - 1) * size1; j < k + size2;
           j += size1) {
        A = separated[j - size1] % A;
        A = A * separated[j - size2 + (length - 1) * size1];
        wilson_lines[j] = A;
      }
    }
  }

  return wilson_lines;
}

template <class T>
double wilson_plane(std::vector<T> &wilson_lines_mu,
                    std::vector<T> &wilson_lines_nu, int size_mu1, int size_mu2,
                    int size_nu1, int size_nu2, int length_mu, int length_nu) {
  int data_size = x_size * y_size * z_size * t_size;

  T loops;
  double result = 0;

#pragma omp parallel for collapse(3) private(loops) reduction(+ : result)
  for (int k = 0; k < data_size; k += size_nu2) {
    for (int i = 0; i < size_nu2; i += size_mu2) {
      for (int j = 0; j < size_mu2; j++) {
        if (j < size_mu2 - length_mu * size_mu1)
          loops = wilson_lines_mu[i + k + j] *
                  wilson_lines_nu[i + k + j + length_mu * size_mu1];
        else
          loops = wilson_lines_mu[i + k + j] *
                  wilson_lines_nu[i + k + j - size_mu2 + length_mu * size_mu1];
        if (i + j < size_nu2 - length_nu * size_nu1)
          loops = loops ^ wilson_lines_mu[i + k + j + length_nu * size_nu1];
        else
          loops = loops ^
                  wilson_lines_mu[i + k + j - size_nu2 + length_nu * size_nu1];

        result += loops.multiply_tr(wilson_lines_nu[i + k + j]);
      }
    }
  }

  return result / data_size;
}

template <class T>
double wilson_loop_test_time(std::vector<std::vector<T>> &wilson_lines,
                             int length_R, int length_T) {

  double result = 0;
  result += wilson_plane(wilson_lines[0], wilson_lines[3], 1, x_size,
                         x_size * y_size * z_size,
                         x_size * y_size * z_size * t_size, length_R, length_T);

  result += wilson_plane(wilson_lines[1], wilson_lines[3], x_size,
                         x_size * y_size, x_size * y_size * z_size,
                         x_size * y_size * z_size * t_size, length_R, length_T);

  result += wilson_plane(wilson_lines[2], wilson_lines[3], x_size * y_size,
                         x_size * y_size * z_size, x_size * y_size * z_size,
                         x_size * y_size * z_size * t_size, length_R, length_T);

  return result / 3;
}

template <class T>
std::map<std::tuple<int, int>, double>
wilson_parallel(std::vector<std::vector<T>> conf, int r_min, int r_max,
                int time_min, int time_max) {

  std::vector<int> steps = {1, x_size, x_size * y_size,
                            x_size * y_size * z_size,
                            x_size * y_size * z_size * t_size};

  std::map<std::tuple<int, int>, double> wilson_loops;

  std::vector<std::vector<T>> time_lines(time_max - time_min + 1);
  std::vector<T> space_lines;
  for (int t = time_min; t <= time_max; t++) {
    time_lines[t - time_min] = wilson_lines(conf[3], t, steps[3], steps[4]);
  }
  T A;
  for (int r = r_min; r <= r_max; r++) {
    for (int mu = 0; mu < 3; mu++) {
      space_lines = wilson_lines(conf[mu], r, steps[mu], steps[mu + 1]);
      for (int t = time_min; t <= time_max; t++) {

        wilson_loops[std::tuple<int, int>(t, r)] +=
            wilson_plane(space_lines, time_lines[t - time_min], steps[mu],
                         steps[mu + 1], steps[3], steps[4], r, t);
      }
    }
  }

  for (auto it = wilson_loops.begin(); it != wilson_loops.end(); it++) {
    it->second = it->second / 3;
  }

  return wilson_loops;
}

// su2
template double plaket_time(const std::vector<su2> &array);
template double plaket_space(const std::vector<su2> &array);
template double plaket(const std::vector<su2> &array);
template std::vector<double> wilson(const std::vector<su2> &array, int r_min,
                                    int r_max, int time_min, int time_max);
template std::vector<su2> wilson_lines(const std::vector<su2> &array, int mu,
                                       int length);
template std::vector<su2> wilson_line_increase(const std::vector<su2> &array,
                                               const std::vector<su2> &lines,
                                               int mu, int length);
template std::vector<wilson_result>
wilson_offaxis(const std::vector<su2> &array,
               const std::vector<std::vector<int>> directions, double r_min,
               double r_max, int time_min, int time_max);
template std::vector<wilson_result>
wilson_offaxis_adjoint(const std::vector<su2> &array,
                       const std::vector<std::vector<int>> directions,
                       double r_min, double r_max, int time_min, int time_max);
template double
calculate_wilson_loop_offaxis(const std::vector<su2> &time_lines, int time,
                              const std::vector<su2> &space_lines,
                              const std::vector<int> &direction);
template double calculate_wilson_loop_offaxis_adjoint(
    const std::vector<su2> &time_lines, int time,
    const std::vector<su2> &space_lines, const std::vector<int> &direction);

template std::vector<su2> wilson_lines_offaxis(const std::vector<su2> &array,
                                               const std::vector<int> pattern);
template std::vector<su2> wilson_lines_offaxis_increase(
    const std::vector<su2> &array, const std::vector<su2> &lines1,
    const std::vector<int> pattern, const std::vector<int> direction);
template double polyakov_loop(const std::vector<su2> &array);
template std::vector<double>
calculate_polyakov_loops_tr(const std::vector<su2> &array);
template std::map<int, double>
polyakov_loop_correlator(const std::vector<su2> &conf, int D_min, int D_max);
template std::vector<su2> wilson_lines_single(const std::vector<su2> &array,
                                              int length);
template std::vector<su2>
calculate_polyakov_loops(const std::vector<su2> &array);

template std::map<int, double>
polyakov_loop_correlator_singlet(const std::vector<su2> &conf, int D_min,
                                 int D_max);
template double wilson_loop_single_size(std::vector<su2> lines1,
                                        std::vector<su2> lines2, int mu, int nu,
                                        int r1, int r2);
template std::map<std::tuple<int, int>, double>
wilson_spatial(const std::vector<su2> &array,
               std::map<std::tuple<int, int>, std::vector<su2>> smeared,
               int time_min, int time_max, int r_min, int r_max);

template std::vector<std::vector<su2>> separate_wilson(std::vector<su2> &conf);

template std::vector<su2>
merge_wilson(std::vector<std::vector<su2>> &conf_separated);

template std::vector<su2> wilson_lines(std::vector<su2> separated, int length,
                                       int size1, int size2);

template double wilson_plane(std::vector<su2> &wilson_lines_mu,
                             std::vector<su2> &wilson_lines_nu, int size_mu1,
                             int size_mu2, int size_nu1, int size_nu2,
                             int length_mu, int length_nu);

template double
wilson_loop_test_time(std::vector<std::vector<su2>> &wilson_lines, int length_R,
                      int length_T);

template std::map<std::tuple<int, int>, double>
wilson_parallel(std::vector<std::vector<su2>> conf, int r_min, int r_max,
                int time_min, int time_max);

template double plaket_plane(std::vector<su2> &conf_mu,
                             std::vector<su2> &conf_nu, int size_mu1,
                             int size_mu2, int size_nu1, int size_nu2);

template double plaket_time_parallel(std::vector<std::vector<su2>> conf);

template double plaket_space_parallel(std::vector<std::vector<su2>> conf);

template double plaket_parallel(std::vector<std::vector<su2>> conf);

// abelian
template double plaket_time(const std::vector<abelian> &array);
template double plaket_space(const std::vector<abelian> &array);
template double plaket(const std::vector<abelian> &array);
template std::vector<double> wilson(const std::vector<abelian> &array,
                                    int r_min, int r_max, int time_min,
                                    int time_max);
template std::vector<abelian> wilson_lines(const std::vector<abelian> &array,
                                           int mu, int length);
template std::vector<abelian>
wilson_line_increase(const std::vector<abelian> &array,
                     const std::vector<abelian> &lines, int mu, int length);
template std::vector<wilson_result>
wilson_offaxis(const std::vector<abelian> &array,
               const std::vector<std::vector<int>> directions, double r_min,
               double r_max, int time_min, int time_max);
template std::vector<wilson_result>
wilson_offaxis_adjoint(const std::vector<abelian> &array,
                       const std::vector<std::vector<int>> directions,
                       double r_min, double r_max, int time_min, int time_max);

template double
calculate_wilson_loop_offaxis(const std::vector<abelian> &time_lines, int time,
                              const std::vector<abelian> &space_lines,
                              const std::vector<int> &direction);
template double calculate_wilson_loop_offaxis_adjoint(
    const std::vector<abelian> &time_lines, int time,
    const std::vector<abelian> &space_lines, const std::vector<int> &direction);

template std::vector<abelian>
wilson_lines_offaxis(const std::vector<abelian> &array,
                     const std::vector<int> pattern);

template std::vector<abelian> wilson_lines_offaxis_increase(
    const std::vector<abelian> &array, const std::vector<abelian> &lines1,
    const std::vector<int> pattern, const std::vector<int> direction);

template double polyakov_loop(const std::vector<abelian> &array);
template std::vector<double>
calculate_polyakov_loops_tr(const std::vector<abelian> &array);
template std::map<int, double>
polyakov_loop_correlator(const std::vector<abelian> &conf, int D_min,
                         int D_max);
template std::vector<abelian>
calculate_polyakov_loops(const std::vector<abelian> &array);
template std::map<int, double>
polyakov_loop_correlator_singlet(const std::vector<abelian> &conf, int D_min,
                                 int D_max);
template std::vector<abelian>
wilson_lines_single(const std::vector<abelian> &array, int length);
template double wilson_loop_single_size(std::vector<abelian> lines1,
                                        std::vector<abelian> lines2, int mu,
                                        int nu, int r1, int r2);
template std::map<std::tuple<int, int>, double>
wilson_spatial(const std::vector<abelian> &array,
               std::map<std::tuple<int, int>, std::vector<abelian>> smeared,
               int time_min, int time_max, int r_min, int r_max);

template std::vector<std::vector<abelian>>
separate_wilson(std::vector<abelian> &conf);

template std::vector<abelian>
merge_wilson(std::vector<std::vector<abelian>> &conf_separated);

template std::vector<abelian> wilson_lines(std::vector<abelian> separated,
                                           int length, int size1, int size2);

template double wilson_plane(std::vector<abelian> &wilson_lines_mu,
                             std::vector<abelian> &wilson_lines_nu,
                             int size_mu1, int size_mu2, int size_nu1,
                             int size_nu2, int length_mu, int length_nu);

template double
wilson_loop_test_time(std::vector<std::vector<abelian>> &wilson_lines,
                      int length_R, int length_T);

template std::map<std::tuple<int, int>, double>
wilson_parallel(std::vector<std::vector<abelian>> conf, int r_min, int r_max,
                int time_min, int time_max);

template double plaket_plane(std::vector<abelian> &conf_mu,
                             std::vector<abelian> &conf_nu, int size_mu1,
                             int size_mu2, int size_nu1, int size_nu2);

template double plaket_time_parallel(std::vector<std::vector<abelian>> conf);

template double plaket_space_parallel(std::vector<std::vector<abelian>> conf);

template double plaket_parallel(std::vector<std::vector<abelian>> conf);

// su3
template double plaket_time(const std::vector<su3> &array);
template double plaket_space(const std::vector<su3> &array);
template double plaket(const std::vector<su3> &array);
template std::vector<double> wilson(const std::vector<su3> &array, int r_min,
                                    int r_max, int time_min, int time_max);
template std::vector<su3> wilson_lines(const std::vector<su3> &array, int mu,
                                       int length);
template std::vector<su3> wilson_line_increase(const std::vector<su3> &array,
                                               const std::vector<su3> &lines,
                                               int mu, int length);
template std::vector<wilson_result>
wilson_offaxis(const std::vector<su3> &array,
               const std::vector<std::vector<int>> directions, double r_min,
               double r_max, int time_min, int time_max);
template std::vector<wilson_result>
wilson_offaxis_adjoint(const std::vector<su3> &array,
                       const std::vector<std::vector<int>> directions,
                       double r_min, double r_max, int time_min, int time_max);
template double
calculate_wilson_loop_offaxis(const std::vector<su3> &time_lines, int time,
                              const std::vector<su3> &space_lines,
                              const std::vector<int> &direction);
template double calculate_wilson_loop_offaxis_adjoint(
    const std::vector<su3> &time_lines, int time,
    const std::vector<su3> &space_lines, const std::vector<int> &direction);
template std::vector<su3> wilson_lines_offaxis(const std::vector<su3> &array,
                                               const std::vector<int> pattern);
template std::vector<su3> wilson_lines_offaxis_increase(
    const std::vector<su3> &array, const std::vector<su3> &lines1,
    const std::vector<int> pattern, const std::vector<int> direction);
template double polyakov_loop(const std::vector<su3> &array);
template std::vector<double>
calculate_polyakov_loops_tr(const std::vector<su3> &array);
template std::map<int, double>
polyakov_loop_correlator(const std::vector<su3> &conf, int D_min, int D_max);
template std::vector<su3>
calculate_polyakov_loops(const std::vector<su3> &array);
template std::map<int, double>
polyakov_loop_correlator_singlet(const std::vector<su3> &conf, int D_min,
                                 int D_max);
template std::vector<su3> wilson_lines_single(const std::vector<su3> &array,
                                              int length);
template double wilson_loop_single_size(std::vector<su3> lines1,
                                        std::vector<su3> lines2, int mu, int nu,
                                        int r1, int r2);
template std::map<std::tuple<int, int>, double>
wilson_spatial(const std::vector<su3> &array,
               std::map<std::tuple<int, int>, std::vector<su3>> smeared,
               int time_min, int time_max, int r_min, int r_max);

template std::vector<std::vector<su3>> separate_wilson(std::vector<su3> &conf);

template std::vector<su3>
merge_wilson(std::vector<std::vector<su3>> &conf_separated);

template std::vector<su3> wilson_lines(std::vector<su3> separated, int length,
                                       int size1, int size2);

template double wilson_plane(std::vector<su3> &wilson_lines_mu,
                             std::vector<su3> &wilson_lines_nu, int size_mu1,
                             int size_mu2, int size_nu1, int size_nu2,
                             int length_mu, int length_nu);

template double
wilson_loop_test_time(std::vector<std::vector<su3>> &wilson_lines, int length_R,
                      int length_T);

template std::map<std::tuple<int, int>, double>
wilson_parallel(std::vector<std::vector<su3>> conf, int r_min, int r_max,
                int time_min, int time_max);

template double plaket_plane(std::vector<su3> &conf_mu,
                             std::vector<su3> &conf_nu, int size_mu1,
                             int size_mu2, int size_nu1, int size_nu2);

template double plaket_time_parallel(std::vector<std::vector<su3>> conf);

template double plaket_space_parallel(std::vector<std::vector<su3>> conf);

template double plaket_parallel(std::vector<std::vector<su3>> conf);

// su3_abelian
template double plaket_time(const std::vector<su3_abelian> &array);
template double plaket_space(const std::vector<su3_abelian> &array);
template double plaket(const std::vector<su3_abelian> &array);
template std::vector<double> wilson(const std::vector<su3_abelian> &array,
                                    int r_min, int r_max, int time_min,
                                    int time_max);
template std::vector<su3_abelian>
wilson_lines(const std::vector<su3_abelian> &array, int mu, int length);
template std::vector<su3_abelian>
wilson_line_increase(const std::vector<su3_abelian> &array,
                     const std::vector<su3_abelian> &lines, int mu, int length);
template std::vector<wilson_result>
wilson_offaxis(const std::vector<su3_abelian> &array,
               const std::vector<std::vector<int>> directions, double r_min,
               double r_max, int time_min, int time_max);
template std::vector<wilson_result>
wilson_offaxis_adjoint(const std::vector<su3_abelian> &array,
                       const std::vector<std::vector<int>> directions,
                       double r_min, double r_max, int time_min, int time_max);
template double
calculate_wilson_loop_offaxis(const std::vector<su3_abelian> &time_lines,
                              int time,
                              const std::vector<su3_abelian> &space_lines,
                              const std::vector<int> &direction);
template double calculate_wilson_loop_offaxis_adjoint(
    const std::vector<su3_abelian> &time_lines, int time,
    const std::vector<su3_abelian> &space_lines,
    const std::vector<int> &direction);
template std::vector<su3_abelian>
wilson_lines_offaxis(const std::vector<su3_abelian> &array,
                     const std::vector<int> pattern);
template std::vector<su3_abelian>
wilson_lines_offaxis_increase(const std::vector<su3_abelian> &array,
                              const std::vector<su3_abelian> &lines1,
                              const std::vector<int> pattern,
                              const std::vector<int> direction);
template double polyakov_loop(const std::vector<su3_abelian> &array);
template std::vector<double>
calculate_polyakov_loops_tr(const std::vector<su3_abelian> &array);
template std::map<int, double>
polyakov_loop_correlator(const std::vector<su3_abelian> &conf, int D_min,
                         int D_max);
template std::vector<su3_abelian>
calculate_polyakov_loops(const std::vector<su3_abelian> &array);
template std::map<int, double>
polyakov_loop_correlator_singlet(const std::vector<su3_abelian> &conf,
                                 int D_min, int D_max);
template std::vector<su3_abelian>
wilson_lines_single(const std::vector<su3_abelian> &array, int length);
template double wilson_loop_single_size(std::vector<su3_abelian> lines1,
                                        std::vector<su3_abelian> lines2, int mu,
                                        int nu, int r1, int r2);
template std::map<std::tuple<int, int>, double>
wilson_spatial(const std::vector<su3_abelian> &array,
               std::map<std::tuple<int, int>, std::vector<su3_abelian>> smeared,
               int time_min, int time_max, int r_min, int r_max);

template std::vector<std::vector<su3_abelian>>
separate_wilson(std::vector<su3_abelian> &conf);

template std::vector<su3_abelian>
merge_wilson(std::vector<std::vector<su3_abelian>> &conf_separated);

template std::vector<su3_abelian>
wilson_lines(std::vector<su3_abelian> separated, int length, int size1,
             int size2);

template double wilson_plane(std::vector<su3_abelian> &wilson_lines_mu,
                             std::vector<su3_abelian> &wilson_lines_nu,
                             int size_mu1, int size_mu2, int size_nu1,
                             int size_nu2, int length_mu, int length_nu);

template double
wilson_loop_test_time(std::vector<std::vector<su3_abelian>> &wilson_lines,
                      int length_R, int length_T);

template std::map<std::tuple<int, int>, double>
wilson_parallel(std::vector<std::vector<su3_abelian>> conf, int r_min,
                int r_max, int time_min, int time_max);

template double plaket_plane(std::vector<su3_abelian> &conf_mu,
                             std::vector<su3_abelian> &conf_nu, int size_mu1,
                             int size_mu2, int size_nu1, int size_nu2);

template double
plaket_time_parallel(std::vector<std::vector<su3_abelian>> conf);

template double
plaket_space_parallel(std::vector<std::vector<su3_abelian>> conf);

template double plaket_parallel(std::vector<std::vector<su3_abelian>> conf);
