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

#include "../include/basic_observables.h"
#include "../include/flux_tube.h"
#include "../include/link.h"

template <class T> FLOAT plaket_time(const std::vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  result res(0);
  FLOAT aver[2];
  for (int dir = 0; dir < 3; dir++) {
    link.move_dir(dir);
    SPACE_ITER_START;
    res.array.push_back(link.plaket_mu(array, 3).tr());
    SPACE_ITER_END;
  }
  res.average(aver);
  return aver[0];
}

template <class T> FLOAT plaket_space(const std::vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  result res(0);
  FLOAT aver[2];
  SPACE_ITER_START;
  for (int mu = 0; mu < 3; mu++) {
    for (int nu = mu + 1; nu < 3; nu++) {
      link.move_dir(nu);
      res.array.push_back(link.plaket_mu(array, mu).tr());
    }
  }
  SPACE_ITER_END;
  res.average(aver);
  return aver[0];
}

template <class T> FLOAT plaket(const std::vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  result res(DATA_SIZE / 4 * 6);
  FLOAT aver[2];
  SPACE_ITER_START;
  for (int mu = 0; mu < 4; mu++) {
    for (int nu = mu + 1; nu < 4; nu++) {
      link.move_dir(nu);
      res.array.push_back(link.plaket_mu(array, mu).tr());
    }
  }
  SPACE_ITER_END;
  res.average(aver);
  return aver[0];
}

// fast wilson_loop
template <class T>
std::vector<FLOAT> wilson(const std::vector<T> &array, int r_min, int r_max,
                          int time_min, int time_max) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<FLOAT> wilson((time_max - time_min + 1) * (r_max - r_min + 1));
  std::vector<std::vector<T>> time_lines(time_max - time_min + 1);
  std::vector<T> space_lines;
  for (int i = time_min; i <= time_max; i++) {
    time_lines[i - time_min] = wilson_lines(array, 3, i);
  }
  T A;
  for (int dir = 0; dir < 3; dir++) {
    for (int r = r_min; r <= r_max; r++) {
      if (r == r_min)
        space_lines = wilson_lines(array, dir, r);
      else
        space_lines = wilson_line_increase(array, space_lines, dir, r - 1);
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
    wilson[i] = wilson[i] / (DATA_SIZE / 4 * 3);
  }
  return wilson;
}

template <class T>
std::vector<T> wilson_lines(const std::vector<T> &array, int mu, int length) {
  link1 link(x_size, y_size, z_size, t_size);
  link.move_dir(mu);
  std::vector<T> vec;
  vec.reserve(DATA_SIZE / 4);
  SPACE_ITER_START;
  vec[link.place / 4] = link.wilson_line(array, length);
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
  lines_new[PLACE1_NODIR] = lines[PLACE1_NODIR] * link.get_matrix(array);
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
               const std::vector<std::vector<int>> directions, FLOAT r_min,
               FLOAT r_max, int time_min, int time_max) {
  // std::vector of resulting wilson loops and their sizes
  std::vector<wilson_result> wilson;

  std::vector<std::vector<T>> time_lines(time_max - time_min + 1);
  std::vector<T> space_lines;

  // calculate time lines
  for (int i = time_min; i <= time_max; i++) {
    time_lines[i - time_min] = wilson_lines(array, 3, i);
  }

  T A;
  FLOAT length_initial;
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
    std::vector<std::vector<result>> wilson_tmp(
        round(r_max / length_initial - 0.5) -
            round(r_min / length_initial + 0.5) + 1,
        std::vector<result>(time_max - time_min + 1));

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

          // TODO: fix increase of wilson line and check if it reduces time of
          // calculations
          // calculate space line in this direction
          // if (round(r_min / length_initial + 0.5) == length_multiplier) {
          //   space_lines = wilson_lines_offaxis(array, pattern_prolonged);
          // } else {
          //   space_lines = wilson_lines_offaxis_increase(
          //       array, space_lines, pattern, direction_prolonged);
          // }
          space_lines = wilson_lines_offaxis(array, pattern_prolonged);

          for (int i = 0; i < 3; i++) {
            direction_prolonged[i] += direction_reflected[i];
          }

          // claculate wilson loop in this direction for different
          // time sizes
          for (int time = time_min; time <= time_max; time++) {

            wilson_tmp[length_multiplier - round(r_min / length_initial + 0.5)]
                      [time - time_min]
                          .array.push_back(calculate_wilson_loop_offaxis(
                              time_lines[time - time_min], time, space_lines,
                              direction_prolonged));
          }
        }
      }
    }
    // push back the result
    FLOAT aver[2];
    wilson_result result;
    for (int length_multiplier = round(r_min / length_initial + 0.5);
         length_multiplier <= round(r_max / length_initial - 0.5);
         length_multiplier++) {

      for (int time = time_min; time <= time_max; time++) {

        wilson_tmp[length_multiplier - round(r_min / length_initial + 0.5)]
                  [time - time_min]
                      .average(aver);
        result.statistics_size =
            wilson_tmp[length_multiplier - round(r_min / length_initial + 0.5)]
                      [time - time_min]
                          .array.size();
        result.wilson_loop = aver[0];
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
FLOAT calculate_wilson_loop_offaxis(const std::vector<T> &time_lines, int time,
                                    const std::vector<T> &space_lines,
                                    const std::vector<int> &direction) {
  T A;
  link1 link(x_size, y_size, z_size, t_size);
  FLOAT result = 0;

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

  return result / (DATA_SIZE / 4);
}

// calculate space lines in particular direction on a lattice
template <class T>
std::vector<T> wilson_lines_offaxis(const std::vector<T> &array,
                                    const std::vector<int> pattern) {
  link1 link(x_size, y_size, z_size, t_size);
  std::vector<T> lines;
  lines.reserve(DATA_SIZE / 4);

  SPACE_ITER_START;

  lines[link.place / 4] = link.wilson_line_offaxis(array, pattern);

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
  std::vector<T> lines;
  lines.reserve(DATA_SIZE / 4);
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

// chaecks if two directions are the same
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

template <class T> FLOAT polyakov(const std::vector<T> &array) {
  link1 link(x_size, y_size, z_size, t_size);
  result res(0);
  FLOAT aver[2];
  link.move_dir(3);
  SPACE_ITER_START;
  res.array.push_back(link.polyakov_loop(array).tr());
  SPACE_ITER_END;
  res.average(aver);
  return aver[0];
}

template <class T>
FLOAT polyakov_loop_corelator(const std::vector<T> &array, int D) {
  std::vector<T> polyakov_loop = calculate_polyakov_loop(array);
  link1 link(x_size, y_size, z_size, t_size);
  result vec(0);
  FLOAT aver[2];
  FLOAT a;
  for (int mu = 0; mu < 3; mu++) {
    SPACE_ITER_START;
    a = polyakov_loop[PLACE1_LINK_NODIR].tr();
    link.move(mu, D);
    a *= polyakov_loop[PLACE1_LINK_NODIR].conj().tr();
    vec.array.push_back(a);
    SPACE_ITER_END;
  }
  vec.average(aver);
  return aver[0];
}

FLOAT MAG_functional_su2(const std::vector<su2> &array) {
  FLOAT result = 0;
  for (int i = 0; i < DATA_SIZE; i++) {
    result += (array[i].sigma3_mult() * array[i].conj()).tr();
  }
  return result / (2 * DATA_SIZE);
}

// su2
template FLOAT plaket_time(const std::vector<su2> &array);
template FLOAT plaket_space(const std::vector<su2> &array);
template FLOAT plaket(const std::vector<su2> &array);
template std::vector<FLOAT> wilson(const std::vector<su2> &array, int r_min,
                                   int r_max, int time_min, int time_max);
template std::vector<su2> wilson_lines(const std::vector<su2> &array, int mu,
                                       int length);
template std::vector<su2> wilson_line_increase(const std::vector<su2> &array,
                                               const std::vector<su2> &lines,
                                               int mu, int length);
template std::vector<wilson_result>
wilson_offaxis(const std::vector<su2> &array,
               const std::vector<std::vector<int>> directions, FLOAT r_min,
               FLOAT r_max, int time_min, int time_max);
template FLOAT
calculate_wilson_loop_offaxis(const std::vector<su2> &time_lines, int time,
                              const std::vector<su2> &space_lines,
                              const std::vector<int> &direction);
template std::vector<su2> wilson_lines_offaxis(const std::vector<su2> &array,
                                               const std::vector<int> pattern);
template std::vector<su2> wilson_lines_offaxis_increase(
    const std::vector<su2> &array, const std::vector<su2> &lines1,
    const std::vector<int> pattern, const std::vector<int> direction);
template FLOAT polyakov(const std::vector<su2> &array);
template FLOAT polyakov_loop_corelator(const std::vector<su2> &array, int D);

// abelian
template FLOAT plaket_time(const std::vector<abelian> &array);
template FLOAT plaket_space(const std::vector<abelian> &array);
template FLOAT plaket(const std::vector<abelian> &array);
template std::vector<FLOAT> wilson(const std::vector<abelian> &array, int r_min,
                                   int r_max, int time_min, int time_max);
template std::vector<abelian> wilson_lines(const std::vector<abelian> &array,
                                           int mu, int length);
template std::vector<abelian>
wilson_line_increase(const std::vector<abelian> &array,
                     const std::vector<abelian> &lines, int mu, int length);
template std::vector<wilson_result>
wilson_offaxis(const std::vector<abelian> &array,
               const std::vector<std::vector<int>> directions, FLOAT r_min,
               FLOAT r_max, int time_min, int time_max);

template FLOAT
calculate_wilson_loop_offaxis(const std::vector<abelian> &time_lines, int time,
                              const std::vector<abelian> &space_lines,
                              const std::vector<int> &direction);

template std::vector<abelian>
wilson_lines_offaxis(const std::vector<abelian> &array,
                     const std::vector<int> pattern);

template std::vector<abelian> wilson_lines_offaxis_increase(
    const std::vector<abelian> &array, const std::vector<abelian> &lines1,
    const std::vector<int> pattern, const std::vector<int> direction);

template FLOAT polyakov(const std::vector<abelian> &array);
template FLOAT polyakov_loop_corelator(const std::vector<abelian> &array,
                                       int D);