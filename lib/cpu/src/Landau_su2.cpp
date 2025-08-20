#include "../include/data.h"
#include "../include/matrix.h"

#include <cmath>
#include <map>
#include <omp.h>
#include <random>
#include <vector>

double Landau_functional_conf(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf,
    const std::vector<su2> &gauge) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  double result = 0;
  int index_gauge;
  int index_conf;
#pragma omp parallel for collapse(4) private(index_gauge, index_conf)          \
    firstprivate(data_pattern) reduction(+ : result)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          index_gauge = data_pattern.get_index_site();
          for (int mu = 0; mu < 4; mu++) {
            index_conf = data_pattern.get_index_link(mu);
            data_pattern.move_forward(1, mu);
            result += ((gauge[index_gauge] * conf[index_conf]) ^
                       gauge[data_pattern.get_index_site()])
                          .tr();
            data_pattern.move_backward(1, mu);
          }
        }
      }
    }
  }
  return result / (data_pattern.get_data_size());
}

double Landau_su2_functional(
    Data::LatticeData<DataPatternLexicographical, su2> &conf_su2) {
  double result = 0;
  for (int i = 0; i < conf_su2.array.size(); i++) {
    result += conf_su2[i].tr();
  }
  return result / conf_su2.array.size();
}

std::vector<su2>
generate_gauge_su2_uniform(DataPatternLexicographical &data_pattern) {
  unsigned seed = time(NULL);
  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);
  int data_size = data_pattern.get_lattice_size();
  std::vector<su2> gauge_complex(data_size);
  double x1, x2, x3, x4;
  double a, b;
  for (int i = 0; i < data_size; i++) {
    do {
      x1 = (double)random_generator() / (random_generator.max() + 1);
      x2 = (double)random_generator() / (random_generator.max() + 1);
      a = x1 * x1 + x2 * x2;
    } while (a >= 1.0);
    do {
      x3 = (double)random_generator() / (random_generator.max() + 1);
      x4 = (double)random_generator() / (random_generator.max() + 1);
      b = x3 * x3 + x4 * x4;
    } while (b >= 1.0);
    gauge_complex[i] =
        su2(x1, x2, x3 * sqrt((1 - a) / b), x4 * sqrt((1 - a) / b));
  }
  return gauge_complex;
}

void apply_gauge_Landau(
    std::vector<su2> &gauge,
    Data::LatticeData<DataPatternLexicographical, su2> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  su2 tmp;
  double norm;
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          for (int mu = 0; mu < 4; mu++) {
            tmp = gauge[data_pattern.get_index_site()] *
                  conf[data_pattern.get_index_link(mu)];
            data_pattern.move_forward(1, mu);
            tmp = tmp ^ gauge[data_pattern.get_index_site()];
            data_pattern.move_backward(1, mu);
            conf[data_pattern.get_index_link(mu)] = tmp;
          }
        }
      }
    }
  }
}

inline void contribution_conj(su2 &contribution, const su2 &A, const su2 &B) {
  contribution.a0 += A.a0 * B.a0 + A.a1 * B.a1 + A.a2 * B.a2 + A.a3 * B.a3;
  contribution.a1 += -A.a0 * B.a1 + B.a0 * A.a1 - A.a3 * B.a2 + A.a2 * B.a3;
  contribution.a2 += -A.a0 * B.a2 + B.a0 * A.a2 - A.a1 * B.a3 + A.a3 * B.a1;
  contribution.a3 += -A.a0 * B.a3 + B.a0 * A.a3 - A.a2 * B.a1 + A.a1 * B.a2;
}

inline void contribution_conj_conj(su2 &contribution, const su2 &A,
                                   const su2 &B) {
  contribution.a0 += A.a0 * B.a0 - A.a1 * B.a1 - A.a2 * B.a2 - A.a3 * B.a3;
  contribution.a1 += A.a3 * B.a2 - B.a3 * A.a2 - A.a1 * B.a0 - A.a0 * B.a1;
  contribution.a2 += A.a1 * B.a3 - B.a1 * A.a3 - A.a2 * B.a0 - A.a0 * B.a2;
  contribution.a3 += A.a2 * B.a1 - B.a2 * A.a1 - A.a3 * B.a0 - A.a0 * B.a3;
}

su2 contribution_site(
    std::vector<su2> &gauge,
    const Data::LatticeData<DataPatternLexicographical, su2> &conf,
    DataPatternLexicographical &data_pattern, int position,
    std::vector<int> &shift) {
  su2 contribution(0, 0, 0, 0);
  // mu = 0
  if (data_pattern.lat_coord[0] < data_pattern.lat_dim[0] - 1) {
    contribution_conj(contribution, conf[position * 4],
                      gauge[position + shift[0]]);
  } else
    contribution_conj(contribution, conf[position * 4],
                      gauge[position + shift[0] - shift[1]]);
  if (data_pattern.lat_coord[0] > 0)
    contribution_conj_conj(contribution, conf[(position - shift[0]) * 4],
                           gauge[position - shift[0]]);
  else
    contribution_conj_conj(contribution,
                           conf[(position - shift[0] + shift[1]) * 4],
                           gauge[position - shift[0] + shift[1]]);
  // mu = 1
  if (data_pattern.lat_coord[1] < data_pattern.lat_dim[1] - 1)
    contribution_conj(contribution, conf[position * 4 + 1],
                      gauge[position + shift[1]]);
  else
    contribution_conj(contribution, conf[position * 4 + 1],
                      gauge[position + shift[1] - shift[2]]);
  if (data_pattern.lat_coord[1] > 0)
    contribution_conj_conj(contribution, conf[(position - shift[1]) * 4 + 1],
                           gauge[position - shift[1]]);
  else
    contribution_conj_conj(contribution,
                           conf[(position - shift[1] + shift[2]) * 4 + 1],
                           gauge[position - shift[1] + shift[2]]);
  // mu = 2
  if (data_pattern.lat_coord[2] < data_pattern.lat_dim[2] - 1)
    contribution_conj(contribution, conf[position * 4 + 2],
                      gauge[position + shift[2]]);
  else
    contribution_conj(contribution, conf[position * 4 + 2],
                      gauge[position + shift[2] - shift[3]]);
  if (data_pattern.lat_coord[2] > 0)
    contribution_conj_conj(contribution, conf[(position - shift[2]) * 4 + 2],
                           gauge[position - shift[2]]);
  else
    contribution_conj_conj(contribution,
                           conf[(position - shift[2] + shift[3]) * 4 + 2],
                           gauge[position - shift[2] + shift[3]]);
  // mu = 3
  if (data_pattern.lat_coord[3] < data_pattern.lat_dim[3] - 1)
    contribution_conj(contribution, conf[position * 4 + 3],
                      gauge[position + shift[3]]);
  else
    contribution_conj(contribution, conf[position * 4 + 3],
                      gauge[position + shift[3] - shift[4]]);
  if (data_pattern.lat_coord[3] > 0)
    contribution_conj_conj(contribution, conf[(position - shift[3]) * 4 + 3],
                           gauge[position - shift[3]]);
  else
    contribution_conj_conj(contribution,
                           conf[(position - shift[3] + shift[4]) * 4 + 3],
                           gauge[position - shift[3] + shift[4]]);
  return contribution;
}

void heat_bath(
    su2 &gauge, su2 &neighbour, double &temperature,
    std::subtract_with_carry_engine<unsigned, 24, 10, 24> &random_generator) {
  double neighbour_det_sq = sqrt(neighbour.determinant());
  neighbour.a0 /= neighbour_det_sq;
  neighbour.a1 /= neighbour_det_sq;
  neighbour.a2 /= neighbour_det_sq;
  neighbour.a3 /= neighbour_det_sq;
  double r1, r2, r3;
  double lambda_sq;
  double a;
  do {
    r1 = (double)random_generator() / (random_generator.max() + 1);
    r2 = (double)random_generator() / (random_generator.max() + 1);
    r3 = (double)random_generator() / (random_generator.max() + 1);
    a = cos(2 * M_PI * r2);
    lambda_sq =
        -(log(r1) + a * a * log(r3)) * temperature / 4 / neighbour_det_sq;
    r1 = (double)random_generator() / (random_generator.max() + 1);
  } while (r1 * r1 > 1 - lambda_sq);
  gauge.a0 = 1 - 2 * lambda_sq;
  double x1, x2, b;
  double c = 2. / random_generator.max();
  do {
    x1 = c * random_generator() - 1;
    x2 = c * random_generator() - 1;
    b = x1 * x1 + x2 * x2;
  } while (b >= 1);
  a = sqrt(1 - b);
  c = sqrt(1 - gauge.a0 * gauge.a0);
  gauge.a1 = 2 * x1 * a * c;
  gauge.a2 = 2 * x2 * a * c;
  gauge.a3 = (1 - 2 * b) * c;
  gauge = gauge ^ neighbour;
}

void heat_bath_update(
    std::vector<su2> &gauge,
    const Data::LatticeData<DataPatternLexicographical, su2> &conf,
    double temperature) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  std::vector<int> shift = {1, data_pattern.lat_dim[0],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2] *
                                data_pattern.lat_dim[3]};
  unsigned seed = time(NULL);
  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);
  su2 contribution;
  int position;
#pragma omp parallel for collapse(4) private(contribution, position)           \
    firstprivate(data_pattern, temperature, random_generator, shift)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          contribution =
              contribution_site(gauge, conf, data_pattern, position, shift);
          heat_bath(gauge[position], contribution, temperature,
                    random_generator);
        }
      }
    }
  }
}

void overrelaxation_update(
    std::vector<su2> &gauge,
    const Data::LatticeData<DataPatternLexicographical, su2> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  std::vector<int> shift = {1, data_pattern.lat_dim[0],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2] *
                                data_pattern.lat_dim[3]};
  int position;
  su2 contribution;
  double contribution_det_sq;
#pragma omp parallel for collapse(4) private(contribution, position,           \
                                                 contribution_det_sq)          \
    firstprivate(data_pattern, shift)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          contribution =
              contribution_site(gauge, conf, data_pattern, position, shift);
          contribution_det_sq = sqrt(contribution.determinant());
          contribution.a0 /= contribution_det_sq;
          contribution.a1 /= contribution_det_sq;
          contribution.a2 /= contribution_det_sq;
          contribution.a3 /= contribution_det_sq;
          gauge[position] =
              (contribution.conj() ^ gauge[position]) ^ contribution;
        }
      }
    }
  }
}

std::tuple<double, double> relaxation_update(
    std::vector<su2> &gauge,
    const Data::LatticeData<DataPatternLexicographical, su2> &conf) {
  DataPatternLexicographical data_pattern(conf.lat_dim);
  std::vector<int> shift = {1, data_pattern.lat_dim[0],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2],
                            data_pattern.lat_dim[0] * data_pattern.lat_dim[1] *
                                data_pattern.lat_dim[2] *
                                data_pattern.lat_dim[3]};
  int position;
  su2 contribution;
  double contribution_det_sq;
  double diff_average = 0;
  double diff_maximal = 0;
  double diff_tmp;
#pragma omp parallel for collapse(4) private(contribution, position, diff_tmp, \
                                                 contribution_det_sq)          \
    firstprivate(data_pattern, shift) reduction(+ : diff_average)              \
    reduction(max : diff_maximal)
  for (int t = 0; t < data_pattern.lat_dim[3]; t++) {
    for (int z = 0; z < data_pattern.lat_dim[2]; z++) {
      for (int y = 0; y < data_pattern.lat_dim[1]; y++) {
        for (int x = 0; x < data_pattern.lat_dim[0]; x++) {
          data_pattern.lat_coord = {x, y, z, t};
          position = data_pattern.get_index_site();
          contribution =
              contribution_site(gauge, conf, data_pattern, position, shift);
          contribution_det_sq = sqrt(contribution.determinant());
          contribution.a0 /= contribution_det_sq;
          contribution.a1 /= contribution_det_sq;
          contribution.a2 /= contribution_det_sq;
          contribution.a3 /= contribution_det_sq;
          diff_tmp = 1 - gauge[position].a0 * contribution.a0 +
                     gauge[position].a1 * contribution.a1 +
                     gauge[position].a2 * contribution.a2 +
                     gauge[position].a3 * contribution.a3;
          if (diff_tmp > diff_maximal)
            diff_maximal = diff_tmp;
          diff_average += diff_tmp;
          gauge[position].a0 = contribution.a0;
          gauge[position].a1 = -contribution.a1;
          gauge[position].a2 = -contribution.a2;
          gauge[position].a3 = -contribution.a3;
        }
      }
    }
  }
  diff_average = diff_average / data_pattern.get_lattice_size();
  return std::tuple<double, double>(diff_maximal, diff_average);
}

std::map<double, double> simulated_annealing_thermalization_test(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf,
    std::vector<su2> &gauge, double T_init, double T_final, double T_step,
    int OR_steps, int thermalization_steps, int local_thermalization_steps) {
  double omp_time;
  std::map<double, double> result;
  for (int i = 0; i < thermalization_steps; i++) {
    heat_bath_update(gauge, conf, T_init);
  }
  double T = T_init;
  while (T > T_final / 2) {
    omp_time = omp_get_wtime();
    for (int i = 0; i < local_thermalization_steps; i++) {
      heat_bath_update(gauge, conf, T);
      for (int i = 0; i < OR_steps; i++) {
        overrelaxation_update(gauge, conf);
      }
    }
    std::cout << "heat_bath_update time: " << omp_get_wtime() - omp_time
              << std::endl;
    result[T] = Landau_functional_conf(conf, gauge);
    std::cout << "T: " << T << " functional: " << result[T] << std::endl;
    if (T <= 3 + T_step && T >= 2.3)
      T -= T_step / 10;
    else
      T -= T_step;
  }
  return result;
}

void make_simulated_annealing(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf,
    std::vector<su2> &gauge, double T_init, double T_final, double T_step,
    int OR_steps, int thermalization_steps) {
  for (int i = 0; i < thermalization_steps; i++) {
    heat_bath_update(gauge, conf, T_init);
  }
  double T = T_init;
  while (T > T_final / 2) {
    heat_bath_update(gauge, conf, T);
    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(gauge, conf);
    }
    if (T <= 3 + T_step && T >= 2.3)
      T -= T_step / 10;
    else
      T -= T_step;
  }
}

void make_maximization_final(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf,
    std::vector<su2> &gauge, int OR_steps, double tolerance_maximal,
    double tolerance_average) {
  bool is_converged = false;
  std::tuple<double, double> difference;
  do {
    difference = relaxation_update(gauge, conf);
    // std::cout << "functional before: " << Landau_functional_conf(conf, gauge)
    //           << std::endl;
    // std::cout << "difference: " << std::get<0>(difference) << " "
    //           << std::get<1>(difference) << std::endl;
    is_converged = (std::get<0>(difference) < tolerance_maximal) &&
                   (std::get<1>(difference) < tolerance_average);
    for (int i = 0; i < OR_steps; i++) {
      overrelaxation_update(gauge, conf);
    }
  } while (!is_converged);
}
