#include "../include/Landau_gauge.h"
#include "../include/link.h"

#include <math.h>
#include <random>
#include <vector>

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

using namespace std;

double functional_Landau(vector<su2> &conf) {
  double functional = 0;
  for (auto A : conf) {
    functional += A.tr();
  }
  return functional / (x_size * y_size * z_size * t_size * 4);
}

double functional_Landau(vector<abelian> &conf) {
  double functional = 0;
  for (auto A : conf) {
    functional += A.tr();
  }
  return functional / (x_size * y_size * z_size * t_size * 4);
}

double functional_Landau(vector<abelian> &angles_conf,
                         vector<double> &angles_gauge) {

  link1 link(x_size, y_size, z_size, t_size);
  double functional = 0;

  abelian tmp;
  abelian tmp_site;

  SPACE_ITER_START

  tmp_site = abelian(1, angles_gauge[link.place / 4]);

  for (int mu = 0; mu < 4; mu++) {

    tmp = tmp_site * angles_conf[link.place + mu];

    link.move(mu, 1);

    tmp = tmp * abelian(1, -angles_gauge[link.place / 4]);
    functional += tmp.tr();

    link.move(mu, -1);
  }

  SPACE_ITER_END

  return functional / (x_size * y_size * z_size * t_size * 4);
}

vector<abelian> convert_su2_abelian(vector<su2> &conf_su2) {
  vector<abelian> angles_conf(x_size * y_size * z_size * t_size * 4);

  for (int i = 0; i < conf_su2.size(); i++) {
    angles_conf[i] = abelian(
        sqrt(conf_su2[i].a3 * conf_su2[i].a3 + conf_su2[i].a0 * conf_su2[i].a0),
        atan2(conf_su2[i].a3, conf_su2[i].a0));
  }

  return angles_conf;
}

vector<double> generate_angles_gauge_random() {
  vector<double> angles(x_size * y_size * z_size * t_size);

  unsigned seed = time(NULL);
  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);

  for (int i = 0; i < angles.size(); i++) {
    angles[i] = (double)random_generator() / random_generator.max() * 2 * M_PI;
  }

  return angles;
}

vector<double> generate_angles_gauge_trivial() {
  vector<double> angles(x_size * y_size * z_size * t_size);

  for (int i = 0; i < angles.size(); i++) {
    angles[i] = 0;
  }

  return angles;
}

abelian angle_neighbour(vector<double> &angles_gauge,
                        vector<abelian> &angles_conf, link1 link) {

  abelian A(0, 0);
  int place_gauge;
  for (int mu = 0; mu < 4; mu++) {
    link.move(mu, 1);
    place_gauge = link.place / 4;
    link.move(mu, -1);
    A = A +
        angles_conf[link.place + mu] * abelian(1, -angles_gauge[place_gauge]);

    link.move(mu, -1);
    A = A + angles_conf[link.place + mu].conj() *
                abelian(1, -angles_gauge[link.place / 4]);

    link.move(mu, 1);
  }

  return A;
}

void heat_bath_local(
    double *angles_gauge, abelian neighbour, double temperature,
    std::subtract_with_carry_engine<unsigned, 24, 10, 24> &random_generator) {
  double temp_factor;
  long double temp_exponent;

  double neighbour_norm = neighbour.r;
  temp_factor = neighbour.r / temperature;
  double tmp;
  if (temp_factor < 10) {
    temp_exponent = exp(-2 * temp_factor);
    tmp = (double)random_generator() / random_generator.max() *
              (1 - temp_exponent) +
          temp_exponent;
    *angles_gauge = acos(log(tmp) / temp_factor + 1) - neighbour.phi;

  } else {
    do {
      do {
        tmp = (double)random_generator() / random_generator.max();
      } while (tmp == 0);

      tmp = log(tmp) / temp_factor + 1;
    } while (tmp < -1 && tmp > 1);

    *angles_gauge = acos(tmp) - neighbour.phi;
  }
}

void heat_bath_step(vector<double> &angles_gauge, vector<abelian> &angles_conf,
                    double temperature) {

  link1 link(x_size, y_size, z_size, t_size);
  abelian neighbour;

  unsigned seed = time(NULL);
  std::subtract_with_carry_engine<unsigned, 24, 10, 24> random_generator(seed);

  SPACE_ITER_START

  neighbour = angle_neighbour(angles_gauge, angles_conf, link);

  // if (x == 0 && y == 0 && z == 0) {
  //   cout << "neighbour " << neighbour << endl;
  // }

  heat_bath_local(&angles_gauge[link.place / 4], neighbour, temperature,
                  random_generator);

  SPACE_ITER_END
}

void relaxation_step(vector<double> &angles_gauge,
                     vector<abelian> &angles_conf) {

  link1 link(x_size, y_size, z_size, t_size);
  abelian neighbour;

  SPACE_ITER_START

  neighbour = angle_neighbour(angles_gauge, angles_conf, link);

  angles_gauge[link.place / 4] = -neighbour.phi;

  SPACE_ITER_END
}

void overrelaxation_step(vector<double> &angles_gauge,
                         vector<abelian> &angles_conf) {
  link1 link(x_size, y_size, z_size, t_size);
  abelian neighbour;

  SPACE_ITER_START

  neighbour = angle_neighbour(angles_gauge, angles_conf, link);

  angles_gauge[link.place / 4] =
      -angles_gauge[link.place / 4] - 2 * neighbour.phi;

  SPACE_ITER_END
}

vector<double> calculate_binomial_coefficients(int w, int order) {
  vector<double> coefficients(order);

  int nominator = w;
  int denominator = 1;
  coefficients[0] = w;

  for (int i = 1; i < order; i++) {
    nominator = nominator * (w - i);
    denominator = denominator * (i + 1);
    coefficients[i] = (double)nominator / denominator;
  }

  return coefficients;
}

void microcanonical_step(vector<double> &angles_gauge,
                         vector<abelian> &angles_conf, double w, int order) {
  link1 link(x_size, y_size, z_size, t_size);
  abelian neighbour;

  vector<double> binomial_coefficients =
      calculate_binomial_coefficients(w, order);

  abelian gauge_tmp;
  abelian angle_new;

  SPACE_ITER_START

  gauge_tmp = abelian(1, angles_gauge[link.place / 4]) - abelian(1, 0);
  angle_new = abelian(1, 0) + binomial_coefficients[0] * gauge_tmp;

  for (int i = 1; i < order; i++) {
    gauge_tmp = gauge_tmp * gauge_tmp;
    angle_new = angle_new + binomial_coefficients[i] * gauge_tmp;
  }

  angles_gauge[link.place / 4] = angle_new.phi;

  SPACE_ITER_END
}

void thermalize(vector<double> &angles_gauge, vector<abelian> &angles_conf,
                double temperature, int therm_steps) {
  for (int i = 0; i < therm_steps; i++) {
    heat_bath_step(angles_gauge, angles_conf, temperature);

    cout << "thermalization functional: "
         << functional_Landau(angles_conf, angles_gauge) << endl;
  }
}

void simulated_annealing(vector<double> &angles_gauge,
                         vector<abelian> &angles_conf, double temperature_start,
                         double temperature_end, double temperature_step) {
  double temperature = temperature_start;
  while (temperature >= temperature_end) {

    cout << "temperature: " << temperature
         << " functional: " << functional_Landau(angles_conf, angles_gauge)
         << endl;
    heat_bath_step(angles_gauge, angles_conf, temperature);

    temperature -= temperature_step;
  }
}

vector<SA_data> simulated_annealing_test(vector<double> &angles_gauge,
                                         vector<abelian> &angles_conf,
                                         double temperature_start,
                                         double temperature_end,
                                         double temperature_step) {

  vector<SA_data> data;

  SA_data tmp;

  double temperature = temperature_start;
  while (temperature >= temperature_end) {

    heat_bath_step(angles_gauge, angles_conf, temperature);

    tmp.temperature = temperature;

    tmp.functional = functional_Landau(angles_conf, angles_gauge);

    cout << "temperature: " << temperature << " functional: " << tmp.functional
         << endl;

    data.push_back(tmp);

    temperature -= temperature_step;
  }

  return data;
}
