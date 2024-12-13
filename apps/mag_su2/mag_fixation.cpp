#include "../../lib/cpu/include/data.h"
#include "../../lib/cpu/include/mag.h"
#include "../../lib/cpu/include/matrix.h"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

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
  string path_previous;
  string path_conf_output;
  string path_spins_output;
  string path_functional_output;

  double T_step;
  double T_init;
  double T_final;
  int OR_steps;
  int thermalization_steps;
  double tolerance_maximal;
  double tolerance_average;
  int tolerance_digits;
  int gauge_copies;

  int ml5_conf_num = 0;

  bool is_new_trial = true;
  bool is_final = true;
  bool is_compare = false;
  bool is_compare_spins = false;
  bool is_functional_save = false;

  // read parameters
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]) == "-conf_format") {
      conf_format = argv[++i];
    } else if (string(argv[i]) == "-path_conf") {
      path_conf = argv[++i];
    } else if (string(argv[i]) == "-path_previous") {
      path_previous = argv[++i];
    } else if (string(argv[i]) == "-path_conf_output") {
      path_conf_output = argv[++i];
    } else if (string(argv[i]) == "-path_spins_output") {
      path_spins_output = argv[++i];
    } else if (string(argv[i]) == "-path_functional_output") {
      path_functional_output = argv[++i];
    } else if (string(argv[i]) == "-bytes_skip") {
      bytes_skip = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-T_step") {
      T_step = stod(string(argv[++i]));
    } else if (string(argv[i]) == "-T_init") {
      T_init = stod(string(argv[++i]));
    } else if (string(argv[i]) == "-T_final") {
      T_final = stod(string(argv[++i]));
    } else if (string(argv[i]) == "-OR_steps") {
      OR_steps = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-tolerance_maximal") {
      tolerance_maximal = stod(string(argv[++i]));
    } else if (string(argv[i]) == "-tolerance_average") {
      tolerance_average = stod(string(argv[++i]));
    } else if (string(argv[i]) == "-thermalization_steps") {
      thermalization_steps = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-tolerance_digits") {
      tolerance_digits = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-gauge_copies") {
      gauge_copies = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-x_size") {
      x_size = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-y_size") {
      y_size = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-z_size") {
      z_size = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-t_size") {
      t_size = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-ml5_conf_num") {
      ml5_conf_num = stoi(string(argv[++i]));
    } else if (string(argv[i]) == "-is_new_trial") {
      istringstream(string(argv[++i])) >> is_new_trial;
    } else if (string(argv[i]) == "-is_final") {
      istringstream(string(argv[++i])) >> is_final;
    } else if (string(argv[i]) == "-is_compare") {
      istringstream(string(argv[++i])) >> is_compare;
    } else if (string(argv[i]) == "-is_compare_spins") {
      istringstream(string(argv[++i])) >> is_compare_spins;
    } else if (string(argv[i]) == "-is_functional_save") {
      istringstream(string(argv[++i])) >> is_functional_save;
    } else
      cout << "unknown parameter " << argv[i] << endl;
  }

  cout << "path_conf " << path_conf << endl;
  cout << "conf_format " << conf_format << endl;
  cout << "path_conf_output " << path_conf_output << endl;
  cout << "path_spins_output " << path_spins_output << endl;
  cout << "path_previous " << path_previous << endl;
  cout << "path_functional_output " << path_functional_output << endl;
  cout << "bytes_skip " << bytes_skip << endl;
  cout << "OR_steps " << OR_steps << endl;
  cout << "is_new_trial " << is_new_trial << endl;
  if (is_new_trial) {
    cout << "T_step " << T_step << endl;
    cout << "T_init " << T_init << endl;
    cout << "T_final " << T_final << endl;
    cout << "thermalization_steps " << thermalization_steps << endl;
  }
  cout << "tolerance_maximal " << tolerance_maximal << endl;
  cout << "tolerance_average " << tolerance_average << endl;
  cout << "tolerance_digits " << tolerance_digits << endl;
  cout << "gauge_copies " << gauge_copies << endl;
  cout << "x_size " << x_size << endl;
  cout << "y_size " << y_size << endl;
  cout << "z_size " << z_size << endl;
  cout << "t_size " << t_size << endl;
  cout << "ml5_conf_num " << ml5_conf_num << endl;

  cout << "is_new_trial " << is_new_trial << endl;
  cout << "is_final " << is_final << endl;
  cout << "is_compare " << is_compare << endl;
  cout << "is_compare_spins " << is_compare_spins << endl;
  cout << "is_functional_save " << is_functional_save << endl;

  data<su2> conf_su2;

  // for ml5 configuration
  vector<float> ml5_data;

  // read configuration
  get_data(conf_su2, path_conf, conf_format, bytes_skip, 0);

  cout.precision(17);

  vector<spin> spins;
  map<int, double> functional;

  // cycle over gauge copies
  for (int copy = 1; copy <= gauge_copies; copy++) {

    // make new trial
    if (is_new_trial) {

      cout << endl;
      cout << "making new trial" << endl << endl;

      // generate uniformly distributed spins
      spins = generate_spins_uniform();

      cout << endl << "starting simulated annealing" << endl << endl;

      start_time = clock();

      make_simulated_annealing(conf_su2.array, spins, T_init, T_final, T_step,
                               OR_steps, thermalization_steps);

      end_time = clock();
      search_time = end_time - start_time;
      cout << "simulated annealing time: " << search_time * 1. / CLOCKS_PER_SEC
           << endl;
    } else {
      spins = read_spins(path_previous);
    }

    cout << endl << "starting RO process " << endl;

    if (is_final) {

      start_time = clock();

      make_maximization_final(conf_su2.array, spins, OR_steps,
                              tolerance_maximal, tolerance_average);

      end_time = clock();
      search_time = end_time - start_time;
      cout << "maximization final time: " << search_time * 1. / CLOCKS_PER_SEC
           << endl;

      gauge_tranformation_spins(conf_su2.array, spins);

      cout << "final functional is " << MAG_functional_su2(conf_su2.array)
           << endl;

      conf_su2.write_double(path_conf_output);

      break;
    } else {

      start_time = clock();

      make_maximization_approximate(conf_su2.array, spins, OR_steps,
                                    tolerance_digits);

      end_time = clock();
      search_time = end_time - start_time;
      cout << "maximization approximate time: "
           << search_time * 1. / CLOCKS_PER_SEC << endl;

      if (is_functional_save) {
        functional[copy] = MAG_functional_su2_spin(conf_su2.array, spins);
      }

      if (is_compare) {

        double functional_old;
        double functional_new;

        if (is_compare_spins) {
          vector<spin> spins_old = read_spins(path_previous);

          functional_old = MAG_functional_su2_spin(conf_su2.array, spins_old);
        } else {
          data<su2> conf_old;
          conf_old.read_double(path_previous, 0);

          functional_old = MAG_functional_su2(conf_old.array);
        }

        functional_new = MAG_functional_su2_spin(conf_su2.array, spins);

        cout << "new functional is " << functional_new << endl;
        cout << "old functional is " << functional_old << endl;

        if (functional_old < functional_new) {
          cout << "new functional is higher, saving spin configuration" << endl;
          write_spins(path_spins_output, spins);
        } else {
          cout << "new functional is lower, spin configuration is not saved"
               << endl;
        }

      } else
        write_spins(path_spins_output, spins);
    }
  }

  if (is_functional_save) {
    ofstream functional_stream;
    functional_stream.open(path_functional_output);

    functional_stream.precision(17);

    functional_stream << "copy,functional" << endl;

    for (auto it = functional.begin(); it != functional.end(); it++) {
      functional_stream << it->first << "," << it->second << endl;
    }

    functional_stream.close();
  }
}
