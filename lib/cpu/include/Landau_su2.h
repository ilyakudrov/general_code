#include "data.h"
#include "matrix.h"

#include <map>
#include <vector>

void apply_gauge_Landau(
    std::vector<su2> &gauge,
    Data::LatticeData<DataPatternLexicographical, su2> &conf);
double Landau_functional_conf(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf,
    const std::vector<su2> &gauge);
double Landau_su2_functional(
    Data::LatticeData<DataPatternLexicographical, su2> &conf_su2);
std::vector<su2>
generate_gauge_su2_uniform(DataPatternLexicographical &data_pattern);
std::map<double, double> simulated_annealing_thermalization_test(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf,
    std::vector<su2> &gauge, double T_init, double T_final, double T_step,
    int OR_steps, int thermalization_steps, int local_thermalization_steps);
void make_simulated_annealing(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf,
    std::vector<su2> &gauge, double T_init, double T_final, double T_step,
    int OR_steps, int thermalization_steps);
void make_maximization_final(
    const Data::LatticeData<DataPatternLexicographical, su2> &conf,
    std::vector<su2> &gauge, int OR_steps, double tolerance_maximal,
    double tolerance_average);