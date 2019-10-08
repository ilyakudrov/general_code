#include "eigen.h"
#include "result.h"
#include "data.h"

using namespace std;

int x_size;
int y_size;
int z_size;
int t_size;

int main(int argc, char* argv[]){
    x_size = atof(argv[1]);
    y_size = atof(argv[2]);
    z_size = atof(argv[3]);
    t_size = atof(argv[4]);
    int mass = 0.0075;
    int mu_q = 0;
    double tolerance = 1e-2;
    char path_vec[] = "../../eigenvectors/time_32/mu0.00/eigenvectors_neig=60_0001";
    char path_val[] = "../../eigenvectors/time_32/mu0.00/eigenvalues_neig=60_0001";
    char path_conf[] = "../../confs/su2/time_32/mu0.00/conf_0001.fl";
    int vec_size = x_size * y_size * z_size * t_size * 2;
    int n_eig = 60;

    data conf;
    conf.read_float(path_conf);
    result raw_vecs(vec_size * 1/*n_eig*/ * 2);
    result vals(n_eig * 2);
    raw_vecs.read(path_vec, vec_size * 1/*n_eig*/ * 2);
    vals.read(path_val, n_eig * 2);
    vector<scomplex_t> vec;
    scomplex_t tmp_vec;
    scomplex_t tmp_val;
    for(int i = 0;i < 1/*n_eig*/;i++){
        for(int j = 0;j < vec_size;j++){
            tmp_vec.re = raw_vecs.array[i * n_eig + 2 * j];
            tmp_vec.im = raw_vecs.array[i * n_eig + 2 * j + 1];
            vec.push_back(tmp_vec);
        }
        tmp_val.re = vals.array[2 * i];
        tmp_val.im = vals.array[2 * i + 1];
        test_eigenvector(&vec[0], tmp_val, vec_size, &conf.array[0], x_size, y_size, z_size, t_size, mass, mu_q, tolerance);

    }
}