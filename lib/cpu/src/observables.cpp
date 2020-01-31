#define data_size 4*x_size*y_size*z_size*t_size
#define PLACE1 (t) * 3 * x_size*y_size*z_size \
				+ (z) * 3 * x_size*y_size \
				+ (y) * 3 * x_size \
				+ (x) * 3 + dir - 1
#define PLACE1_NODIR (t) * 3 * x_size*y_size*z_size \
				+ (z) * 3 * x_size*y_size \
				+ (y) * 3 * x_size \
				+ (x) * 3 - 1
#define PLACE3_LINK_NODIR (link.coordinate[3]) * 3 * x_size*y_size*z_size \
				+ (link.coordinate[2]) * 3 * x_size*y_size \
				+ (link.coordinate[1]) * 3 * x_size \
				+ (link.coordinate[0]) * 3
#define PLACE_FIELD (t) * x_size*y_size*z_size \
				+ (z) * x_size*y_size \
				+ (y) * x_size \
				+ (x)
#define PLACE3 (link.coordinate[3]) * 3 * x_size*y_size*z_size \
				+ (link.coordinate[2]) * 3 * x_size*y_size \
				+ (link.coordinate[1]) * 3 * x_size \
				+ (link.coordinate[0]) * 3 + link.direction - 1
#define PLACE3_NODIR (link.coordinate[3]) * 3 * x_size*y_size*z_size \
				+ (link.coordinate[2]) * 3 * x_size*y_size \
				+ (link.coordinate[1]) * 3 * x_size \
				+ (link.coordinate[0]) * 3 - 1
#define PLACE_FIELD1 (link.coordinate[3]) * x_size*y_size*z_size \
				+ (link.coordinate[2]) * x_size*y_size \
				+ (link.coordinate[1]) * x_size \
				+ (link.coordinate[0])
#define SPACE_ITER_START for (int t = 0; t < t_size; t++) { \
            				for (int z = 0; z < z_size; z++) { \
                				for (int y = 0; y < y_size; y++) { \
                    				for (int x = 0; x < x_size; x++) {
#define SPACE_ITER_END }}}}

#include "observables.h"
#include "link.h"

double plaket_time(const data& conf){
	link1 link(x_size, y_size, z_size, t_size);
	result res(0);
  	double aver[2];
    for (int dir = 1; dir < 4; dir++) {
        link.move_dir(dir);
		SPACE_ITER_START;
        link.go(x, y, z, t);
        res.array.push_back(link.plaket(conf).tr());
		SPACE_ITER_END;
    }
    res.average(aver);
    return aver[0];
}

double plaket_space(const data& conf){
	link1 link(x_size, y_size, z_size, t_size);
	result res(0);
  	double aver[2];
    SPACE_ITER_START;
    link.go(x, y, z, t);
    for(int mu = 1;mu < 4;mu++){
        for(int nu = mu + 1;nu < 4;nu++){
            link.move_dir(nu);
            res.array.push_back(link.plaket_mu(conf, mu).tr());
        }
    }
    SPACE_ITER_END;
    res.average(aver);
    return aver[0];
}

double plaket(const data& conf){
	link1 link(x_size, y_size, z_size, t_size);
	result res(0);
  	double aver[2];
    SPACE_ITER_START;
	link.go(x, y, z, t);
    for(int mu = 1;mu < 5;mu++){
    	for(int nu = mu + 1;nu < 5;nu++){
    		link.move_dir(nu);
        	res.array.push_back(link.plaket_mu(conf, mu).tr());
        }
    }
    SPACE_ITER_END;
    res.average(aver);
    return aver[0];
}

double polyakov(const data& conf){
	link1 link(x_size, y_size, z_size, t_size);
	result res(0);
  	double aver[2];
  	link.move_dir(4);
    SPACE_ITER_START;
    link.go(x, y, z, t);
    res.array.push_back(link.polyakov_loop(conf).tr());
    SPACE_ITER_END;
    res.average(aver);
    return aver[0];
}

double wilson(const data& conf, int R, int T) {
	link1 link(x_size, y_size, z_size, t_size);
	result vec(3 * (data_size/4));
	double aver[2];
	for (int dir = 1; dir < 4; dir++) {
		link.move_dir(dir);
		SPACE_ITER_START;
		link.go(x, y, z, t);
		vec.array[PLACE1] = link.wilson_loop(conf, R, T).tr();
	}
	SPACE_ITER_END;
	vec.average(aver);
	return aver[0];
}

double wilson_abelian(const data& conf, int R, int T/*, double wil[2]*/) {
	link1 link(x_size, y_size, z_size, t_size);
	result vec_cos(3 * (data_size/4));
	double aver[2];
	for (int dir = 1; dir < 4; dir++) {
		link.move_dir(dir);
		SPACE_ITER_START;
		link.go(x, y, z, t);
		vec_cos.array[PLACE1] = link.wilson_loop_abelian(conf, R, T);
	}
	SPACE_ITER_END;
	vec_cos.average(aver);
	return aver[0];
}

void polyakov_abelian(const data& conf, double pol[2]) {
	link1 link(x_size, y_size, z_size, t_size);
	result vec_cos(data_size/4);
	result vec_sin(data_size/4);
	double aver[2];
	double a;
	link.move_dir(4);
	SPACE_ITER_START;
	link.go(x, y, z, t);
	a = link.polyakov_loop_abelian(conf);
	vec_cos.array[PLACE_FIELD] = cos(a);
	vec_sin.array[PLACE_FIELD] = sin(a);
	SPACE_ITER_END;
	vec_cos.average(aver);
	pol[0] = aver[0];
	vec_sin.average(aver);
	pol[1] = aver[0];
}

double wilson_dir(const data& conf, int R, int T, int dir) {
	link1 link(x_size, y_size, z_size, t_size);
	result vec(data_size/4);
	double aver[2];
	link.move_dir(dir);
	SPACE_ITER_START;
	link.go(x, y, z, t);
	vec.array[PLACE_FIELD] = link.wilson_loop(conf, R, T).tr()/2;
	SPACE_ITER_END;
	vec.average(aver);
	return aver[0];
}

void fields(const vector<vector<matrix> >& schwinger_line, const vector<matrix>& plaket, const vector<matrix>& polyakov_loop, vector<vector<result> >& field1, vector<vector<result> >& field2, vector<result>& field3, int d, int D, int x_trans) {
	link1 link(x_size, y_size, z_size, t_size);
	for (int dir = 1; dir < 4; dir++) {
		link.move_dir(dir);
		SPACE_ITER_START;
		link.go(x, y, z, t);
		field3[dir - 1].array[PLACE_FIELD] = link.field3(polyakov_loop, D, x_trans);
		for (int direct = 1; direct < 4; direct++) {
			if (dir != direct) {
				field1[dir - 1][direct - 1].array[PLACE_FIELD] = link.field1(schwinger_line, plaket, polyakov_loop, d, D, direct, x_trans);
				field2[dir - 1][direct - 1].array[PLACE_FIELD] = link.field2(plaket, polyakov_loop, d, D, direct, x_trans);
			}
		}
		SPACE_ITER_END;
	}
}

void field1_average(const vector<vector<matrix> >& schwinger_line, const vector<matrix>& plaket, const vector<matrix>& polyakov_loop, vector<vector<result> >& field1, int d, int D, int x_trans) {
	link1 link(x_size, y_size, z_size, t_size);
	double aver[2];
	for (int dir = 1; dir < 4; dir++) {
		link.move_dir(dir);
		SPACE_ITER_START;
		link.go(x, y, z, t);
		for (int direct = 1; direct < 4; direct++) {
			if (dir != direct) {
				field1[dir - 1][direct - 1].array[PLACE_FIELD] = link.field1(schwinger_line, plaket, polyakov_loop, d, D, direct, x_trans);
			}
		}
		SPACE_ITER_END;
	}
}

vector<vector<matrix> > calculate_schwinger_line(const data& conf, int d, int x_trans) {
	vector<vector<matrix> > vec(3, vector<matrix>(data_size));
	link1 link(x_size, y_size, z_size, t_size);
	for (int dir = 1; dir < 4; dir++) {
		link.move_dir(dir);
		SPACE_ITER_START;
		link.go(x, y, z, t);
		for (int nu = 1; nu < 4; nu++) {
			if (nu != dir) vec[nu - 1][PLACE1] = link.schwinger_line(conf, d, nu, x_trans);
		}
		SPACE_ITER_END;
	}
	return vec;
}

vector<matrix> calculate_plaket(const data& conf) {
	vector<matrix> vec(data_size/4*3);
	link1 link(x_size, y_size, z_size, t_size);
	for (int dir = 1; dir < 4; dir++) {
		link.move_dir(dir);
		SPACE_ITER_START;
		link.go(x, y, z, t);
		vec[PLACE1] = link.plaket(conf);
		SPACE_ITER_END;
	}
	return vec;
}

vector<double> calculate_plaket_time_tr(const data& conf) {
	vector<double> vec(data_size/4*3);
	link1 link(x_size, y_size, z_size, t_size);
	SPACE_ITER_START;
	link.go(x, y, z, t);
	for (int dir = 1; dir < 4; dir++) {
		link.move_dir(dir);
		vec[PLACE1] = link.plaket(conf).tr();
	}
	SPACE_ITER_END;
	return vec;
}

vector<double> calculate_plaket_space_tr(const data& conf) {
	vector<double> vec(data_size/4*3);
	link1 link(x_size, y_size, z_size, t_size);
	int place_dir;
	SPACE_ITER_START;
	link.go(x, y, z, t);
	for (int dir = 1; dir < 4; dir++) {
		for(int j = dir + 1;j < 4;j++){
			link.move_dir(dir);
			vec[PLACE1_NODIR + dir + j - 2] = link.plaket_mu(conf, j).tr();
		}
	}
	SPACE_ITER_END;
	return vec;
}

double plaket4_time_optimized(const vector<double>& plaket_tr, link1& link){
	double a = plaket_tr[PLACE3];
	link.move(link.direction, -1);
	a += plaket_tr[PLACE3];
	link.move(4, 1);
	a += plaket_tr[PLACE3];
	link.move(link.direction, 1);
	a += plaket_tr[PLACE3];
	link.move(4, -1);
	return a / 4;
}

double plaket4_space_optimized(const vector<double>& plaket_tr, link1& link, int nu){
	double a = plaket_tr[PLACE3_NODIR + link.direction + nu - 2];
	link.move(link.direction, -1);
	a += plaket_tr[PLACE3_NODIR + link.direction + nu - 2];
	link.move(nu, 1);
	a += plaket_tr[PLACE3_NODIR + link.direction + nu - 2];
	link.move(link.direction, 1);
	a += plaket_tr[PLACE3_NODIR + link.direction + nu - 2];
	link.move(nu, -1);
	return a / 4;
}

vector<matrix> calculate_polyakov_loop(const data& conf) {
	vector<matrix> vec((data_size) / 4);
	link1 link(x_size, y_size, z_size, t_size);
	int dir = 4;
	link.move_dir(4);
	SPACE_ITER_START;
	link.go(x, y, z, t);
	vec[PLACE_FIELD] = link.polyakov_loop(conf);
	SPACE_ITER_END;
	return vec;
}

vector<matrix> calculate_wilson_loop(const data& conf, int R, int T) {
	vector<matrix> vec(data_size/4*3);
	link1 link(x_size, y_size, z_size, t_size);
	for (int dir = 1; dir < 4; dir++) {
		link.move_dir(dir);
		SPACE_ITER_START;
		link.go(x, y, z, t);
		vec[PLACE1] = link.wilson_loop(conf, R, T);
		SPACE_ITER_END;
	}
	return vec;
}

vector<double> calculate_wilson_loop_tr(const data& conf, int R, int T) {
	vector<double> vec(data_size/4*3);
	link1 link(x_size, y_size, z_size, t_size);
	for (int dir = 1; dir < 4; dir++) {
		link.move_dir(dir);
		SPACE_ITER_START;
		link.go(x, y, z, t);
		vec[PLACE1] = link.wilson_loop(conf, R, T).tr();
		SPACE_ITER_END;
	}
	return vec;
}

double polyakov_loop_corelator(const data& conf, int D){
	vector<matrix> polyakov_loop = calculate_polyakov_loop(conf);
	link1 link(x_size, y_size, z_size, t_size);
	result vec(0);
	double aver[2];
	double a;
	for(int mu = 1;mu < 4;mu++){
		SPACE_ITER_START;
		link.go(x, y, z, t);
		a = polyakov_loop[PLACE_FIELD1].tr();
		link.move(mu, D);
		a *= polyakov_loop[PLACE_FIELD1].conj().tr();
		vec.array.push_back(a);
		SPACE_ITER_END;
	}
	vec.average(aver);
	return aver[0];
}

double plaket_correlator(const vector<matrix>& plaket, int dist) {
	link1 link(x_size, y_size, z_size, t_size);
	matrix A;
	double b;
	result a(data_size);
	double aver[2];
	SPACE_ITER_START;
	link.go(x, y, z, t);
	for (int mu = 1; mu < 4; mu++) {
		link.move_dir(mu);
		A = link.get_matrix(plaket);
		b = A.tr();
		link.move(4, dist);
		A = link.get_matrix(plaket);
		b = b * A.tr();
		link.move(-4, dist);
		a.array.push_back(b);
	}
	SPACE_ITER_END;
	a.average(aver);
	return aver[0];
}

double plaket_correlator_space(const vector<matrix>& plaket, int dist) {
	link1 link(x_size, y_size, z_size, t_size);
	matrix A;
	double b;
	result a(6 * (data_size));
	double aver[2];
	SPACE_ITER_START;
	link.go(x, y, z, t);
	for (int mu = 1; mu < 4; mu++) {
		for (int nu = 1; nu < 4; nu++) {
			if (mu != nu) {
				link.move_dir(mu);
				A = link.get_matrix(plaket);
				b = A.tr();
				link.move(nu, dist);
				A = link.get_matrix(plaket);
				b = b * A.tr();
				link.move(-nu, dist);
				a.array.push_back(b);
			}
		}
	}
	SPACE_ITER_END;
	a.average(aver);
	return aver[0];
}

result wilson_plaket_correlator_electric_optimized(const data& conf, const vector<double>& wilson_loop_tr, const vector<double>& plaket_tr, int R, int T, int x_trans, int d_min, int d_max){
	link1 link(x_size, y_size, z_size, t_size);
	double vec[d_max - d_min + 1];
	for(int i = 0;i < d_max - d_min + 1;i++){
		vec[i] = 0;
	}
    result final(0);
    double aver[2];
    double a;
    for (int dir = 1; dir < 4; dir++) {
    	SPACE_ITER_START;
        link.go(x, y, z, t);
        link.move_dir(dir);
        a = wilson_loop_tr[PLACE1];
        link.move(4, T/2);
		link.move(dir, d_min);
		for(int d = d_min;d <= d_max;d++){
			if(x_trans == 0){
				for(int mu = 1;mu < 4;mu++){
					link.move_dir(mu);
					vec[d - d_min] += a * plaket4_time_optimized(plaket_tr, link);
				}
			}
			else{
        		for(int nu = 1;nu < 4;nu++){
        			if(nu != dir){
        				link.move(nu, x_trans);
        				for(int mu = 1;mu < 4;mu++){
        					link.move_dir(mu);
        					vec[d - d_min] += a * plaket4_time_optimized(plaket_tr, link);
        				}
        				link.move(nu, -2*x_trans);
        				for(int mu = 1;mu < 4;mu++){
        					link.move_dir(mu);
        					vec[d - d_min] += a * plaket4_time_optimized(plaket_tr, link);
        				}
        				link.move(nu, x_trans);
        			}
        		}
			}
			link.move(dir, 1);
		}
    	SPACE_ITER_END;
    }
	int count;
	if(x_trans == 0) count = data_size / 4 * 9;
	else count = data_size / 4 * 36;
	for(int d = d_min;d <= d_max;d++){
    	final.array.push_back(vec[d - d_min]/count);
	}
    return final;
}

result wilson_plaket_correlator_electric_new(const data& conf, const vector<double>& wilson_loop_tr, int R, int T, int x_trans, int d_min, int d_max){
	link1 link(x_size, y_size, z_size, t_size);
	result vec(0);
    result final(0);
    double aver[2];
    double a;
    for(int d = d_min;d <= d_max;d++){
    	for (int dir = 1; dir < 4; dir++) {
        	SPACE_ITER_START;
            link.go(x, y, z, t);
            link.move_dir(dir);
            a = wilson_loop_tr[PLACE1];
            link.move(4, T/2);
            link.move(dir, d);
            for(int nu = 1;nu < 4;nu++){
            	if(nu != dir){
            		link.move(nu, x_trans);
            		for(int mu = 1;mu < 4;mu++){
            			link.move_dir(mu);
            			if(T%2 == 0) vec.array.push_back(a * link.plaket_implement4(conf, 4).tr());
            			if(T%2 == 1) vec.array.push_back(a * link.plaket_implement2(conf, 4).tr());
            		}
            		if(x_trans != 0){
            		link.move(nu, -2*x_trans);
            		for(int mu = 1;mu < 4;mu++){
            			link.move_dir(mu);
            			if(T%2 == 0) vec.array.push_back(a * link.plaket_implement4(conf, 4).tr());
            			if(T%2 == 1) vec.array.push_back(a * link.plaket_implement2(conf, 4).tr());
            		}
            		link.move(nu, x_trans);
            		}
            		else link.move(nu, -x_trans);
            	}
            }
       		SPACE_ITER_END;
    	}
    	vec.average(aver);
    	final.array.push_back(aver[0]);
    	vec.array.clear();
    }
    return final;
}

result wilson_plaket_correlator_electric_x_new(const data& conf, vector<double> wilson_loop_tr, int R, int T, int x_trans_min, int x_trans_max, int d){
	link1 link(x_size, y_size, z_size, t_size);
	result vec(0);
    result final(0);
    double aver[2];
    double a;
    for(int x_trans = x_trans_min;x_trans <= x_trans_max;x_trans++){
    	for (int dir = 1; dir < 4; dir++) {
        	SPACE_ITER_START;
            link.go(x, y, z, t);
            //link.move_dir(dir);
            a = wilson_loop_tr[PLACE1];
            link.move(4, T/2);
            link.move(dir, d);
            for(int nu = 1;nu < 4;nu++){
            	if(nu != dir){
            		link.move(nu, x_trans);
            		for(int mu = 1;mu < 4;mu++){
            			link.move_dir(mu);
            			if(T%2 == 0) vec.array.push_back(a * link.plaket_implement4(conf, 4).tr());
            			if(T%2 == 1) vec.array.push_back(a * link.plaket_implement2(conf, 4).tr());
            		}
            		link.move(nu, -x_trans);
            	}
            }
       		SPACE_ITER_END;
    	}
    	vec.average(aver);
    	final.array.push_back(aver[0]);
    	vec.array.clear();
    }
    return final;
}

result polyakov_plaket_correlator_electric(const data& conf, const data& smeared, int R, int x_trans, int d_min, int d_max){
	vector<matrix> polyakov_loop = calculate_polyakov_loop(smeared);
	link1 link(x_size, y_size, z_size, t_size);
	result vec(0);
    result final(0);
    double aver[2];
    double a;
    for(int d = d_min;d <= d_max;d++){
    	for (int dir = 1; dir < 4; dir++) {
        	SPACE_ITER_START;
        	link.go(x, y, z, t);
        	/*if(x == 1 && y == 1 && z == 1 && t == 1 && dir == 1){
        		cout<<"polyakov1 ";
        	 	link.print_link();
        	}*/
        	a = polyakov_loop[PLACE_FIELD1].tr();
        	link.move(dir, R);
        	/*if(x == 1 && y == 1 && z == 1 && t == 1 && dir == 1){
        		cout<<"polyakov2 ";
        	 	link.print_link();
        	}*/
        	a *= polyakov_loop[PLACE_FIELD1].conj().tr();
        	link.move(dir, d-R);
        	for(int nu = 1;nu < 4;nu++){
        		if(nu != dir){
        			link.move(nu, x_trans);
        			for(int mu = 1;mu < 4;mu++){
        				link.move_dir(mu);
        				/*if(x == 1 && y == 1 && z == 1 && t == 1 && dir == 1){
        					cout<<"plaket ";
        	 				link.print_link();
        				}*/
        				vec.array.push_back(a * link.plaket_implement4(conf, 4).tr());
        				//vec.array.push_back(a * link.plaket_mu(conf, 4).tr());
        			}
        			if(x_trans != 0){
        			link.move(nu, -2*x_trans);
        			for(int mu = 1;mu < 4;mu++){
        				link.move_dir(mu);
        				vec.array.push_back(a * link.plaket_implement4(conf, 4).tr());
        			}
        			link.move(nu, x_trans);
        			}
        			else link.move(nu, -x_trans);
        		}
        	}
       		SPACE_ITER_END;
    	}
    	vec.average(aver);
    	final.array.push_back(aver[0]);
    	vec.array.clear();
    }
    return final;
}

result wilson_plaket_correlator_magnetic_optimized(const data& conf, const vector<double>& wilson_loop_tr, const vector<double>& plaket_tr, int R, int T, int x_trans, int d_min, int d_max){
	link1 link(x_size, y_size, z_size, t_size);
	double vec[d_max - d_min + 1];
	for(int i = 0;i < d_max - d_min + 1;i++){
		vec[i] = 0;
	}
    result final(0);
    double aver[2];
    double a;
    for (int dir = 1; dir < 4; dir++) {
    	SPACE_ITER_START;
        link.go(x, y, z, t);
        link.move_dir(dir);
        a = wilson_loop_tr[PLACE1];
        link.move(4, T/2);
		link.move(dir, d_min);
		for(int d = d_min;d <= d_max;d++){
			if(x_trans == 0){
				for(int mu = 1;mu < 4;mu++){
            		for(int j = mu + 1;j < 4;j++){
            			link.move_dir(mu);
            			vec[d - d_min] += a * plaket4_space_optimized(plaket_tr, link, j);
            		}
            	}
			}
			else{
        		for(int nu = 1;nu < 4;nu++){
        			if(nu != dir){
        				link.move(nu, x_trans);
        				for(int mu = 1;mu < 4;mu++){
            				for(int j = mu + 1;j < 4;j++){
            					link.move_dir(mu);
            					vec[d - d_min] += a * plaket4_space_optimized(plaket_tr, link, j);
            				}
            			}
        				link.move(nu, -2*x_trans);
        				for(int mu = 1;mu < 4;mu++){
            				for(int j = mu + 1;j < 4;j++){
            					link.move_dir(mu);
            					vec[d - d_min] += a * plaket4_space_optimized(plaket_tr, link, j) ;
            				}
            			}
        				link.move(nu, x_trans);
        			}
        		}
			}
			link.move(dir, 1);
		}
    	SPACE_ITER_END;
    }
	int count;
	if(x_trans == 0) count = data_size / 4 * 9;
	else count = data_size / 4 * 36;
	for(int d = d_min;d <= d_max;d++){
    	final.array.push_back(vec[d - d_min]/count);
	}
    return final;
}

result wilson_plaket_correlator_magnetic_new(const data& conf, vector<double> wilson_loop_tr, int R, int T, int x_trans, int d_min, int d_max){
	link1 link(x_size, y_size, z_size, t_size);
	result vec(0);
    result final(0);
    double aver[2];
    double a;
    for(int d = d_min;d <= d_max;d++){
    	for (int dir = 1; dir < 4; dir++) {
        	SPACE_ITER_START;
            link.go(x, y, z, t);
            link.move_dir(dir);
            a = wilson_loop_tr[PLACE1];
            link.move(4, T/2);
            link.move(dir, d);
            for(int nu = 1;nu < 4;nu++){
            	if(nu != dir){
            		link.move(nu, x_trans);
            		for(int mu = 1;mu < 4;mu++){
            			for(int j = mu + 1;j < 4;j++){
            				link.move_dir(mu);
            				vec.array.push_back(a * link.plaket_implement4(conf, j).tr());
            			}
            		}
            		if(x_trans != 0){
            		link.move(nu, -2*x_trans);
            		for(int mu = 1;mu < 4;mu++){
            			for(int j = mu + 1;j < 4;j++){
            				link.move_dir(mu);
            				vec.array.push_back(a * link.plaket_implement4(conf, j).tr());
            			}
            		}
            		link.move(nu, x_trans);
            		}
            		else link.move(nu, x_trans);
            	}
            }
            SPACE_ITER_END;
       	}
       	vec.average(aver);
    	final.array.push_back(aver[0]);
    	vec.array.clear();
    }
    return final;
}

result polyakov_plaket_correlator_magnetic(const data& conf, const data& smeared, int R, int x_trans, int d_min, int d_max){
	vector<matrix> polyakov_loop = calculate_polyakov_loop(smeared);
	link1 link(x_size, y_size, z_size, t_size);
	result vec(0);
    result final(0);
    double aver[2];
    double a;
    for(int d = d_min;d <= d_max;d++){
    	for (int dir = 1; dir < 4; dir++) {
    	//for (int dir = 1; dir < 2; dir++) {
        	SPACE_ITER_START;
            link.go(x, y, z, t);
            a = polyakov_loop[PLACE_FIELD1].tr();
            link.move(dir, R);
            a *= polyakov_loop[PLACE_FIELD1].conj().tr();
            link.move(dir, d-R);
            for(int nu = 1;nu < 4;nu++){
            	if(nu != dir){
            		link.move(nu, x_trans);
            		for(int mu = 1;mu < 4;mu++){
            			if(dir != mu){
            				for(int j = 1;j < 4;j++){
            					if(j != dir && mu < j){
            						link.move_dir(mu);
            						vec.array.push_back(a * link.plaket_implement4(conf, j).tr());
            					}
            				}
            			}
            		}
            		link.move(nu, -2*x_trans);
            		for(int mu = 1;mu < 4;mu++){
            			if(dir != mu){
            				for(int j = 1;j < 4;j++){
            					if(j != dir && mu < j){
            						link.move_dir(mu);
            						vec.array.push_back(a * link.plaket_implement4(conf, j).tr());
            					}
            				}
            			}
            		}
            		link.move(nu, x_trans);
            	}
            }
            SPACE_ITER_END;
       	}
       	vec.average(aver);
    	final.array.push_back(aver[0]);
    	vec.array.clear();
    }
    return final;
}

result wilson_plaket_correlator_magnetic_x_new(const data& conf, vector<double> wilson_loop_tr, int R, int T, int x_trans_min, int x_trans_max, int d){
	link1 link(x_size, y_size, z_size, t_size);
	result vec(0);
    result final(0);
    double aver[2];
    double a;
    for(int x_trans = x_trans_min;x_trans <= x_trans_max;x_trans++){
    	for (int dir = 1; dir < 4; dir++) {
        	SPACE_ITER_START;
            link.go(x, y, z, t);
            link.move_dir(dir);
            a = wilson_loop_tr[PLACE1];
            link.move(4, T/2);
            link.move(dir, d);
            for(int nu = 1;nu < 4;nu++){
            	if(nu != dir){
            		link.move(nu, x_trans);
            		for(int mu = 1;mu < 4;mu++){
            			for(int i = 1;i < 4;i++){
            				if(i != mu){
            					for(int j = 1;j < 4;j++){
            						if(i != j){
            							link.move_dir(i);
            							vec.array.push_back(a * link.plaket_implement4(conf, j).tr());
            						}
            					}
            				}
            			}
            		}
            		link.move(nu, -x_trans);
            	}
            }
            SPACE_ITER_END;
       	}
       	vec.average(aver);
    	final.array.push_back(aver[0]);
    	vec.array.clear();
    }
    return final;
}

void push_result(result& values1_out, result& values2_out, result& values3_out, vector<vector<result> >& field1_values, vector<vector<result> >& field2_values, vector<result>& field3_values){
	double aver1[3][3][2];
        double aver2[3][3][2];
        double aver3[3][2];
	for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                        if (i != j) {
                                field1_values[i][j].average(aver1[i][j]);
                                field2_values[i][j].average(aver2[i][j]);
                        }
                }
                field3_values[i].average(aver3[i]);
        }
        for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                        if (i != j) {
                                values1_out.array.push_back(aver1[i][j][0]);
                                values2_out.array.push_back(aver2[i][j][0]);
                        }
                }
                values3_out.array.push_back(aver3[i][0]);
        }
}

void calculate_and_push(data& conf, result& values1_out, result& values2_out, result& values3_out, vector<vector<result> >& field1_values, vector<vector<result> >& field2_values, vector<result>& field3_values, vector<vector<matrix> >& schwinger_line, vector<matrix>& plaket, vector<matrix>& polyakov_loop, int d, int D, int x_trans){
        schwinger_line = calculate_schwinger_line(conf, d, x_trans);
        plaket = calculate_plaket(conf);
        polyakov_loop = calculate_polyakov_loop(conf);

        fields(schwinger_line, plaket, polyakov_loop, field1_values, field2_values, field3_values, d, D, x_trans);
        push_result(values1_out, values2_out, values3_out, field1_values, field2_values, field3_values);
}

//monopoles

void length(loop* ll, int& ss1){
    if(ll->link[0]==NULL) return ;
    int i=0;
    int Dir=0;
    do {
       length(ll->link[i], ss1);
       Dir=ll->get_dir(i+1);

       if( Dir!=0 ) ss1=ss1+1;

       i++;
    } while((ll->link[i]!=NULL)&&(i<=6));
}

result calculate_cluster_lengths(vector<loop*>& LL, int& max_number){
	int n = 0;
	result res(0);
	int count = 0;
	for(int i = 0;i < LL.size();i++){
		n = 0;
		length(LL[i], n);
		res.array.push_back(n);
		if(n > count) {
			count = n;
			max_number = i;
		}
	}
	return res;
}

void length_mu(loop* ll, int mu, int& s){
	if(ll->link[0]==NULL) return ;
	int i=0;
	int dir=0;
	do {
		length_mu(ll->link[i], mu, s);
		dir=ll->get_dir(i+1);
		if(dir == mu) s+=1;
		if(dir == -mu) s-=1;
		i++;
	} while((ll->link[i]!=NULL)&&(i<=6));
}

void calculate_t_clusters(vector<loop*>& LL, vector<loop*>& t_clusters, int max_number){
	int s = 0;
	for(int i = 0;i < LL.size();i++){
		if(i != max_number){
			s = 0;
			length_mu(LL[i], 4, s);
			if(s != 0) t_clusters.push_back(LL[i]);
		}
	}
}

void calculate_t_clusters_n(vector<loop*>& LL, vector<loop*>& t_clusters_n, int max_number, int n){
	int s = 0;
	for(int i = 0;i < LL.size();i++){
		if(i != max_number){
			s = 0;
			length_mu(LL[i], 4, s);
			if(abs(s/t_size) == n) t_clusters_n.push_back(LL[i]);
		}
	}
}

void calculate_s_clusters(vector<loop*>& LL, vector<loop*>& s_clusters, int max_number){
	int s = 0;
	for(int i = 0;i < LL.size();i++){
		if(i != max_number){
			for(int j = 1;j < 4;j++){
				s = 0;
				length_mu(LL[i], j, s);
				if(s != 0) s_clusters.push_back(LL[i]);
			}
		}
	}
}

double t_density_n(vector<loop*>& t_clusters, int n){
	int s = 0;
	int count = 0;
	for(int i = 0;i < t_clusters.size();i++){
		s = 0;
		length_mu(t_clusters[i], 4, s);
		if(abs(s/t_size) == n) count++;
	}
	return (double)count;
}

double time_length_portion(vector<loop*>& t_clusters){
	result res(0);
	int s1 = 0;
	int s2 = 0;
	for(int i = 0;i < t_clusters.size();i++){
		s1 = 0;
		s2 = 0;
		length(t_clusters[i], s1);
		length_mu(t_clusters[i], 4, s2);
		res.array.push_back(fabs(1.*s1/s2));
	}
	double aver[2];
	res.average(aver);
	return aver[0];
}

void sites_unique(loop* ll, vector<loop*>& sites){
	int a = 0;
	for(int r=0;r < sites.size();r++){
		if(sites[r]->node.coordinate[0] == ll->node.coordinate[0]
			&& sites[r]->node.coordinate[1] == ll->node.coordinate[1]
			&& sites[r]->node.coordinate[2] == ll->node.coordinate[2]
			&& sites[r]->node.coordinate[3] == ll->node.coordinate[3]) a = 1;
	}
	if(a != 1) sites.push_back(ll);
	int i = 0;
	while (ll->link[i]!=NULL && i<=6){
		sites_unique(ll->link[i], sites);
		i++;
	}
}

void aver_r(vector<loop*> sites, double* aver_coord){
	int size = sites.size();
	aver_coord[0] = 0; aver_coord[1] = 0; aver_coord[2] = 0;
	for(int k = 0;k < size;k++){
		aver_coord[0] += 1.*sites[k]->node.coordinate[0]/size;
		aver_coord[1] += 1.*sites[k]->node.coordinate[1]/size;
		aver_coord[2] += 1.*sites[k]->node.coordinate[2]/size;
	}
}

double distance_shortest(double a, double b){
	if(fabs(a - b) <= (t_size - fabs(a - b))) return fabs(a - b);
	else return (t_size - fabs(a - b));
}

double disp_r(vector<loop*>& sites, double* aver_coord){
	double disp = 0;
	double dist_x = 0; double dist_y = 0; double dist_z = 0;
	for(int k = 0;k < sites.size();k++){
		dist_x = distance_shortest(sites[k]->node.coordinate[0], aver_coord[0]);
		dist_y = distance_shortest(sites[k]->node.coordinate[1], aver_coord[1]);
		dist_z = distance_shortest(sites[k]->node.coordinate[2], aver_coord[2]);
		disp += dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
	}
	disp = disp/sites.size();
	return disp;
}

double calculate_disp_r(vector<loop*>& t_clusters){
	result res(0);
	vector<loop*> sites(0);
	double aver_coord[3];
	for(int i = 0;i < t_clusters.size();i++){
		sites_unique(t_clusters[i], sites);
		aver_r(sites, aver_coord);
		res.array.push_back(1./disp_r(sites, aver_coord));
	}
	sites.clear();
	double aver[2];
	res.average(aver);
	return aver[0];
}

bool sites_close(loop* l, loop* ll){
	int x = distance_shortest(ll->node.coordinate[0], l->node.coordinate[0]);
	int y = distance_shortest(ll->node.coordinate[1], l->node.coordinate[1]);
	int z = distance_shortest(ll->node.coordinate[2], l->node.coordinate[2]);
	int t = distance_shortest(ll->node.coordinate[3], l->node.coordinate[3]);
	return ((x*x + y*y + z*z + t*t) == 1);
}

double dimension(vector<loop*> sites) {
	int count = 0;
	for(int i = 0;i < sites.size();i++){
		for(int j = 0;j < sites.size();j++){
			if(sites_close(sites[i], sites[j])) count++;
		}
	}
	return 1.*count/sites.size();
}

double charge_difference(vector<loop*>& t_clusters_1){
	int count1 = 0;
	int count2 = 0;
	int t_length = 0;
	for(int i = 0;i < t_clusters_1.size();i++){
		t_length = 0;
		length_mu(t_clusters_1[i], 4, t_length);
		if(t_length > 0) count1++;
		if(t_length < 0) count2++;
	}
	return (double)(count1 - count2);
}
