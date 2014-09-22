#ifndef __POD_HH__
#define __POD_HH__

// Standard Template Library
# include <cmath>
# include <array>
# include <vector>
# include <string>
# include <sstream>
# include <fstream>
# include <iomanip>
# include <iostream>
# include <algorithm>
# include <lapacke.h>
# include <mpi.h>

// CC_lib
# include "Operators.hpp"

using namespace std;

template <class T, class T1>
class POD {
	public:
	    int dsize;
	    int npod;
	    vector<T> volume;
	    vector<T> filter;
	    vector<T> mean;
	    vector<vector<double>> lambda;
	    vector<vector<vector<double>>> coeff;
	    vector<vector<T>> phi;
	    void calc_mean(const bool &);
	    void calc_modes(const bool &);
	    void calc_projection_matrix(const bool &);
	    void calc_reconstruction(const bool &);
	protected:
	    POD();
	    T1 job();
	    void calc_coefficients();
	    void calc_correlation(vector<vector<vector<double>>> &);
	    vector<double> inner_product(vector<vector<double>> &, vector<vector<double>> &);
	    vector<double> inner_product(vector<vector<vector<double>>> &, vector<vector<vector<double>>> &);
	    vector<vector<vector<double>>> ass_proj_matrix(vector<vector<double>> &);
	    vector<vector<vector<double>>> ass_proj_matrix(vector<vector<vector<double>>> &);
	    vector<vector<T>> ass_phi(vector<vector<double>> &);
	    vector<vector<T>> ass_phi(vector<vector<vector<double>>> &);
	    vector<vector<double>> project(vector<vector<double>> &);
	    vector<vector<double>> project(vector<vector<vector<double>>> &);
	    void reconstruct(vector<vector<double>> &, vector<vector<double>> &);
	    void reconstruct(vector<vector<double>> &, vector<vector<vector<double>>> &);


};

template <class T>
class interface : public POD< T, interface<T> > {
	friend class POD<interface, T>;
	public:
		int nsnap;
		int nmode;
		int nrecon;
		vector<T> input;
		interface(vector<T> &, int &, int &);
		interface(vector<T> &, int &, int &, int &);		
	    	void get_fields(int &, string, vector<T> &);
	    	void get_input(vector<double> &, string);
	    	void get_input(vector<vector<double>> &, string);
	    	void get_input(vector<vector<vector<double>>> &, string);
	    	void give_fields(int &, string, vector<T> &);
	    	void give_output(vector<double> &, string);
	    	void give_output(vector<vector<double>> &, string);
	    	void give_output(vector<vector<vector<double>>> &, string);


};

# include "libPOD.tpp"

#endif
