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

// CC_lib
# include "Operators.hpp"

using namespace std;

/*template <class T>
class POD {
	public:
	    int dsize;
	    int nfields;
	    vector<T> input;
	    POD();
	    POD(int &, int &);
	    vector<T> calc_mean(const bool &, const int &);
	    //virtual void get_fields(int, string, vector<T> &);
};*/

template <class T, class T1>
class POD {
	public:
	    int dsize;
	    int npod;
	    int frows;
	    int ftype;
	    vector<T> mean;
	    vector<vector<double>> lambda;
	    vector<vector<vector<double>>> coeff;
	    vector<T> calc_mean(const bool &);
	    void calc_modes(const bool &);
	    void calc_projection_matrix(const bool &);
	protected:
	    POD();
	    //POD(int &, int &);
	    T1 job();
	    void calc_coefficients();
	    void calc_correlation(vector<vector<vector<double>>> &);
	    void inner_product(vector<T> &, vector<T> &, vector<double> &);

};

template <class T>
class interface : public POD< T, interface<T> > {
	friend class POD<interface, T>;
	public:
		int nsnap;
		int nmode;
		int dimension;
		string flag;
		vector<T> input;
		interface();
		interface(vector<T> &);
		interface(vector<T> &, int &, int &);
		interface(vector<T> &, int &, int &, string);		
	    	void get_fields(int &, string, vector<T> &);
	    	void get_input(vector<T> &, string);
	    	void get_input(vector<vector<T>> &, string);
	    	void give_fields(int &, string, vector<T> &);
	    	void give_output(vector<T> &, string);
	    	void give_output(vector<vector<T>> &, string);


};

# include "libPOD.tpp"

#endif
