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
	    int nfields;
	    vector<T> mean;
	    vector<T> calc_mean(const bool &);
	    void calc_modes(const bool &);
	protected:
	    POD();
	    POD(int &, int &);
	    T1 job();
	    //void calc_coefficients();
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
	    	void get_fields(int, string, vector<T> &);

};

# include "libPOD.tpp"

#endif
