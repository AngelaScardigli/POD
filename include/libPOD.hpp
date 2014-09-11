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
	    vector<T> calc_mean(const bool &, const int &);
	protected:
	    POD();
	    POD(int &, int &);
	    T1 job();
	    //void get_fields(int, string, vector<T> &);


};

template <class T>
class interface : public POD< T, interface<T> > {
	friend class POD<interface, T>;
	public:
		vector<T> input;
		interface();
		interface(vector<T> &);
	    	void get_fields(int, string, vector<T> &);

};

# include "libPOD.tpp"

#endif
