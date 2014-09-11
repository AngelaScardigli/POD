template <class T, class T1>
POD <T, T1>::POD(){
}

template <class T, class T1>
POD <T, T1>::POD(int &mysize, int &myfields){
    dsize=mysize;
    nfields=myfields;
}

template <class T, class T1>
T1 POD <T, T1>::job(){
	return *static_cast<T1*>(this);
}


// ================================================================================== //
// CALC MEAN	                                                                      //
// ================================================================================== //

template <class T, class T1>
vector<T> POD<T, T1>::calc_mean(const bool &tf, const int &n){

    vector< T >    mean;
    vector< T >	   snap;

    if (tf == 1){

    	for (int i = 1; i <= n; i++){
	    job().get_fields(i,"snap", snap);
    	    if(i==1){
		mean=snap*0.0; //guarda se c'e` un modo piu` furbo
		
	    }
	    mean=mean+snap;
	}

	mean=mean/(double(n));

    }
    else{
	job().get_fields(0,"mean", mean);
    }
    return (mean); 

}

