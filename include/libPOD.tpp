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
vector<T> POD<T, T1>::calc_mean(const bool &tf){

    vector< T >    cmean;
    vector< T >	   snap;

    //nsnap=n;

    if (tf == 1){

    	for (int i = 0; i < job().nsnap; i++){
	    job().get_fields(i,"snap", snap);
    	    if(i==0){
		cmean=snap*0.0; //guarda se c'e` un modo piu` furbo
		
	    }
	    cmean=cmean+snap;
	}

	cmean=cmean/(double(job().nsnap));

    }
    else{
	job().get_fields(0,"mean", cmean);
    }
    mean=cmean;
    return (mean); 

}

// ================================================================================== //
// CALC MODES	                                                                      //
// ================================================================================== //

template <class T, class T1>
void POD<T, T1>::calc_modes(const bool &tf){
   
    vector<T> filter;
    vector<T> volume;
    vector<T> s1;
    vector<T> s2;
    vector<double> corr;
    vector<vector<vector<double>>> Mcorr;   
  
    int n;
    if (job().nmode<job().nsnap){
	n=job().nmode;
    }
    else {
	n=job().nsnap-1;
    }

    if (tf == 1){

	job().get_fields(0,"filter", filter);
	job().get_fields(0,"volume", volume);

	job().get_fields(0,"snap", s1);
	job().get_fields(1,"snap", s2);

	calc_correlation(Mcorr);
	//solve eigenproblem
	//calc_coefficients();
    }
    else {
	//read eigenproblem
    }

}
// calc correlation ================================================================= //

template <class T, class T1>
void POD <T, T1>::calc_correlation(vector<vector<vector<double>>> &Mcorr){
	vector<T> snap1;
	vector<T> snap2;
	vector<double> corr;

	Mcorr.resize(job().nsnap,vector<vector<double>>(job().nsnap,vector<double>(corr.size(),0.0)));
	for (int i=0; i<job().nsnap; ++i){

	    job().get_fields(i,"snap",snap1);
	    snap1=snap1-mean;
	
	    for (int j=i; j<job().nsnap; ++j){

		job().get_fields(j,"snap",snap2);
		snap2=snap2-mean;

		inner_product(snap1, snap2, corr);

		Mcorr[i][j]=corr;
		Mcorr[j][i]=corr;

	    }

	}



}

// inner product ==================================================================== //

template <class T, class T1>
void POD <T, T1>::inner_product(vector<T> &s1, vector<T> &s2, vector<double> &corr){
	vector<T> filter;
	vector<T> volume;
	vector<T> factor;
	double corr1;

	job().get_fields(0,"filter",filter);
	job().get_fields(0,"volume",volume);
	if (job().flag.compare("vectorial")==0){
            corr.resize(s1.size()/3,0.0);
	
	    for (int i=0; i<corr.size(); ++i)
	    {
		sum((Dot_Product(s1[0+3*i],s2[0+3*i])+Dot_Product(s1[1+3*i],s2[1+3*i])+Dot_Product(s1[2+3*i],s2[2+3*i]))*volume,corr1);
		corr[i]=corr1;
	    }
	}
	else {
	    vector<T> vcorr=s1*s2; 
	    vector<T> vfactor=volume*filter;
	    corr.resize(s1.size(),0.0);

	    for (int i=0; i<corr.size(); ++i){

		vcorr[i]=vcorr[i]*vfactor[0];
		sum(vcorr[i],corr1);
		corr[i]=corr1;
	    }

	}

	//return(corr);

}



