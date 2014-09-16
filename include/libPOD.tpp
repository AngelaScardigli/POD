template <class T, class T1>
POD <T, T1>::POD(){
}

/*template <class T, class T1>
POD <T, T1>::POD(int &mysize, int &myfields){
    dsize=mysize;
    nfields=myfields;
}*/

template <class T, class T1>
T1 POD <T, T1>::job(){
	return *static_cast<T1*>(this);
}


// ================================================================================== //
// CALC MEAN	                                                                      //
// ================================================================================== //

template <class T, class T1>
vector<T> POD<T, T1>::calc_mean(const bool &tf){

    //vector< T >    mean;
    vector< T >	   snap;
    int fi = 0;

    //nsnap=n;

    if (tf == 1){

    	for (int i = 0; i < job().nsnap; i++){
	    job().get_fields(i,"snap", snap);
    	    if(i==0){
		mean=snap*0.0; //guarda se c'e` un modo piu` furbo
		
	    }
	    mean=mean+snap;
	}

	mean=mean/(double(job().nsnap));
	job().give_fields(fi,"mean",mean);

    }
    else{
	job().get_fields(fi,"mean", mean);
    }

    if (job().flag.compare("vectorial")==0){
	ftype=3;
    }
    else {
	ftype=1;
    }

    dsize=mean.size()/ftype;
    frows=mean.size();
    return (mean); 

}

// ================================================================================== //
// CALC MODES	                                                                      //
// ================================================================================== //

template <class T, class T1>
void POD<T, T1>::calc_modes(const bool &tf){
   
    int ft=0;
    vector<T> filter;
    vector<T> volume;
    vector<T> s1;
    vector<T> s2;
    vector<T> val;
    vector<T> mode;
    vector<vector<vector<double>>> phi;
    vector<vector<double>> help; 
  
    if (job().nmode<job().nsnap){
	npod=job().nmode;
    }
    else {
	npod=job().nsnap-1;
    }

    help.resize(frows,vector<double>(mean[0].size(),0.0));    
    phi.resize(frows,vector<vector<double>>(npod,vector<double>(mean[0].size(),0.0)));  

    if (tf == 1){

	job().get_fields(ft,"filter",filter);
	
	//solve eigenproblem
	calc_coefficients();

	for (int i=0; i<job().nsnap; ++i){
	    cout << "projecting snapshot " << i+1 << " of " << job().nsnap << endl;
	    job().get_fields(i,"snap",s1);
	    val=s1-mean;

	    for (int j=0; j<npod; ++j){
		for (int n=0; n<frows; ++n){
			phi[n][j]=phi[n][j]+((val[n]*filter[0])*coeff[n/ftype][j][i])/abs(lambda[n/ftype][j]);
			
		}

	    }	
	    
	}
	for (int pd=0; pd<npod; ++pd){
	    for (int n=0; n<frows; ++n){
		help[n]=phi[n][pd];
	    }

	    job().give_fields(pd,"modes",help);
	}
    }
    else {
	for ( int pd=0; pd<npod; ++pd){
	    for (int n=0; n<frows; ++n){
		job().get_fields(pd,"mode",mode);
		phi[n][pd]=mode[n];
	    }
	}
	job().get_input(coeff,"coefficients");
	job().get_input(lambda,"lambda");
    }
	
}

// calc coefficients ================================================================ //

template <class T, class T1>
void POD <T , T1>::calc_coefficients(){
	vector<vector<vector<double>>> Mcorr;
	vector<vector<vector<double>>> Mcorr_reduce;
		
	lapack_int info, n;

	double alambda[job().nsnap][1];
	double Marr[job().nsnap][job().nsnap];

	lambda.resize(dsize, vector<double>(job().nsnap,0.0));
	coeff.resize(dsize, vector<vector<double>>(npod, vector<double>(job().nsnap,0.0)));	


	n=job().nsnap;

	calc_correlation(Mcorr);

	//Mcorr_reduce.resize(Mcorr.size(),vector<vector<double>>(Mcorr[0].size(),vector<double>(Mcorr[0][0].size(),0.0)));
	Mcorr_reduce=Mcorr;

	for (int i=0; i<dsize; ++i){
	    for (int k=0; k<Mcorr_reduce[i].size(); ++k){
		for (int j=k; j<Mcorr_reduce[i].size(); ++j){
		    Marr[k][j]=Mcorr_reduce[i][k][j];
		    Marr[j][k]=Marr[k][j];
		}
	    }
	    cout << "solving eigenproblem nr " << i+1 << endl;
	    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, *Marr, n, *alambda);

	    for (int q=0; q<job().nsnap; ++q){
		lambda[i][q]=alambda[-q-1+job().nsnap][0];
	    }

	    for (int t=0; t<job().nsnap; ++t){
		for ( int m=0; m<npod; ++m){
			
		    coeff[i][m][t]= Mcorr_reduce[i][t][job().nsnap-m-1]*sqrt(abs(lambda[i][m]));
	
		}
	    }

	}

	job().give_output(coeff,"coefficients");
	job().give_output(lambda,"lambda");

}

// calc correlation ================================================================= //

template <class T, class T1>
void POD <T, T1>::calc_correlation(vector<vector<vector<double>>> &Mcorr){
	vector<T> snap1;
	vector<T> snap2;
	vector<double> corr;

	Mcorr.resize(dsize,vector<vector<double>>(job().nsnap,vector<double>(job().nsnap,0.0)));

	for (int i=0; i<job().nsnap; ++i){

	    cout << "row number " << i+1 << " of " << job().nsnap << endl;
	    job().get_fields(i,"snap",snap1);
	    snap1=snap1-mean;
	
	    for (int j=i; j<job().nsnap; ++j){

		job().get_fields(j,"snap",snap2);
		snap2=snap2-mean;

		inner_product(snap1, snap2, corr);
		
		for (int k=0; k<dsize; ++k){

		    Mcorr[k][i][j]=corr[k];
		    Mcorr[k][j][i]=corr[k];
		}

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
	int fi=0;
	corr.resize(dsize,0.0);

	job().get_fields(fi,"filter",filter);
	job().get_fields(fi,"volume",volume);

	if (job().flag.compare("vectorial")==0){

	    for (int i=0; i<dsize; ++i)
	    {
		sum((Dot_Product(s1[0+3*i],s2[0+3*i])+Dot_Product(s1[1+3*i],s2[1+3*i])+Dot_Product(s1[2+3*i],s2[2+3*i]))*volume,corr1);
		corr[i]=corr1;
	    }

	}
	else {

	    vector<T> vcorr=s1*s2; 
	    vector<T> vfactor=volume*filter;

	    for (int i=0; i<dsize; ++i){

		vcorr[i]=vcorr[i]*vfactor[0];
		sum(vcorr[i],corr1);
		corr[i]=corr1;
	    }

	}

}

// ================================================================================== //
// CALC PROJECTION MATRIX                                                             //
// ================================================================================== //

template <class T, class T1>
void POD<T, T1>::calc_projection_matrix(const bool &tf){
	
	lapack_int info, n, ipiv;
	n=npod;
	
	double mat[dsize][npod][npod];
	double idty[npod][npod];
	//lapack_int ipiv[npod];
	
	if (tf==1){
	    
		//info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, npod, npod, "double* mat", npod,"lapack_int* ipiv", "double* b", npod);

	}

}
