template <class T, class T1>
POD <T, T1>::POD(){
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

    vector< T >	   snap;
    int fi = 0;

    if (tf == 1){

    	for (int i = 0; i < job().nsnap; i++){
	    job().get_fields(i,"snap", snap);
    	    if(i==0){
		mean=snap;
	    }
	    else{
	        mean=mean+snap;
	    }
	}

	for (int d=0; d<mean.size(); ++d){
	    mean[d]=mean[d]/(double(job().nsnap));
	}
	job().give_fields(fi,"mean",mean);
    }
    else{
	job().get_fields(fi,"mean", mean);
    }

    dsize=mean.size();

    return (mean); 

}

// ================================================================================== //
// CALC MODES	                                                                      //
// ================================================================================== //

template <class T, class T1>
void POD<T, T1>::calc_modes(const bool &tf){
   
    int i0=0;
    vector<vector<vector<double>>> Mcorr;
    vector<T> help;
    vector<T> mode;
  
    if (job().nmode<job().nsnap){
	npod=job().nmode;
    }
    else {
	npod=job().nsnap-1;
    }


    job().get_fields(i0,"filter",filter);
    job().get_fields(i0,"volume",volume);

    if (tf == 1){


	//solve eigenproblem
	calc_coefficients();
	phi=ass_phi(mean);

	help=mean;

	for (int pd=0; pd<npod; ++pd){
	    for (int d=0; d<dsize; ++d){
		help[d]=phi[pd][d];
	    }

	    job().give_fields(pd,"modes",help);
	}
    }
    else {
	for ( int pd=0; pd<npod; ++pd){
	for (int d=0; d<dsize; ++d){
	    job().get_fields(pd,"mode",mode);
	    phi[pd][d]=mode[d];
	}
	}

	job().get_input(coeff,"coefficients");
	job().get_input(lambda,"lambda");
    }
	
}

// ass_phi 2D ======================================================================= //

template <class T, class T1>
vector<vector<T>> POD <T , T1>::ass_phi(vector<vector<double>> &in){

	vector<T> s1;
	vector<T> val;
	vector<vector<T>> mphi;
	mphi.resize(npod,vector<vector<double>>(dsize,vector<double>(in[0].size(),0.0)));

	for (int i=0; i<job().nsnap; ++i){
	    cout << "projecting snapshot " << i+1 << " of " << job().nsnap << endl;
	    job().get_fields(i,"snap",s1);
	    val=s1-mean;
	    for (int j=0; j<npod; ++j){
	    for (int d=0; d<dsize; ++d){
	    for (int n=0; n<filter[0].size(); ++n)
		if (filter[0][n]==1){
		    mphi[j][d][n]=mphi[j][d][n]+(val[d][n]*coeff[d][j][i])/abs(lambda[d][j]);
		}	
	    }
	    } 
	}
	return(mphi);
} 

// ass_phi 3D ======================================================================= //

template <class T, class T1>
vector<vector<T>> POD <T , T1>::ass_phi(vector<vector<vector<double>>> &in){

	vector<T> s1;
	vector<T> val;
	vector<vector<T>> mphi;
	mphi.resize(npod,vector<vector<vector<double>>>(dsize,vector<vector<double>>(in[0].size(),vector<double>(in[0][0].size(),0.0))));

	for (int i=0; i<job().nsnap; ++i){
	    cout << "projecting snapshot " << i+1 << " of " << job().nsnap << endl;
	    job().get_fields(i,"snap",s1);
	    val=s1-mean;

	    for (int j=0; j<npod; ++j){
	    for (int d=0; d<dsize; ++d){
	    for (int n=0; n<filter[0].size(); ++n){
		if (filter[0][n][0]==1){
		    mphi[j][d][n]=mphi[j][d][n]+(val[d][n]*coeff[d][j][i])/abs(lambda[d][j]);
		}	
	    }
	    }
	    } 
	}
	return(mphi);
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
		
		    Mcorr_reduce[i][t][job().nsnap-m-1]=Marr[t][job().nsnap-m-1];	
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

		corr=inner_product(snap1, snap2);
		
		for (int k=0; k<dsize; ++k){

		    Mcorr[k][i][j]=corr[k];
		    Mcorr[k][j][i]=corr[k];
		}

	    }

	}
}

// inner product 2D ================================================================= //

template <class T, class T1>
vector<double> POD <T, T1>::inner_product(vector<vector<double>> &s1, vector<vector<double>> &s2){

	vector<double> corr;
	corr.resize(dsize,0.0);

	for (int d=0; d<dsize; ++d){	
	for (int n=0; n<filter[0].size(); ++n){
	    if (filter[0][n]==1){

		corr[d]=corr[d]+(s1[d][n]*s2[d][n]*volume[0][n]);

	    }
	}
	}

	return(corr);
}

// inner product 3D ================================================================= //

template <class T, class T1>
vector<double> POD <T, T1>::inner_product(vector<vector<vector<double>>> &s1, vector<vector<vector<double>>> &s2){

	vector<double> corr; 
	corr.resize(dsize,0.0);

	for (int d=0; d<dsize; ++d){	
	for (int n=0; n<filter[0].size(); ++n){
	    if (filter[0][n][0]==1){

		corr[d]=corr[d]+(Dot_Product(s1[d][n],s2[d][n])*volume[0][n][0]);

	    }
	}
	}

	return(corr);
}


// ================================================================================== //
// CALC PROJECTION MATRIX                                                             //
// ================================================================================== //

template <class T, class T1>
void POD<T, T1>::calc_projection_matrix(const bool &tf){

	vector<vector<vector<double>>> proj_mat;
	vector<vector<vector<double>>> proj_mat_reduce;
	lapack_int info;
	
	double mat[npod][npod];
	double idty[npod][npod];
	int ipiv[npod];

	if (tf==1){

	    proj_mat=ass_proj_matrix(mean);
	    proj_mat_reduce=proj_mat;
	    
	    for (int d=0; d<dsize; ++d){
		for (int i=0; i<npod; ++i){
		for (int j=0; j<npod; ++j){
		    mat[i][j]=0.0;
		    idty[i][j]=0.0;
		    mat[i][j]=proj_mat_reduce[d][i][j];
		}
		    idty[i][i]=1.0;
		}
		
		info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, npod, npod, *mat, npod, ipiv, *idty, npod);

		for (int i=0; i<npod; ++i){
		for (int j=0; j<npod; ++j){

		    proj_mat[d][i][j]=idty[i][j];

		}
		}

		job().give_output(proj_mat[d],"proj_matrix");

	    }
	    

	}

}

// ass_proj_matrix 2D =============================================================== //

template <class T, class T1>
vector<vector<vector<double>>> POD <T, T1>::ass_proj_matrix(vector<vector<double>> &s0){

	vector<vector<vector<double>>> proj_mat;
	proj_mat.resize(dsize,vector<vector<double>>(npod,vector<double>(npod,0.0)));

	for (int n=0; n<filter[0].size(); ++n){
	    if (filter[0][n]==1){
		for (int i=0; i<npod; ++i){
		for (int j=0; j<npod; ++j){
		for (int d=0; d<dsize; ++d){
		    proj_mat[d][i][j]=proj_mat[d][i][j]+phi[i][d][n]*phi[j][d][n]*volume[0][n];				
		}
		}
		}
	    }
	}

	return(proj_mat);
}

// ass_proj_matrix 3D =============================================================== //

template <class T, class T1>
vector<vector<vector<double>>> POD <T, T1>::ass_proj_matrix(vector<vector<vector<double>>> &s0){

	vector<vector<vector<double>>> proj_mat;
	proj_mat.resize(dsize,vector<vector<double>>(npod,vector<double>(npod,0.0)));

	for (int n=0; n<filter[0].size(); ++n){
	    if (filter[0][n][0]==1){
		for (int i=0; i<npod; ++i){
		for (int j=0; j<npod; ++j){
		for (int d=0; d<dsize; ++d){
		    proj_mat[d][i][j]=proj_mat[d][i][j]+Dot_Product(phi[i][d][n],phi[j][d][n])*volume[0][n][0];				
		}
		}
		}
	    }
	}

	return(proj_mat);
}

// ================================================================================== //
// CALC RECONSTRUCTION								      //
// ================================================================================== //

template <class T, class T1>
void POD<T, T1>::calc_reconstruction(const bool &tf){

	vector<T> snapo;
	vector<T> snapr;
	vector<vector<double>> coeff_pod;
	coeff_pod.resize(dsize,vector<double>(npod,0.0));
	int rs=job().nrecon;

	job().get_fields(rs,"snap",snapo);
	snapo=snapo-mean;

	coeff_pod=project(snapo);

	for (int d=0; d<dsize; ++d){

	    job().give_output(coeff_pod[d],"coeff pod");	
		
	}

	snapr=mean;
	reconstruct(coeff_pod,snapr);

	job().give_fields(rs,"reconstruction",snapr);	

}

// project 2D ======================================================================= //

template <class T, class T1>
vector<vector<double>> POD<T, T1>::project(vector<vector<double>> &so){

	vector<vector<double>> coeff_pod;
	vector<vector<double>> coeff_myid;
	coeff_pod.resize(dsize,vector<double>(npod,0.0));
	coeff_myid.resize(dsize,vector<double>(npod,0.0));

	for (int d=0; d<dsize; ++d){
	for (int n=0; n<filter[0].size(); ++n)
	{
	   if (filter[0][n]==1){
		for (int i=0; i<npod; ++i){
		    coeff_myid[d][i]=coeff_myid[d][i]+so[d][n]*phi[i][d][n]*volume[0][n];
		}
	   }
	}	    
	}

	coeff_pod=coeff_myid;
	return(coeff_pod);
}

// project 3D ======================================================================= //

template <class T, class T1>
vector<vector<double>> POD<T, T1>::project(vector<vector<vector<double>>> &so){

	vector<vector<double>> coeff_pod;
	vector<vector<double>> coeff_myid;
	coeff_pod.resize(dsize,vector<double>(npod,0.0));
	coeff_myid.resize(dsize,vector<double>(npod,0.0));

	for (int d=0; d<dsize; ++d){
	for (int n=0; n<filter[0].size(); ++n)
	{
	   if (filter[0][n][0]==1){
		for (int i=0; i<npod; ++i){
		    coeff_myid[d][i]=coeff_myid[d][i]+Dot_Product(so[d][n],phi[i][d][n])*volume[0][n][0];
		}
	   }
	}	    
	}

	coeff_pod=coeff_myid;
	return(coeff_pod);
}

// reconstruct 2D ===================================================================== //

template <class T, class T1>
void POD<T, T1>::reconstruct(vector<vector<double>> &cpod, vector<vector<double>> &sr){

	for (int d=0; d<dsize; ++d){
	for (int i=0; i<npod; ++i){
	for (int n=0; n<filter[0].size(); ++n){
	    if (filter[0][n]==1){
		sr[d][n]=sr[d][n]+cpod[d][i]*phi[i][d][n];
	    }
	}
	}	
	}

}

// reconstruct 3D ===================================================================== //

template <class T, class T1>
void POD<T, T1>::reconstruct(vector<vector<double>> &cpod, vector<vector<vector<double>>> &sr){

	for (int d=0; d<dsize; ++d){
	for (int i=0; i<npod; ++i){
	for (int n=0; n<filter[0].size(); ++n){
	    if (filter[0][n][0]==1){
		sr[d][n]=sr[d][n]+cpod[d][i]*phi[i][d][n];
	    }
	}
	}	
	}

}



