/*---- k Givens Roations used in GMRES ----*/
void givrot(int k, double *c, double *s, double *v)
{
	int i;
	double w1, w2;
	for(i=0; i<k; i++){
		w1 = c[i] * v[i] - s[i] * v[i+1];
		w2 = s[i] * v[i] + c[i] * v[i+1];
		v[i] = w1;
		v[i+1] = w2;
	}
}

/*--------------------------------------------------------
 Linear Equation Subroutine
	GMRES (Generalized Minimum Residual) 
	Cf. C.T.Kelley: Iterative Methods for Linear and Nonlinear Equations
						T.Ohtsuka  '00/01/24
	axfunc(int n, double *x, double *ax): a function that gives A*x
	n: dim x
	b: the right-hand-side of the equation, A*x = b
	x: initial guess, to be replaced with the solution
	kmax: number of iteration
	err: residual norms, |b - A*x_i| (i=1,...,n+1)
----------------------------------------------------------*/

/// gmres
template<class F, class X, int kmax>
//void nfgmres(F& axfunc, const X& b, X& x, vector_t<kmax+1>& err){
void nfgmres(F&& axfunc, const X& b, X& x, vector_t<kmax+1>& err){
	int i,j,k;
	double rho, nu; 
	
	vector_t<kmax+1> cvec; 
	vector_t<kmax+1> svec; 
	vector_t<kmax+1> gvec;
	vector_t<kmax+1> hmat[kmax+1]; 
	X tmpvec; 
	X vmat[kmax+1]; 

	// b - A * x0
	axfunc(x, tmpvec); 
	tmpvec = b - tmpvec;
	
	rho = tmpvec.norm();

	//try{
	//rho = tmpvec.norm();
	//if(rho > std::pow(10.0,100.0)){			
	//if(!_finite(rho))
	//	throw std::pow(10.0,300.0);
	//}}
	//catch(double _rho){
	//	rho = _rho;
	//}
	gvec.clear();
	gvec[0] = rho; 
	err[0] = rho;

	// �΍��𐳋K�������x�N�g��
	vmat[0] = tmpvec / rho;

	// kmax�񔽕�
	for(k = 0; k < kmax; k++){
		// vmat[k+1] = A * vmat[k];
		axfunc(vmat[k], vmat[k+1]); 

		/* Modified Gram-Schmidt */
		for(j = 0; j <= k; j++){
			hmat[k][j] = vmat[j] * vmat[k+1];
			vmat[k+1] -= vmat[j] * hmat[k][j];
		}
		hmat[k][k+1] = vmat[k+1].norm();
		
		/* No Breakdown? */
		if( hmat[k][k+1] != 0){
			vmat[k+1] /= hmat[k][k+1];
		}
		else{
			printf("gmress() : breakdown \n");
		}
		
		/* Givens Rotation */
		givrot(k, cvec, svec, hmat[k]); 
		nu = vec2_t(hmat[k][k], hmat[k][k+1]).norm();
		if( nu != 0 ){
			cvec[k] = hmat[k][k] / nu;
			svec[k] = - hmat[k][k+1] / nu;
			hmat[k][k] = cvec[k] * hmat[k][k] - svec[k] * hmat[k][k+1];
			hmat[k][k+1] = 0;
			givrot(1, cvec+k, svec+k, gvec+k);
		}
		else printf("nu is zero!\n");
		
		/* Residual Update */
		rho = fabs(gvec[k+1]);
		err[k+1] = rho;
	}
	
	/* Solve hmat * y = gvec (hmat: upper triangular) */
	for(i=k-1; i>=0; i--){
		for(nu=gvec[i], j=i+1; j<k; j++){
			nu -= hmat[j][i] * cvec[j];
		}
		cvec[i] = nu / hmat[i][i] ; 
	}

	/* Ans. */
	for(i=0; i<(int)x.size(); i++){
		for(nu=0, j=0; j<k; j++){
			nu += vmat[j][i] * cvec[j];
		}
		x[i] += nu;
	}
}
