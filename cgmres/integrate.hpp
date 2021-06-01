/// integration
typedef void (*model_func_t)(double t, double x[], double u[], double xprime[]);

template<class X, class U>
struct model_t{
	model_func_t func;
	void operator()(double t, const X& x, const U& u, X& xp){
		func(t, (double*)&x, (double*)&u, (double*)&xp);
	}
	model_t(model_func_t f){ func = f; }
};

template<class X, class U>
void nfrkginpex(model_func_t func, double t, const X& x, const U& u, double h, X& ans){
	X fval;
	X k1, k2, k3, k4;
	X xp1, xp2, xp3;
	X q1, q2, q3;

	const double c1 = 0.2928932188134528;
	const double c2 = 0.1213203435596426;
	const double c3 = 0.5857864376269054;
	const double c4 = 1.707106781186548;
	const double c5 = 4.121320343559643;
	const double c6 = 3.414213562373097;

	model_t<X, U> F(func);

	F(t, x, u, fval);

	xp1 = x + (0.5 * h) * fval;
	q1  = h * fval;

	F(t + 0.5 * h, xp1, u, fval);

	xp2 = xp1 + c1 * (h * fval - q1);
	q2 = c2 * q1 + c3 * h * fval;

	F(t + 0.5 * h, xp2, u, fval);

	xp3 = xp2 + c4 * (h * fval - q2);
	q3 = -c5 * q2 + c6 * h * fval;

	F(t + h, xp3, u, fval);

	ans = xp3 + h * fval / 6.0 - q3 / 3.0;
}

template<class X, class U>
void nfeulerinp(model_func_t func, double t, const X& x, const U& u, double h, X& ans){
	X fval;
	model_t<X, U> F(func);
	F(t, x, u, fval);
	ans = x + h * fval;
}
