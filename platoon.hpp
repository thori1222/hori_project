#ifndef PLATOON_HPP_
#define PLATOON_HPP_

#include <math.h>
#include <float.h>
#include "model.hpp"

#define STEP 200

template<int NCAR>
class Platoon : public Model<3*NCAR, 1*NCAR>{
	using model_t = Model<3*NCAR, 1*NCAR>;
public:
	static const int NCar = NCAR;
	
	int		refmode;
	double  cal_lim;

	bool    lagrange;    // true or false //
	double  mu[STEP];    // Lagrange multiplier //
	
	double  a;		    // Acceleration of the leading car
	double	v;			// Velocity of the leading car

	double	Ds;			// Standstill reference gap
	double  thw[NCAR];	// Time headway	
	double  Dmin;			// Minimum gap as safety constraint

	double  alpha[NCAR]; //(G(s) = beta/(s+alpha), a(t+1) = (1/dt - alpha)*a(t) + beta*dt*u(t))
	double  beta[NCAR]; //(a(t+1) = (1/dt - 1/tau)*a(t) + dt*K/tau*u(t), a_dot(t) = -1/tau*a(t) + K/tau*u(t) ,tau = 1/alpha, K = beta/alpha)

	double  umax[NCAR];		// Maximum acceleration
	double  umin[NCAR];		// Minimum acceleration

	int		sf    [NCAR*3];	// Weight on the final cost 
	int		q     [NCAR*3];	// Weight on the process cost
	int		r     [NCAR];	// Weight on the control input
	int		su	  [NCAR];   // Weight on the input constraint
	int		sd	  [NCAR];	// Weight on the distnace constraint

	
	Platoon(){
		refmode = 1;
		cal_lim = pow(10.0, 50.0);

		lagrange  = false; 

		a         = 0.0;
		v		  = 0.0;

		Ds        = 10.0;
		Dmin	  = 5.0;
		
		for(int i = 0; i < NCAR; i++){
			sf[3 * i + 0] = 1000;
			sf[3 * i + 1] = 100;
			sf[3 * i + 2] = 0;
			q[3*i+0]      = 10;
			q[3*i+1]      = 1;
			q[3*i+2]      = 0;
			r[i] = 50;
			su[i]		  = 100;
			sd[i]		  = 100;
			thw[i]        = 1.2;
			umax[i]	  = 3.0;
			umin[i]	  = -3.0;
			alpha[i] = 10.0;
			beta[i]  = 10.0; //tau = 0.1, K = 1;
		}

		for(int i = 0; i < NCAR; i++){
			this->x0[3*i + 0] =  Ds;
			this->x0[3*i + 1] =  0.0;
			this->x0[3*i + 2] =  0.0;
			this->u0[i]       =  0.000001;
		};
	}

	/*-------------- dPhi/dx -------------- */
	void phix(double t, const typename model_t::x_t& x, typename model_t::x_t& phx1)
	{
		for(int i = 0; i < NCAR; i++){
			phx1[3*i+0] = sf[3*i+0] * (x[3*i+0] - (Ds + thw[i]*x[3*i+1]));
			phx1[3*i+1] = phx1[3*i+0] * (-thw[i]) + sf[3*i+1]*(x[3*i+1]-(i==0 ? v : (refmode==0?v:x[3*(i-1)+1]))) + (refmode==0?0:(i == NCAR-1 ? 0 : -sf[3*(i+1)+1]*(x[3*(i+1)+1]-x[3*i+1])));
			phx1[3*i+2] = sf[3*i+2]*(x[3*i+2] - (i==0?0:(refmode==0?0:x[3*(i-1)+2]))) + refmode==0?0:(i==NCAR-1?0:-sf[3*(i+1)+1]*(x[3*(i+1)+2]-x[3*i+2]));
		}
	}

	/*-------------- State Equation -------------- */
	void xpfunc(double t, const typename model_t::x_t& x, const typename model_t::u_t& u, typename model_t::x_t& xprime)
	{
		for(int i = 0; i < NCAR; i++){
			xprime[3*i+0] = (i==0 ? v : x[3*(i-1)+1]) - x[3*i+1]; // + 0.5*0.01*(i==0?a:x[3*(i-1)+2]-x[3*i+2]); 
			xprime[3*i+1] = x[3*i+2];
			xprime[3*i+2] = beta[i]*u[i] - alpha[i]*x[3*i+2];
		}
	}

	/*-------------- Costate Equation -------------- */
	void lpfunc(double t, const typename model_t::x_t& lmd, const typename model_t::xu_t& linp, typename model_t::x_t& lprime, int j)
	{
		typename model_t::x_t    _x;
		typename model_t::u_t    _u;
		double tmp;
 
		_x = linp.x();
		_u = linp.u();
		
		// H_xにマイナスかけたもの
		for(int i = 0; i < NCAR; i++){
			// auto
			lprime[3*i+0] = -(q[3*i+0]*(_x[3*i+0]-(Ds+thw[i]*_x[3*i+1])));
			double tmp = 0;
			int gain = 10;
			if (_x[3 * i + 0] < Dmin + 2.0) {
					try {
					tmp = -(sd[i] * -gain * exp(-gain * (_x[3 * i + 0] - Dmin))); //Exponential function
					//tmp = -(sd[i] * (_x[3*i+0]<Dmin?-gain:0)); //Piece-wise linear function
					//tmp = -(sd[i] * (_x[3*i+0]<Dmin?4*(_x[3*i+0]-Dmin)*(_x[3*i+0]-Dmin)*(_x[3*i+0]-Dmin):0));

					if (tmp > cal_lim || tmp < -cal_lim)
						throw((tmp > 0 ? 1 : -1)*cal_lim);
				}
				catch (double _tmp) {
					tmp = _tmp;
				}
				lprime[3 * i + 0] += tmp;
			}
			if(lagrange) lprime[3 * i + 0] += (-mu[j]);   //とりあえず1 car per 1 groupのみ//
			lprime[3*i+1] = -(q[3*i+0]*(_x[3*i+0]-(Ds+thw[i]*_x[3*i+1]))*(-thw[i]) 
							- q[3*i+1]*((i==0 ? v : (refmode==0?v:_x[3*(i-1)+1])) -_x[3*i+1])
							+ (refmode==0?0:q[3*(i+1)+1]*(i==NCAR-1?0:_x[3*i+1] - _x[3*(i+1)+1])) 
							- lmd[3*i+0] + (i==NCAR-1?0:lmd[3*(i+1)+0]));
			lprime[3*i+2] = -(q[3*i+2]*(_x[3*i+2]-(i==0?a:refmode==0?0:_x[3*(i-1)+2])) + (refmode==0?0:q[3*(i+1)+2]*(i==NCAR-1?0:_x[3*i+2]-_x[3*(i+1)+2]))                     
							// - lmd[3*i+0]*0.5*0.01 + (i==NCAR-1?0:lmd[3*(i+1)+0]*0.5*0.01) 
							+ lmd[3*i+1] - lmd[3*i+2]*alpha[i]);
		}
	}

	/*-------------- Error in Optimality Condition, Hu -------------- */
	void hufunc(double t,
            const typename model_t::x_t& x,
            const typename model_t::x_t& lmd,
            const typename model_t::u_t& u,
            typename model_t::u_t& hui)
	{
		for(int i = 0; i < NCAR; i++){
			hui[i] = r[i]*u[i] + lmd[3*i+2]*beta[i];

			double tmp = 0;
			int gain = 10;
			try{
				tmp = su[i]*(gain*exp(gain*(u[i]-umax[i]))-gain*exp(-gain*(u[i]-umin[i]))); //Exponential function
				//tmp = su[i]*(u[i]>umax?gain:(u[i]<-umax?-gain:0)); //Piece-wise linear function
				//tmp = su[i]*4*u[i]*u[i]*u[i]*gain; 
				//tmp = su[i]*(u[i]>umax?4*u[i]*u[i]*u[i]:(u[i]>0?gain:(u[i]>-umax?-gain:4*u[i]*u[i]*u[i]))); 
				if(tmp>cal_lim || tmp<-cal_lim)
					throw((tmp>0?1:-1)*cal_lim);
			}
			catch(double _tmp){
				tmp = _tmp;
			}
			hui[i] += tmp;
		}
	}

};

template<int NCAR>
class PlatoonController : public Controller< Platoon<NCAR>, STEP, 10>{   //(ステップ数，GMRESの反復回数) 3くらいかも
public:
	PlatoonController(){
		this->tf     = 2.0;
		this->ht     = HT;
		this->alpha  = 1.5;
		this->zeta   = 50;
		this->hdir   = 1.e-8;//0.002;
		this->rtol   = 1.e-6;
		this->dstep  = 1;
	}
};
#endif