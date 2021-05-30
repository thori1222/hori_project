#pragma once

/**
	シンボル一覧：
	- agex.c
		- DIMX		状態変数の次元
		- DIMUC		制御入力の次元
				
		- tsim0		制御開始時刻
		- tsim		制御終了時刻

		- ht			制御周期
		- hdir		前進差分近似のステップ幅
		- tf			ホライズン長さ
		- dv			ホライズン分割数
		- zeta		continuationのフィードバックゲイン（論文のAs = -zeta * 単位行列）
		- kmax		GMRES反復回数

		- x0			初期状態
		- u0			制御入力初期値

		- tstart		計算開始時刻
		- tend		計算終了時刻
 **/
/// Springheadのベクトル・行列クラス
#include <Base/TQuaternion.h>

#include <boost/bind.hpp>

#include "vector.hpp"
#include "gmres.hpp"
#include "integrate.hpp"

template<int DIMX, int DIMUC>
class Model{
public:
	static const int DimX = DIMX;
	static const int DimU = DIMUC;

	using x_t = vector_t<DIMX>;		///< 状態
    using u_t =  vector_t<DIMUC>;	///< 制御入力
	
	struct xu_t : vector_t<DIMX + DIMUC>{	///< 状態と入力の複合型
		x_t& x(){ return *(x_t*)this; }
		u_t& u(){ return *(u_t*)((double*)this + DIMX); }
		const x_t& x()const{ return *(const x_t*)this; }
		const u_t& u()const{ return *(const u_t*)((const double*)this + DIMX); }
	};

	x_t    x0;		///< 初期状態
	u_t    u0;		///< 初期入力

	x_t    x;		///< 現在状態
	x_t    x1;		///< 次時刻の状態
	u_t    u;		///< 制御入力

	/*-------------- dPhi/dx -------------- */
	virtual void phix(double t, const x_t& x, x_t& phx1, int i) = 0;

	/*-------------- State Equation -------------- */
	virtual void xpfunc(double t, const x_t& x, const u_t& u, x_t& xprime, int i) = 0;

	/*-------------- Costate Equation -------------- */
	virtual void lpfunc(double t, const x_t& lmd, const xu_t& linp, x_t& lprime, int i) = 0; 

	/*-------------- Error in Optimality Condition, Hu -------------- */
	virtual void hufunc(double t, const x_t& x, const x_t& lmd, const u_t& u, u_t& hui) = 0;
};

template<class MODEL, int N, int KMAX>
class Controller{
public:
	static const int DimX = MODEL::DimX;
	static const int DimU = MODEL::DimU;
	static const int dv   = N;
	static const int kmax = KMAX;

	using controller_t = Controller<MODEL, N, KMAX>;
	using model_t = MODEL;
    using x_t = typename MODEL::x_t;
    using u_t = typename MODEL::u_t;
    using xu_t = typename MODEL::xu_t;

	/// 状態の時系列
	struct X_t : vector_t<DimX * (N+1)>{	
		template<class T> X_t& operator=(const T& v){ assign(v); return *this; }
		x_t& elem(int i){ return ((x_t*)this)[i]; }
		const x_t& elem(int i)const{ return ((const x_t*)this)[i]; }
	};
	/// 制御入力の時系列
	struct U_t : vector_t<DimU * N>{
		template<class T> U_t& operator=(const T& v){ this->assign(v); return *this; }
		u_t& elem(int i){ return ((u_t*)this)[i]; }
		const u_t& elem(int i)const { return ((const u_t*)this)[i]; }
	};

	double tf;		//< 予測ホライズンの長さ
	double ht;		//< 制御周期
	double alpha;
	double zeta;
	double hdir;
	double rtol;
	int    dstep;

	double htau;
	double onemhdir;
	double onephdir;
	double hdirbht;
	double onemzetahdir; 
	double ts;
	double tauf;

	x_t x1s;
	x_t	lmd0;
	u_t hu0;
	X_t xtau;
	X_t xtau1;
	X_t ltau;
	U_t duvec;
	U_t utau;
	U_t utau1;
	U_t hutau;
	U_t hutau1;
	U_t hutau2;
	U_t bvec;
	U_t utautmp;
	
	vector_t<dv  +1> tau; 
	vector_t<kmax+1> errvec;

	model_t*	model;

public:
	Controller(){}

	void init(model_t* m){
		model = m;

		onemhdir     = 1 - hdir / ht;
		onephdir	 = 1 + hdir / ht;
		hdirbht      = hdir / ht;
		onemzetahdir = 1 - zeta * hdir;
	
		int i;
		vector_t<DimU>   b;
		vector_t<DimU>   du0;
		vector_t<DimU+1> erru0;

		model->x = model->x0;
	
		// phi_x(x, t)
		model->phix(tsim0, model->x0, lmd0, i);

		// Hu( 
		model->hufunc(tsim0, model->x0, lmd0, model->u0, hu0);

		// 全時刻のuをu0で埋める
		model->u = model->u0;
		
		for(i=0; i<dv; i++){
			utau .elem(i) = model->u;
			hutau.elem(i) = hu0;
		}
		duvec.clear();
	}

	/*-------------- Control Update -------------- */

	// F(U,x,t)
	void errfunc(double t, const x_t& x, const U_t& u, U_t& hu)
	{
		int i;
		double taut;
		xu_t linp;

		// 予測ホライズンの時間刻み
		tauf = tf;
		htau = tauf / dv;

		// 初期時刻のx
		xtau.elem(0) = x;
	
		// 前進計算でxの時系列を計算（論文(1)式）
		x_t xd;
		for(taut = t, i=0; i < dv; taut += htau, i++){
			model->xpfunc(taut, xtau.elem(i), u.elem(i), xd, i);
			xtau.elem(i+1) = xtau.elem(i) + htau * xd;
			tau[i] = taut; 
		}
		tau[i] = taut; 

		// lambdaの終端条件（論文(7)式）
		model->phix(taut, xtau.elem(dv), ltau.elem(dv), dv);

		// 後退計算でlambdaの時系列を計算（論文(6)式）
		x_t ld;
		for(i = dv-1; i >= 0; i--){
			linp.x() = xtau.elem(i);
			linp.u() = u.elem(i);
			model->lpfunc(taut, ltau.elem(i+1), linp, ld, i);
			ltau.elem(i) = ltau.elem(i+1) - htau * ld;
			taut -= htau; 
			// 得られたxとlambdaよりHuを計算（論文(5)式）
			model->hufunc(taut, xtau.elem(i), ltau.elem(i+1), ((U_t&)u).elem(i), hu.elem(i));
		}
	}

	void adufunc(const U_t& du, U_t& adu){
		// U + h dU
		utau1 = utau + du * hdir;
		// F(U + hdU, x + h xd, t + h)
		errfunc(ts, x1s, utau1, hutau2);
		// F_U dU approx [F(U + h dU, x + h xd, t + h) - F(U, x + h xd, t + h)] / h
		adu = (hutau2 - hutau1) / hdir;
	}
	struct call_adufunc{
		controller_t* ctrl;
		void operator()(const U_t& du, U_t& adu){
			ctrl->adufunc(du, adu);
		}
		call_adufunc(controller_t* c):ctrl(c){}
	};

	void unew(double t, const x_t& x, const x_t& x1, u_t& u){
		// t + h
		ts = t + hdir;
		
		 //前進差分計算のためにxの時間微分xdが必要だが，実際の計算ではxdを求めないため，ここで差分で近似している
		x1s = hdirbht * x1 + x * onemhdir;
		
		// F(t, x, u)
		errfunc(t, x, utau, hutau);

		// F(t + h, x + h xd, u)
		errfunc(ts, x1s, utau, hutau1);

		// 論文(11)式のb
		// 論文では b = As F(U,x,t) - DhF(U,x,t:0,xd,1) であるが，
		// DhF(U,x,t:0,xd,1) = (F(U, x + h xd, t+h) - F(U,x,t))/h より
		// b = [(1 + h As) F(U,x,t) - F(U, x + h xd, t+h)] / h
		bvec = (hutau * onemzetahdir - hutau1) / hdir;

		// 連立方程式 F_U Ud = b をGMRESで解く
		// adufuncは渡されたvに関してF_U vを前進差分で求める関数
		nfgmres<call_adufunc, U_t, kmax>(call_adufunc(this), bvec, duvec, errvec);
		
		// U += dU * ht;
		utau += ht * duvec;

		// 更新されたu系列を入力としたときの予測状態の計算（次のステップでerrfuncで計算するが、保存＋比較のためにここでも計算しておく）
		//errfunc(t,x1,utau,hutau1);
		
		// model_t modeltmp;
		// for(int i = 0; i < N; i++)for(int j = 0; j < modeltmp.NCar; j++){
		// 	if(utau.elem(i)[j] > modeltmp.umax[j]){
		// 			utau.elem(i)[j] = modeltmp.umax[j];
		// 	}
		// 	else if(utau.elem(i)[j] < modeltmp.umin[j]){
		// 			utau.elem(i)[j] = modeltmp.umin[j];
		// 	}
		// 	duvec = (utau - utautmp)/ht;
		// }
		errfunc(t,x1,utau,hutau1);
		
		// MPCなので求めたuの時系列の先頭要素を返す
		u = utau.elem(0);
		utautmp = utau;
	}

	private:
        const double tsim0 = 0.0;
};
