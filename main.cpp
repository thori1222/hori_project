#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <random>
#include <stdio.h>
#include <string.h>

using namespace std;

#include "platoon.hpp"


/*-------------- Global Variables -------------- */
double tsim0     = 0;			//< Simulation start time
double tsim      = 40;			//< Simulation end time

double J_1,J_2,J_3,J_4, J = 0;				//<Cost functions

#define NCAR   1	   		//< #cars in a CACC group
#define NGROUP 4			//< #groups
#define STEP   200          //< ステップ数(platoon.hppに合わせる！)

#define NTHW   1			//< #THW
#define NR	   1 			//< #R
#define NQ1	   1			//< #Q1
#define NQ2	   1 			//< #Q2
#define NQ3	   1			//< #Q3
#define NSF	   1			//< #Sf

double vl;		    		//< Leading vehicle speed
double dl;		    		//< Leading vehicle distance
double al;		    		//< Leading vehicle acceleration
double vf[NGROUP][NCAR];	//< Following vehicle speed
double df[NGROUP][NCAR];	//< Following vehicle distance
double af[NGROUP][NCAR];	//< Following vehicle acceleration

//G(s) = beta/(s+alpha), a(t+1) = (1/dt - alpha)*a(t) + beta*dt*u(t)
//(a(t+1) = (1/dt - 1/tau)*a(t) + dt*K/tau*u(t), a_dot(t) = -1/tau*a(t) + K/tau*u(t) ,tau = 1/alpha, K = beta/alpha)
double tau[NCAR];
double K[NCAR];	
double alpha[NCAR];
double beta[NCAR];	

double THW[NTHW] = {0.0};
int R[NR] = {80};
int Q1[NQ1] = {30};
int Q2[NQ2] = {1};
int Q3[NQ3] = {0};
int SF[NSF] = {100};

bool   Lagrange = false;   //ラグランジュ双対を行うかどうか//
double GAMMA    = 10.0;     //ラグランジュ双対のステップ幅//
double Rc       = 40.0;    //ラグランジュ双対の隊列長制約//
double epsilon  = 1.0;     //ラグランジュ双対の更新判定//

Platoon<NCAR>			platoon[NGROUP];
PlatoonController<NCAR>	controller[NGROUP];


//Velocity patterns
// double vel_pattern[21] = {
// 	  0.0, 20.0, 40.0, 60.0, 80.0, 80.0, 60.0, 40.0, 40.0, 40.0,
// 	 50, 60.0, 70.0, 80.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0, 30.0
// };
// double vel_pattern[21] = {
// 	  60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0,
// 	  60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0, 60.0
// };
// double vel_pattern[21] = {
// 	  15.0, 15.0, 25.0, 20.0, 30.0, 35.0, 30.0, 30.0, 20.0, 15.0,
// 	 10.0, 10.0, 15.0, 20.0, 20.0, 10.0, 25.0, 30.0, 30.0, 20.0, 20.0
// };
// double vel_pattern[21] = {
// 	 25.0, 25.0, 20.0, 30.0, 50.0, 70.0, 80.0, 80.0, 60.0, 30.0, 20.0,
// 	 10.0, 10.0, 20.0, 50.0, 60.0, 65.0, 60.0, 40.0, 30.0, 30.0
// };
double vel_pattern[21] = {
	  60.0, 60.0, 70.0, 80.0, 80.0, 90.0, 110.0, 110.0, 100.0, 95.0,
	 90.0, 80.0, 80.0, 60.0, 60.0, 70.0, 80.0, 70.0, 70.0, 60.0, 60.0
};

//Function for calculating velocity and acceleration of the leading vehicle
void calc_vel(double t, double* vel, double* acc){
	double h = 10.0; //----10.0----
	double k = 1000.0 / 3600.0; //< km/h -> m/s
	int n = 21;
	int idx = (int)floor(t/h);
	
	if(idx < 0){
		*vel = k * vel_pattern[0];
		*acc = 0.0;
	}
	else if(idx >= n){
		*vel = k * vel_pattern[n];
		*acc = 0.0;
	}
	else{
		double tau = t - (double)idx * h;
		*vel = k * ((1.0 - tau/h) * vel_pattern[idx] + tau/h * vel_pattern[idx+1]);
		*acc = k * (vel_pattern[idx+1] - vel_pattern[idx]) / h;
	}
}

void calc_cosvel(double t, double* vel, double* acc);
void calc_randvel(double t, double* vel, double* acc);
float random(int pct);

int sum_loop = NQ1*NQ2*NQ3*NSF*NR*NTHW;
int idx_loop = 0;
int pct_loop = 0;

/*----------------------main loop------------------------*/
int main(void){
	for(int ic = 0; ic < NCAR; ic++){      //対象の一次遅れ系モデルのパラメータ
		tau[ic] = 0.3;
		K[ic] = 0.9;

		alpha[ic] = 1/tau[ic];
		beta[ic] = K[ic]/tau[ic];
	}

	//ラグランジュ乗数の初期化
	if (Lagrange){
		for(int k = 0; k < STEP; k++) {
			for(int ig = 0; ig < NGROUP; ig++) platoon[ig].mu[k] = 0.0;
		}
	}

	for(int ig = 0; ig < NGROUP; ig++) platoon[ig].lagrange = Lagrange;      
	
	for(int iq1 = 0; iq1 < NQ1; iq1++){
	for(int iq2 = 0; iq2 < NQ2; iq2++){
	for(int iq3 = 0; iq3 < NQ3; iq3++){
	for(int isf = 0; isf < NSF; isf++){
	for(int ir = 0; ir < NR; ir++){
	for(int ithw = 0; ithw < NTHW; ithw++){
	
	pct_loop = (int)(idx_loop*100/sum_loop);
	idx_loop++;

	// //乱数生成関数の初期化（時刻で初期化を行うことで常に異なるパターンが生成）
	// srand((unsigned int)time(NULL));

	int    nsim = (int)((tsim - tsim0) / controller[0].ht);
	int	   isim;
	double t;

	clock_t start, end, calc_start;

	for(int ig = 0; ig < NGROUP; ig++)for(int ic = 0; ic < NCAR; ic++){
		platoon[ig].q[3*ic+0] = Q1[iq1];
		platoon[ig].q[3*ic+1] = Q2[iq2];
		platoon[ig].q[3*ic+2] = Q3[iq3];
		platoon[ig].sf[3*ic+0] = SF[isf]*Q1[iq1];
		platoon[ig].sf[3*ic+1] = SF[isf]*Q2[iq2];
		platoon[ig].sf[3*ic+2] = SF[isf]*Q3[iq3];
		platoon[ig].thw[ic] = THW[ithw];
		platoon[ig].r[ic] = R[ir];
		platoon[ig].alpha[ic] = alpha[ic];
		platoon[ig].beta[ic] = beta[ic];
	}
	cout<<"--------------------------"<<endl;


	printf(" Start  \n");
	start = clock();

	// 状態と制御入力の初期値を設定
	for(int i = 0; i < NGROUP; i++){
		controller[i].init(&platoon[i]);
	}
	
	// シミュレーションモデルの初期値
	dl = 0.0;
	calc_vel(0, &vl, &al);

	for(int ig = 0; ig < NGROUP; ig++)for(int ic = 0; ic < NCAR; ic++){
		df[ig][ic] = ((ig == 0 && ic == 0) ? dl : (ic == 0 ? df[ig-1][NCAR-1] : df[ig][ic-1])) - (platoon[ig].Ds + platoon[ig].thw[ic]*vl);
		platoon[ig].x[3*ic+1] = vl;
		vf[ig][ic] = platoon[ig].x[3*ic+1];
		af[ig][ic] = al;
	}

	//ファイル名指定
	char fn[200]="File";
	char temp[50];

	// sprintf(temp, "_Auto3carsCacc");
	// strcat(fn, temp);
	
	// sprintf(temp, "_tau%.2f",tau[0]);
	// strcat(fn, temp);

	// sprintf(temp, "_K%.2f",K[0]);
	// strcat(fn, temp);

	// sprintf(temp, "_THW%.1f", THW[ithw]);
	// strcat(fn, temp);

	// sprintf(temp,"_Q%i%i%i",Q1[0],Q2[0],Q3[0]);
	// strcat(fn, temp);

	// sprintf(temp, "_R%i_",R[ir]);
	// strcat(fn, temp);

	// sprintf(temp, "%iSf", SF[isf]);
	// strcat(fn, temp);

	// sprintf(temp, "_dconstEXPgain10");
	// strcat(fn, temp);

	// sprintf(temp, "_weight%i_dmin%.1f", platoon[0].sd[0], platoon[0].Dmin);
	// strcat(fn, temp);

	// sprintf(temp, "_uconstEXPgain10");
	// strcat(fn, temp);

	// sprintf(temp, "_weight%i_umax%.1f", platoon[0].su[0], platoon[0].umax[0]);
	// strcat(fn, temp);

	// sprintf(temp, "_zeta100_kmax10");
	// strcat(fn, temp);

	sprintf(temp, ".csv");
	strcat(fn, temp);

	//ファイル作成（t dl vl al df vf al x0 x1 x2 u J_1 J_2 J_3 J_4 J)
	ofstream fs(fn);
	//
	fs<<"t, dl, vl, al"<<",";
	for(int ig = 0; ig< NGROUP; ig++)for(int ic = 0; ic < NCAR; ic++){
		fs<<"df, vf, af"<<",";
	}
	for(int ig = 0; ig< NGROUP; ig++)for(int ic = 0; ic < NCAR; ic++){
		fs<<"x0, x1, x2, u"<<",";
	}
	fs << "J_1, J_2, J_3, J_4, J, F";

	// for (int ig = 0; ig < NGROUP; ig++)for (int ih = 0; ih < controller[ig].dv; ih++)for (int ic = 0; ic < NCAR; ic++) {
	// 	fs << "," << "h[" << ih << "]c[" << ig << "][" << ic << "]_x0";
	// 	fs << "," << "h[" << ih << "]c[" << ig << "][" << ic << "]_x1";
	// 	fs << "," << "h[" << ih << "]c[" << ig << "][" << ic << "]_x2";
	// 	fs << "," << "h[" << ih << "]c[" << ig << "][" << ic << "]_l0";
	// 	fs << "," << "h[" << ih << "]c[" << ig << "][" << ic << "]_l1";
	// 	fs << "," << "h[" << ih << "]c[" << ig << "][" << ic << "]_l2";
	// 	fs << "," << "h[" << ih << "]c[" << ig << "][" << ic << "]_u";
	// 	fs << "," << "h[" << ih << "]c[" << ig << "][" << ic << "]_hu";
	// }

	fs << endl;

	bool step_flag = false;
	/*----------------------Simulation Loop----------------------*/
    for(t = tsim0, isim = 0; isim < nsim; t += controller[0].ht, isim++){

		// 速度パターンより前方車速度と加速度を取得
		calc_vel(t, &vl, &al);

		// 走行距離更新
		dl += vl * controller[0].ht + 0.5*al*controller[0].ht*controller[0].ht; //等速度モデル

		// 自車速度，走行距離更新
		for(int ig = 0; ig < NGROUP; ig++)for(int ic = 0; ic < NCAR; ic++){
			df[ig][ic] +=  vf[ig][ic] * controller[ig].ht+ 0.5*af[ig][ic]*controller[0].ht*controller[0].ht;
			vf[ig][ic] += af[ig][ic]*controller[ig].ht;
			af[ig][ic] = beta[ic]*platoon[ig].u[ic]*controller[0].ht + (1-alpha[ic]*controller[0].ht)*af[ig][ic];
		
			//vfがマイナスになったらゼロにする（後進しないよう）<- 制約として入れるべき！
			if(vf[ig][ic] < 0) vf[ig][ic] = 0;
		}

		// モデルの状態設定
		for(int ig = 0; ig < NGROUP; ig++)for(int ic = 0; ic < NCAR; ic++){
			platoon[ig].x1[3*ic+0] = ((ig == 0 && ic == 0) ? dl : (ic == 0 ? df[ig-1][NCAR-1] : df[ig][ic-1])) - df[ig][ic];
			platoon[ig].x1[3*ic+1] = vf[ig][ic];
			platoon[ig].x1[3*ic+2] = af[ig][ic]>5.0?5.0:af[ig][ic]<-5.0?-5.0:af[ig][ic];
		}
		
		// MPCモデルの前方車加速度を設定
		for(int ig = 0; ig < NGROUP; ig++){
			//等速モデル
			//platoon[ig].a = 0; 
			//platoon[ig].v = (ig == 0 ? vl : platoon[ig-1].x[3*(NCAR-1)+1]);

			//等加速度モデル
			platoon[ig].a = (ig == 0 ? al : platoon[ig-1].x[3*(NCAR-1)+2]);
			platoon[ig].v = (ig == 0 ? vl : platoon[ig-1].x[3*(NCAR-1)+1]);
		}

		// ラグランジュパラメータ初期化（毎回μ=0スタート）
		if (Lagrange){
			for(int k = 0; k < STEP; k++) {
				for(int ig = 0; ig < NGROUP; ig++) platoon[ig].mu[k] = 0.0;
			}
		}
		int i = 0;
		double Gamma = GAMMA;

		// C/GMRESで制御入力を更新
		// 双対問題
		if(Lagrange){
			calc_start = clock();
			while(1) {
				double sum[STEP] = {0.0};

				//　入力を更新
				for(int ig = 0; ig < NGROUP; ig++){
					controller[ig].unew(t, platoon[ig].x, platoon[ig].x1, platoon[ig].u);
				}

				// xの時系列ごとの和を出す
				for(int k = 0; k < STEP; k++) {
					for(int ig = 0; ig < NGROUP; ig++) sum[k] += controller[ig].xtau.elem(k)[0];
				}
				
				// μの更新
				int j = 0;
				for(int ig = 0; ig < NGROUP; ig++){
					for(int k = 0; k < STEP; k++) {
						if (sum[k] - Rc > epsilon) {
							platoon[ig].mu[k] += Gamma * (sum[k] - Rc);
							j++;
						}
					}
				}
				if (j == 0) break;      //更新判定

				// 反復制限
				i++;
				// if (i >= 100) break;                                                      //最大反復回数
				if ((double)(clock() - calc_start)/CLOCKS_PER_SEC >= 0.01 * NGROUP) break;   //最大反復時間(s)

				// γの更新
				Gamma = GAMMA * 1 / sqrt(i);

				//確認用
				// std::cout << "i = " << i << " , γ = " << Gamma << std::endl;
				// double Sum = 0.0;
				// for(int k = 0; k < STEP; k++) 
				// Sum += sum[k] - Rc;
				// std::cout << "mu[" << k << "] = " << platoon[0].mu[k] << std::endl;
				// std::cout << "残差[" << k << "] = " << sum[k] - Rc << std::endl;
				// std::cout << "Sum = " << Sum << std::endl;

			}
		}
		// 双対でない場合
		else{
			for(int ig = 0; ig < NGROUP; ig++){
				controller[ig].unew(t, platoon[ig].x, platoon[ig].x1, platoon[ig].u);
			}
		}
		
		// 状態の更新
		for(int ig = 0; ig < NGROUP; ig++){
			platoon[ig].x = platoon[ig].x1;
		}

		//評価関数
		J_1 = 0;
		J_2 = 0;
		J_3 = 0;
		J_4 = 0;
		J   = 0;
		for(int ig = 0; ig < NGROUP; ig++)for(int ic =0; ic < NCAR; ic++){
			J_1 += 0.5* pow(platoon[ig].x[3*ic+0]-(platoon[ig].Ds+platoon[ig].thw[ic]*platoon[ig].x[3*ic+1]), 2.0); 
			J_2	+= 0.5* pow(platoon[ig].x[3*ic+1] - (ic == 0? (ig == 0? vl: vf[ig-1][NCAR]): vf[ig][ic-1]), 2.0);
			J_3 += 0.5* pow(platoon[ig].x[3*ic+2] - (ic == 0? (ig == 0? al: af[ig-1][NCAR]): af[ig][ic-1]), 2.0);
			J_4 += 0.5* pow(platoon[ig].u[ic], 2.0);
		}

		J += J_1 + J_2 + J_3 + J_4;

		//ファイル書き込み（t dl vl al df vf af x0 x1 x2 u y, w lap J_1 J_2 J_3 J_4 J)
		fs << t << ", " << dl << ", " << vl << "," << al;
		for(int ig = 0; ig < NGROUP; ig++)for(int ic = 0; ic < NCAR; ic++){
			fs << ", " << df[ig][ic] << ", " << vf[ig][ic] << ", " << af[ig][ic];
		}

		for(int ig = 0; ig < NGROUP; ig++)for(int ic = 0; ic < NCAR; ic++){
			fs << ", " << platoon[ig].x[3*ic+0] << ", " << platoon[ig].x[3*ic+1] << ", "<< platoon[ig].x[3*ic+2] << ", " << platoon[ig].u[ic];
		}

		fs << ", " << J_1 << ", " << J_2 << ", " << J_3 << ", " << J_4 << ", " << J; //endl;
		fs << "," << controller[0].hutau1.norm();

		// //予測ホライゾンの状態を保存
		// for (int ig = 0; ig < NGROUP; ig++)for (int ih = 0; ih < controller[ig].dv; ih++)for (int ic = 0; ic < NCAR; ic++) {
		// 	fs << " , " << controller[ig].xtau.elem(ih)[3 * ic + 0] <<
		// 		" , " << controller[ig].xtau.elem(ih)[3 * ic + 1] <<
		// 		" , " << controller[ig].xtau.elem(ih)[3 * ic + 2] <<
		// 		" , " << controller[ig].ltau.elem(ih)[3 * ic + 0] <<
		// 		" , " << controller[ig].ltau.elem(ih)[3 * ic + 1] <<
		// 		" , " << controller[ig].ltau.elem(ih)[3 * ic + 2] <<
		// 		" , " << (platoon[ig].manual[0] ? 0 : controller[ig].utau.elem(ih)[ic]) <<
		// 		" , " << controller[ig].hutau.elem(ih)[ic] <<
		// 		" , " << controller[ig].ptau.elem(ih)[ic];
		// }
	    fs<<endl;

		// 確認用 car1におけるcar0の速度予測の系列
		// std::cout << "-------------------------------------------" << std::endl;
		// cout << "-------------------------t = " << t << "-------------------------------" << endl;
		// double v1 = df[0][0];
		// double v2 = df[0][0];
		// for (int ih = 0; ih < controller[0].dv; ih++) {
		// 	v2 += controller[0].xtau.elem(ih)[1]*0.01;
		// 	cout << "ih : " << ih << "     =     " << ((controller[0].xtau.elem(ih+1)[0] + v2) - (controller[0].xtau.elem(ih)[0] + v1))/0.01 << endl;
		// 	v1 = v2;
		// }

		//時間の確認
		if(fmod(isim,100) == 0)
			cout<<"["<<pct_loop<<"% | "<<idx_loop<<"/"<<sum_loop<<"] Simulation time = "<<t<<"/"<<tsim<<endl;
		}
	end = clock();
	
	printf("%.2f秒かかりました\n",(double)(end-start)/CLOCKS_PER_SEC);
	}
	}
	}
	}
	}
	}
	getchar();
}



//正弦波で速度パターンを生成する関数
//< @2018.10.02
double v0 = 100.0/3.6; //< 初期速度
double amax = 1.0; //< 加速度振幅
double T = 50.0; //< 周期[sec]

void calc_cosvel(double t, double* vel, double* acc){
//速度関数：v(t) = v0 - amax*T/(2*pi)*cos(2*pi/T*t)
//加速度関数：a(t) = amax*sin(2*pi/T*t)
	/*if((int)t%10 == 0)
		T = (int)(100 - (9.0/20.0)*t);*/
	*vel = v0 - amax*T/(2*Spr::M_PI)*cos(2*Spr::M_PI/T*t);
	*acc = amax*sin(2*Spr::M_PI/T*t);
} //> @2018.10.02

//乱数バイナリでの加速度パターン生成関数
//< @2018.10.15
std::mt19937 engine(0);
std::uniform_real_distribution<> dist;
void calc_randvel(double t, double* vel, double* acc){
	double tmpvel = *vel;

	if(int(t*100)%10 == 0){
	double tmp = dist(engine);
		if(tmp > 0.66)
			*acc = 1.0;
		else if(tmp < 0.33)
			*acc = -1.0;
		else
			*acc = 0.0;
	}
	*vel = tmpvel + *acc * controller[0].ht;
} //> @2018.10.15

//乱数パーセンテージ生成；結果は1.02とか0.93など
float random(int pct){
	int level = 1;	//<	範囲の桁（１００なら０．０１％のレベルで生成）
	return 1 + (pct-(rand() % (pct*2*level) + 1)/(float)level) * 0.01;
}
