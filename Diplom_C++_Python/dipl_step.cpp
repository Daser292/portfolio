#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <math.h>

#define PI 3.14159265 

using namespace std;

double W (double x,double h){
	float f = abs(x) / h;
	float N = 0.7 * PI * pow(h, 2);
	if(f < 1){
		return (1 - 3 * pow(f, 2) / 2 + 3 * pow(f, 3) / 4) / N;
	}
	else if(f < 2){
		return pow(2 - f, 3) / (4 * N);
	}
	else{
		return 0;
	}
}

double gradW (double f_xy, double h){
	float f = abs(f_xy) / h;
	float N = 28 * PI * pow(h, 3);
	if(f < 1 and f > 0){
		return (-12 * f + 9 * pow(f, 2)) / N;
	}
	else if(f < 2 and f >= 1){
		return (-3 * pow(2 - f, 2)) / N;
	}
	else{
		return 0;
	}
}

double skalar (double a_x, double a_y, double b_x, double b_y){
	return a_x * b_x + a_y * b_y;
}

double rast (double a_x, double b_x, double a_y, double b_y){
	return sqrt (pow (a_x - b_x, 2) + pow (a_y - b_y, 2));
}

void files (string f, int k, int num_all, int num_all_f, double mass[], double mass_f[]){
	ofstream fl;
	fl.open ("./out/" + f + to_string (k) +".txt", fstream::out);
	for(int a = 0; a < num_all_f; a++){
		fl << to_string (mass_f [a]) << '\n';
	}
	for (int a = 0; a < num_all; a++){
		fl << to_string (mass [a]) << '\n';
	}
	fl.close ();
}

double sreda (double h, double v_x_1, double v_x_2, double v_y_1, double v_y_2, double k_x_1, double k_x_2, double k_y_1, double k_y_2, double p_1, double p_2){
	double alfa = 1.0, betta = 0.5, eps = 0.01, M;
	M = h * skalar (v_x_1 - v_x_2, v_y_1 - v_y_2, k_x_1 - k_x_2, k_y_1 - k_y_2) / (pow (rast (k_x_1, k_x_2, k_y_1, k_y_2), 2) + (eps * h * h));
	if (M < 0){
		return ((-alfa) * 10 * M + betta * M * M) / (0.5 * (p_1  + p_2));
	}
    else{
        return 0;
	} 
}

double energy (double en [], double d [], double d_f [], double m [], double p [], double v_x [], double v_y [], double k_x [], double k_y [], double v_x_f [], double v_y_f [], double k_x_f [], double k_y_f [], double m_f [], double p_f [], int i, double t, double h, int num_all, int num_all_f, int s){
	double e = 0;
	for(int j = 0; j < num_all; j++){
		if(rast (k_x [j], k_x [i], k_y [j], k_y [i]) < (h * 2) && 
		   rast (k_x [j], k_x [i], k_y [j], k_y [i]) > 0){
			e = e + ((d [i]) / (p [i] * p [i]) + 
					 (d [j]) / (p [j] * p [j]) + 
					  0.5 * s * sreda (h, v_x [i], v_x [j], 
					  	                  v_y [i], v_y [j], 
					  	                  k_x [i], k_x [j], 
					  	                  k_y [i], k_y [j], 
					  	                  p [i], p [j])) * 
					        0.5 * m [j] * 
							skalar (v_x [i] - v_x [j], 
							       v_y [i] - v_y [j], 
							       k_x [i] - k_x [j],
							       k_y [i] - k_y [j]) *
							gradW (rast (k_x [i], k_x [j], k_y [i], k_y [j]), h) *
							(1.0 / rast (k_x [i], k_x [j], k_y [i], k_y [j]));
		}
	}
	for(int j = 0; j < num_all_f; j++){
		if(rast(k_x_f [j], k_x [i], k_y_f [j], k_y [i]) < (h * 2)){
			e = e + ((d [i]) / (p [i] * p [i]) + 
					 (d_f [j]) / (p_f [j] * p_f [j]) +
					  0.5 * s * sreda (h, v_x [i], v_x_f [j], 
					  	                  v_y [i], v_y_f [j], 
					  	                  k_x [i], k_x_f [j], 
					  	                  k_y [i], k_y_f [j], 
					  	                  p [i], p_f [j])) *
					  	         0.5 * m_f [j] * 
								skalar (v_x [i] - v_x_f [j], 
								       v_y [i] - v_y_f [j], 
								       k_x [i] - k_x_f [j],
								       k_y [i] - k_y_f [j]) *
								gradW (rast (k_x [i], k_x_f [j], k_y [i], k_y_f [j]), h) *
								(1.0 / rast (k_x [i], k_x_f [j], k_y [i], k_y_f [j]));
		}
	}
	return e * t + en [i];
}

double skotost_x (double d[], double m[], double p[], double v_x[], double v_y[], double k_x[], double k_y[], double d_f[], double v_x_f [], double v_y_f [], double k_x_f [], double k_y_f [], double m_f [], double p_f [], int i, double t, double h, int num_all, int num_all_f, int s){
	double v = 0;
	for(int j = 0; j < num_all; j++){
		if(rast (k_x [j], k_x [i], k_y [j], k_y [i]) < (h * 2) && 
		   rast (k_x [j], k_x [i], k_y [j], k_y [i]) > 0 &&
		   abs(k_x [j] - k_x [i]) > 0){
			v = v - m [j] * 
					(d [i] / (p [i] * p [i]) + 
					 d [j] / (p [j] * p [j]) +
					 s * sreda (h, v_x [i], v_x [j], 
		  	                   v_y [i], v_y [j], 
		  	                   k_x [i], k_x [j], 
		  	                   k_y [i], k_y [j], 
		  	                   p [i], p [j])) *  
					gradW (rast (k_x [j], k_x [i], k_y [j], k_y [i]), h) *
					(k_x [i] - k_x [j]) / rast (k_x [j], k_x [i], k_y [j], k_y [i]);
		}
	}
	for(int j = 0; j < num_all_f; j++){
		if(rast(k_x_f [j], k_x [i], k_y_f [j], k_y [i]) < (h * 2) &&
			abs(k_x_f [j] - k_x [i]) > 0){
			v = v - m_f [j] * 
					(d [i] / (p [i] * p [i]) + 
					 d_f [j] / (p_f [j] * p_f [j]) +
					 s * sreda (h, v_x [i], v_x_f [j], 
			  	                   v_y [i], v_y_f [j], 
			  	                   k_x [i], k_x_f [j], 
			  	                   k_y [i], k_y_f [j], 
			  	                   p [i], p_f [j])) *  
					gradW (rast(k_x_f [j], k_x [i], k_y_f [j], k_y [i]), h) *
					(k_x [i] - k_x_f [j]) / rast(k_x_f [j], k_x [i], k_y_f [j], k_y [i]);
		}
	}
	return v * t + v_x [i];
}

double skotost_y (double d[], double m[], double p[], double v_x[], double v_y[], double k_x[], double k_y[], double d_f[], double v_x_f [], double v_y_f [], double k_x_f [], double k_y_f [], double m_f [], double p_f [], int i, double t, double h, int num_all, int num_all_f, int s){
	double v = 0;
	for(int j = 0; j < num_all; j++){
		if(rast (k_x [j], k_x [i], k_y [j], k_y [i]) < (h * 2) && 
		   rast (k_x [j], k_x [i], k_y [j], k_y [i]) > 0 &&
		   abs(k_y [j] - k_y [i]) > 0){
			v = v - m [j] * 
					(d [i] / (p [i] * p [i]) + 
					 d [j] / (p [j] * p [j]) +
					 s * sreda (h, v_x [i], v_x [j], 
		  	                   v_y [i], v_y [j], 
		  	                   k_x [i], k_x [j], 
		  	                   k_y [i], k_y [j], 
		  	                   p [i], p [j])) * 
					gradW (rast (k_x [j], k_x [i], k_y [j], k_y [i]), h) *
					(k_y [i] - k_y [j]) / rast (k_x [j], k_x [i], k_y [j], k_y [i]);
		}
	}
	
	for(int j = 0; j < num_all_f; j++){
		if(rast(k_x_f [j], k_x [i], k_y_f [j], k_y [i]) < (h * 2) &&
			abs(k_y_f [j] - k_y [i]) > 0){
			v = v - m_f [j] * 
					(d [i] / (p [i] * p [i]) + 
					 d_f [j] / (p_f [j] * p_f [j]) +
					 s * sreda (h, v_x [i], v_x_f [j], 
			  	                   v_y [i], v_y_f [j], 
			  	                   k_x [i], k_x_f [j], 
			  	                   k_y [i], k_y_f [j], 
			  	                   p [i], p_f [j])) *  
					gradW (rast(k_x_f [j], k_x [i], k_y_f [j], k_y [i]), h) *
					(k_y [i] - k_y_f [j]) / rast(k_x_f [j], k_x [i], k_y_f [j], k_y [i]);
		}
	}
	return v * t + v_y [i];
}

double plotnost (double m [], double pl [], double v_x [], double v_y [], double k_x [], double k_y [], double v_x_f [], double v_y_f [], double k_x_f [], double k_y_f [], double m_f [], double p_f [], int i, double t, double h, int num_all, int num_all_f, int s){ 
	double p=0;
	for(int j = 0; j < num_all; j++){
		if(rast (k_x [j], k_x [i], k_y [j], k_y [i]) < (h * 2) && 
		   rast (k_x [j], k_x [i], k_y [j], k_y [i]) > 0){
			p = p + m[j] * (pl [i] / pl [j]) *
					skalar( v_x [i] - v_x [j],
							v_y [i] - v_y [j], 
							k_x [i] - k_x [j],
							k_y [i] - k_y [j]) * 
							gradW (rast (k_x [i], k_x [j], k_y [i], k_y [j]), h) *
							(1.0 / rast (k_x [i], k_x [j], k_y [i], k_y [j]));
		}
	}
	for(int j = 0; j < num_all_f; j++){
		if(rast(k_x_f [j], k_x [i], k_y_f [j], k_y [i]) < (h * 2)){
			p = p + m_f [j] * (pl [i] / p_f [j]) *
					skalar (v_x [i] - v_x_f [j], 
							v_y [i] - v_y_f [j], 
							k_x [i] - k_x_f [j],
							k_y [i] - k_y_f [j]) *
					gradW (rast (k_x [i], k_x_f [j], k_y [i], k_y_f [j]), h) *
					(1.0 / rast (k_x [i], k_x_f [j], k_y [i], k_y_f [j]));
		}
	}
	return p * t + pl [i];
}

double tempreture (double m [], double p [], double k_x [], double k_y [], double k_x_f [], double k_y_f [], double m_f [], double p_f [], double tp [], double tp_f [], double d_f [], double d [], int i, double t, double h, int num_all, int num_all_f){ 
	double tmp = 0;
	double c = 1; //???
	for(int j = 0; j < num_all; j++){
		if(rast (k_x [j], k_x [i], k_y [j], k_y [i]) < (h * 2) && 
		   rast (k_x [j], k_x [i], k_y [j], k_y [i]) > 0){
			tmp = tmp + 1 / (p [i] * c) * (m [j] / p [j]) * 
						(4 * d [i] * d [j]) / (d [i] + d [j]) * (tp [i] - tp [j]) * 
				        gradW (rast (k_x [i], k_x [j], k_y [i], k_y [j]), h) *
						(1.0 / rast (k_x [i], k_x [j], k_y [i], k_y [j]));
		}
	}
	for(int j = 0; j < num_all_f; j++){
		if(rast(k_x_f [j], k_x [i], k_y_f [j], k_y [i]) < (h * 2)){
			tmp = tmp + 1 / (p [i] * c) * (m_f [j] / p_f [j]) * 
						(4 * d [i] * d_f [j]) / (d [i] + d_f [j]) * (tp [i] - tp_f [j]) * 
				        gradW (rast (k_x [i], k_x_f [j], k_y [i], k_y_f [j]), h) *
						(1.0 / rast (k_x [i], k_x_f [j], k_y [i], k_y_f [j]));
		}
	}
	return tmp * t + tp [i];
}

double davlenie (double p, double e){
	return (0.4) * p * e;
}

double en (double p, double d){
	return d / ((0.4) * p);
}

double tmp (double p, double d){
	return d / ((8.31) * p);
}

int main(){
	ofstream x_y_num;

	int py = 1; //запустить питон?
	double h; //расстояние взаимодействия
	double t_step; //шаг по времени
	int t_steps; //количество шагов по времени
	int sreda = 1; //вязкость

	t_step = 0.001;
	h = 0.1;
	t_steps = 5;

	//начальные значения
	int num_y = 50;
	int num_left = 10;
	int num_x = 3 * num_y;

	double r = 1.0 / double (num_y);

	int num_h = round(h * 2 / r);

	int num_all = num_y * num_y * 3 - round((0.2 / 1) * num_y) * round(2.4 / 3 * num_x) + num_left * num_y;
	int num_all_f = (num_left + num_x + num_h) * (num_y + num_h * 2) - num_all;

	//перенос данных в питон
	x_y_num.open("./out/x_y_num.txt", fstream::out);
	x_y_num << to_string(num_all + num_all_f) << '\n';
	x_y_num.close();

	//основные
	double koord_x_mass [num_all], 	koord_y_mass [num_all], velosity_x_mass [num_all], velosity_y_mass [num_all];
	double davlenie_mass [num_all], plotnost_mass [num_all], energy_mass [num_all], mass_mass [num_all];

	//фиктивные
	double koord_x_mass_f [num_all_f], koord_y_mass_f [num_all_f];
	double energy_mass_f [num_all_f], mass_mass_f [num_all_f];
	double velosity_x_mass_f [num_all_f], velosity_y_mass_f [num_all_f];
	double davlenie_mass_f [num_all_f], plotnost_mass_f [num_all_f];

	//на следущем шаге 
	double koord_x_mass_next [num_all], koord_y_mass_next [num_all], velosity_x_mass_next [num_all], velosity_y_mass_next [num_all];
	double energy_mass_next [num_all], davlenie_mass_next [num_all], plotnost_mass_next [num_all];

	int osn = 0;
	int fict = 0;

	//определение начальных данных
	//
	for(int i = 0; i < num_x; i++){
		for(int j = 0; j < num_y; j++){
			double kor_x = i * r;
			double kor_y = j * r;
			if(kor_x < 0.6 || kor_y >= 0.2){
				koord_x_mass [osn] = kor_x;
				koord_y_mass [osn] = kor_y;
				plotnost_mass [osn] = 1.4;
				davlenie_mass [osn] = 1.0;
				velosity_y_mass [osn] = 0.0;
				velosity_x_mass [osn] = 0.0001;
				energy_mass [osn] = en(1.4, 1.0);
				mass_mass [osn] = 1.4 / (pow(r, 2));
				osn++;
			}
			else{
				koord_x_mass_f [fict] = kor_x;
				koord_y_mass_f [fict] = kor_y;
				plotnost_mass_f [fict] = 1.4;
				davlenie_mass_f [fict] = 1.0;
				velosity_y_mass_f [fict] = 0.0;
				velosity_x_mass_f [fict] = 0.0;
				energy_mass_f [fict] = en(1.4, 1.0);
				mass_mass_f [fict] = 1.4 / (pow(r, 2));
				fict++;
			}
		}
	}

	for(int i = 0; i < num_left + num_x + num_h; i++){
		for(int j = 0; j < num_y + num_h * 2; j++){
			double kor_x = (i - num_left) * r;
			double kor_y = (j - num_h) * r;
			if(kor_x < 0 || kor_y < 0 || kor_x >= 3 || kor_y >= 1){
				if(kor_x < 0 && kor_y < 1 && kor_y >= 0){
					koord_x_mass [osn] = kor_x;
					koord_y_mass [osn] = kor_y;
					plotnost_mass [osn] = 1.4;
					davlenie_mass [osn] = 1.0;
					velosity_y_mass [osn] = 0.0;
					velosity_x_mass [osn] = 0.0001;
					energy_mass [osn] = en(1.4, 1.0);
					mass_mass [osn] = 1.4 / (pow(r, 2));
					osn++;
				}
				else if(kor_x >= 3 && kor_y < 1 && kor_y >= 0.2){
					koord_x_mass_f [fict] = kor_x;
					koord_y_mass_f [fict] = kor_y;
					plotnost_mass_f [fict] = 1.4;
					davlenie_mass_f [fict] = 1.0;
					velosity_y_mass_f [fict] = 0.0;
					velosity_x_mass_f [fict] = 0.0001;
					energy_mass_f [fict] = en(1.4, 1.0);
					mass_mass_f [fict] = 1.4 / (pow(r, 2));
					fict++;
				}
				else{
					koord_x_mass_f [fict] = kor_x;
					koord_y_mass_f [fict] = kor_y;
					plotnost_mass_f [fict] = 1.4;
					davlenie_mass_f [fict] = 1.0;
					velosity_y_mass_f [fict] = 0.0;
					velosity_x_mass_f [fict] = 0.0;
					energy_mass_f [fict] = en(1.4, 1.0);
					mass_mass_f [fict] = 1.4 / (pow(r, 2));
					fict++;
				}
			}
		}
	}

	//шаги
	for(int k = 0; k < t_steps; k++){
		files("k_x_", k, num_all, num_all_f, koord_x_mass, koord_x_mass_f);
		files("k_y_", k, num_all, num_all_f, koord_y_mass, koord_y_mass_f);
		files("e_", k, num_all, num_all_f, energy_mass, energy_mass_f);
		files("v_x_", k, num_all, num_all_f, velosity_x_mass, velosity_x_mass_f);
		files("v_y_", k, num_all, num_all_f, velosity_y_mass, velosity_y_mass_f);
		files("d_", k, num_all, num_all_f, davlenie_mass, davlenie_mass_f);
		files("p_", k, num_all, num_all_f, plotnost_mass, plotnost_mass_f);

		//вычисление по шагам
		for(int i = 0; i < num_all; i++){
			if(koord_x_mass [i] >= 0 && koord_x_mass [i] < 3 && koord_y_mass [i] >= 0 && koord_y_mass [i] < 1){
				plotnost_mass_next [i] = plotnost (mass_mass, plotnost_mass, velosity_x_mass, velosity_y_mass, koord_x_mass, koord_y_mass, velosity_x_mass_f, velosity_y_mass_f, koord_x_mass_f, koord_y_mass_f, mass_mass_f, plotnost_mass_f, i, t_step, h, num_all, num_all_f, sreda);
				energy_mass_next [i] = energy (energy_mass, davlenie_mass, davlenie_mass_f, mass_mass, plotnost_mass, velosity_x_mass, velosity_y_mass, koord_x_mass, koord_y_mass,velosity_x_mass_f, velosity_y_mass_f, koord_x_mass_f, koord_y_mass_f, mass_mass_f, plotnost_mass_f, i, t_step, h, num_all, num_all_f, sreda);
				velosity_x_mass_next [i] = skotost_x (davlenie_mass, mass_mass, plotnost_mass, velosity_x_mass, velosity_y_mass, koord_x_mass, koord_y_mass, davlenie_mass_f, velosity_x_mass_f, velosity_y_mass_f, koord_x_mass_f, koord_y_mass_f, mass_mass_f, plotnost_mass_f, i, t_step, h, num_all, num_all_f, sreda);
				velosity_y_mass_next [i] = skotost_y (davlenie_mass, mass_mass, plotnost_mass, velosity_x_mass, velosity_y_mass, koord_x_mass, koord_y_mass, davlenie_mass_f, velosity_x_mass_f, velosity_y_mass_f, koord_x_mass_f, koord_y_mass_f, mass_mass_f, plotnost_mass_f, i, t_step, h, num_all, num_all_f, sreda);
				davlenie_mass_next [i] = davlenie (plotnost_mass [i], energy_mass [i]);
			}
		}

		//смена координат
		for(int a = 0; a < num_all; a++){
				koord_x_mass_next [a] = velosity_x_mass [a] * t_step + koord_x_mass [a];
				koord_y_mass_next [a] = velosity_y_mass [a] * t_step + koord_y_mass [a];
		}
		for(int a = 0; a < num_all_f; a++){
			koord_x_mass_f [a] = velosity_x_mass_f [a] * t_step + koord_x_mass_f [a];
			koord_y_mass_f [a] = velosity_y_mass_f [a] * t_step + koord_y_mass_f [a];
		}

		//перенос в основные массивы
		for(int i = 0; i < num_all; i++){
			if(koord_x_mass [i] >= 0 && koord_x_mass [i] < 3 && koord_y_mass [i] >= 0 && koord_y_mass [i] < 1){
				plotnost_mass [i] = plotnost_mass_next [i];
				energy_mass [i] = energy_mass_next [i];
				davlenie_mass [i] = davlenie_mass_next [i];
				velosity_x_mass [i] = velosity_x_mass_next [i];
				velosity_y_mass [i] = velosity_y_mass_next [i];
				koord_x_mass [i] = koord_x_mass_next [i];
				koord_y_mass [i] = koord_y_mass_next [i];
			}
		}
	}

	if (py){
		system("python graph_step.py");
	}
 
	return 0;
}
