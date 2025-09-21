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

void files_pole (string f, int k, int num_all, double mass[]){
	ofstream fl;
	fl.open ("./out/" + f + to_string (k) +".txt", fstream::out);
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
	system("del \"/out/\"");

	ofstream x_y_num;

	int py; //запустить питон?
	double h; //расстояние взаимодействия
	double t_step; //шаг по времени
	int t_steps; //количество шагов по времени
	int x_num_pole; //x точек поля
	int y_num_pole; //y точек поля
	int pole; //вывод точек поля
	int change_koordin = 1; //изменение координат
	int sreda = 1; //вязкость
	int d_3; //3х мерный график
	double p_ne, p_nw, p_sw, p_se;
	double d_ne, d_nw, d_sw, d_se;
	double v_x_ne, v_x_nw, v_x_sw, v_x_se;
	double v_y_ne, v_y_nw, v_y_sw, v_y_se;
	int num_f_ne, num_f_nw, num_f_sw, num_f_se;
	int num_ne, num_nw, num_sw, num_se;
	int num_all, num_all_f;

	t_step = 0.001;
	h = 0.1;
	t_steps = 3;

	p_ne = 1.5;
	p_nw = 0.5323;
	p_sw = 0.138;
	p_se = 0.5323;

	d_ne = 1.5;
	d_nw = 0.3;
	d_sw = 0.029;
	d_se = 0.3;

	v_x_ne = 0.0;
	v_x_nw = 1.206;
	v_x_sw = 1.206;
	v_x_se = 0.0;

	v_y_ne = 0.0;
	v_y_nw = 0.0;
	v_y_sw = 1.206;
	v_y_se = 1.206;

	num_all = 10000;

	x_num_pole = 50;
	y_num_pole = 50;

	pole = 0;

	py = 1;
	d_3 = 1;

	double p_all = p_ne + p_nw + p_sw + p_se;

	num_ne = int(sqrt((p_ne / p_all) * num_all));
	num_nw = int(sqrt((p_nw / p_all) * num_all));
	num_sw = int(sqrt((p_sw / p_all) * num_all));
	num_se = int(sqrt((p_se / p_all) * num_all));

	num_all = pow(num_ne, 2) + pow(num_nw, 2) + pow(num_sw, 2) + pow(num_se, 2);

	num_f_ne = int(h * 4 * num_ne);
	num_f_nw = int(h * 4 * num_nw);
	num_f_sw = int(h * 4 * num_sw);
	num_f_se = int(h * 4 * num_se);

	num_all_f = pow(num_ne + num_f_ne, 2) + pow(num_nw + num_f_nw, 2) + pow(num_sw + num_f_sw, 2) + pow(num_se + num_f_se, 2) - num_all;

	cout << num_f_ne << ' ' << num_f_nw << ' ' << num_f_sw << ' ' << num_f_se << endl;
	cout << num_ne << ' ' << num_nw << ' ' << num_sw << ' ' << num_se << endl;
	cout << num_all << ' ' << num_all_f << endl;

	if (pole == 0){
		x_num_pole = 0;
		y_num_pole = 0;
	}

	//перенос данных в питон
	x_y_num.open("./out/x_y_num.txt", fstream::out);
	x_y_num << to_string(num_all + num_all_f) << '\n';
	x_y_num << to_string(x_num_pole * y_num_pole) << '\n';
	x_y_num << to_string(pole) << '\n';
	x_y_num.close();

	//основные
	double koord_x_mass [num_all], koord_y_mass [num_all], velosity_x_mass [num_all], velosity_y_mass [num_all];
	double davlenie_mass [num_all], plotnost_mass [num_all], energy_mass [num_all], mass_mass [num_all];

	//фиктивные
	double koord_x_mass_f [num_all_f], koord_y_mass_f [num_all_f];
	double energy_mass_f [num_all_f], mass_mass_f [num_all_f];
	double velosity_x_mass_f [num_all_f], velosity_y_mass_f [num_all_f];
	double davlenie_mass_f [num_all_f], plotnost_mass_f [num_all_f];

	//на следущем шаге 
	double koord_x_mass_next [num_all], koord_y_mass_next [num_all], velosity_x_mass_next [num_all], velosity_y_mass_next [num_all];
	double energy_mass_next [num_all], davlenie_mass_next [num_all], plotnost_mass_next [num_all];

	//поле
	double koord_x_mass_p [x_num_pole * y_num_pole], koord_y_mass_p [x_num_pole * y_num_pole];
	double energy_mass_p [x_num_pole * y_num_pole];
	double velosity_x_mass_p [x_num_pole * y_num_pole], velosity_y_mass_p [x_num_pole * y_num_pole];
	double davlenie_mass_p [x_num_pole * y_num_pole], plotnost_mass_p [x_num_pole * y_num_pole];

	int osn = 0;
	int fict = 0;

	//определение начальных данных
	//заполение частиц sw
	for(int i = 0; i < num_sw + num_f_sw; i++){
		for(int j = 0; j < num_sw + num_f_sw; j++){
			double tmp_y = double (i - num_f_sw) * (0.5 / num_sw);
			double tmp_x = double (j - num_f_sw) * (0.5 / num_sw);
			if (tmp_x < 0 || tmp_y < 0){
				plotnost_mass_f [fict] = p_sw;
				energy_mass_f [fict] = en(p_sw, d_sw);
				davlenie_mass_f [fict] = d_sw;
				mass_mass_f [fict] = p_sw / (double(num_sw) * 4 * double(num_sw));
				velosity_x_mass_f [fict] = v_x_sw;
				velosity_y_mass_f [fict] = v_y_sw;
				koord_x_mass_f [fict] = tmp_x;
				koord_y_mass_f [fict] = tmp_y;
				fict++;
			}
			else{
				plotnost_mass [osn] = p_sw;
				energy_mass [osn] = en(p_sw, d_sw);
				davlenie_mass [osn] = d_sw;
				mass_mass [osn] = p_sw / (double(num_sw) * 4 * double(num_sw));
				velosity_x_mass [osn] = v_x_sw;
				velosity_y_mass [osn] = v_y_sw;
				koord_x_mass [osn] = tmp_x;
				koord_y_mass [osn] = tmp_y;
				osn++;
			}
		}
	}

	//заполение частиц nw
	for(int i = 0; i < num_nw + num_f_nw; i++){
		for(int j = 0; j < num_nw + num_f_nw; j++){
			double tmp_y = 0.5 + double (i) * (0.5 / num_nw);
			double tmp_x = double (j - num_f_nw) * (0.5 / num_nw);
			if (tmp_x < 0 || tmp_y >= 1){
				plotnost_mass_f [fict] = p_nw;
				energy_mass_f [fict] = en(p_nw, d_nw);
				davlenie_mass_f [fict] = d_nw;
				mass_mass_f [fict] = p_nw / (double(num_nw) * 4 * double(num_nw));
				velosity_x_mass_f [fict] = v_x_nw;
				velosity_y_mass_f [fict] = v_y_nw;
				koord_x_mass_f [fict] = tmp_x;
				koord_y_mass_f [fict] = tmp_y;
				fict++;
			}
			else{
				plotnost_mass [osn] = p_nw;
				energy_mass [osn] = en(p_nw, d_nw);
				davlenie_mass [osn] = d_nw;
				mass_mass [osn] = p_nw / (double(num_nw) * 4 * double(num_nw));
				velosity_x_mass [osn] = v_x_nw;
				velosity_y_mass [osn] = v_y_nw;
				koord_x_mass [osn] = tmp_x;
				koord_y_mass [osn] = tmp_y;
				osn++;
			}
		}
	}

	//заполение частиц se
	for(int i = 0; i < num_se + num_f_se; i++){
		for(int j = 0; j < num_se + num_f_se; j++){
			double tmp_y = double (i - num_f_se) * (0.5 / num_se);
			double tmp_x = 0.5 + double (j) * (0.5 / num_se);
			if (tmp_x >= 1 || tmp_y < 0){
				plotnost_mass_f [fict] = p_se;
				energy_mass_f [fict] = en(p_se, d_se);
				davlenie_mass_f [fict] = d_se;
				mass_mass_f [fict] = p_se / (double(num_se) * 4 * double(num_se));
				velosity_x_mass_f [fict] = v_x_se;
				velosity_y_mass_f [fict] = v_y_se;
				koord_x_mass_f [fict] = tmp_x;
				koord_y_mass_f [fict] = tmp_y;
				fict++;
			}
			else{
				plotnost_mass [osn] = p_se;
				energy_mass [osn] = en(p_se, d_se);
				davlenie_mass [osn] = d_se;
				mass_mass [osn] = p_se / (double(num_se) * 4 * double(num_se));
				velosity_x_mass [osn] = v_x_se;
				velosity_y_mass [osn] = v_y_se;
				koord_x_mass [osn] = tmp_x;
				koord_y_mass [osn] = tmp_y;
				osn++;
			}
		}
	}

	//заполение частиц ne
	for(int i = 0; i < num_ne + num_f_ne; i++){
		for(int j = 0; j < num_ne + num_f_ne; j++){
			double tmp_y = 0.5 + double (i) * (0.5 / num_ne);
			double tmp_x = 0.5 + double (j) * (0.5 / num_ne);
			if (tmp_x >= 1 || tmp_y >= 1){
				plotnost_mass_f [fict] = p_ne;
				energy_mass_f [fict] = en(p_ne, d_ne);
				davlenie_mass_f [fict] = d_ne;
				mass_mass_f [fict] = p_ne / (double(num_ne) * 4 * double(num_ne));
				velosity_x_mass_f [fict] = v_x_ne;
				velosity_y_mass_f [fict] = v_y_ne;
				koord_x_mass_f [fict] = tmp_x;
				koord_y_mass_f [fict] = tmp_y;
				fict++;
			}
			else{
				plotnost_mass [osn] = p_ne;
				energy_mass [osn] = en(p_ne, d_ne);
				davlenie_mass [osn] = d_ne;
				mass_mass [osn] = p_ne / (double(num_ne) * 4 * double(num_ne));
				velosity_x_mass [osn] = v_x_ne;
				velosity_y_mass [osn] = v_y_ne;
				koord_x_mass [osn] = tmp_x;
				koord_y_mass [osn] = tmp_y;
				osn++;
			}
		}
	}

	//определение координат поля
	for(int i = 0; i < y_num_pole; i++){
		for(int j = 0; j < x_num_pole; j++){
			koord_x_mass_p [i * y_num_pole + j] = j / double(x_num_pole);
			koord_y_mass_p [i * y_num_pole + j] = i / double(y_num_pole);
			plotnost_mass_p [i * y_num_pole + j] = 0;
			velosity_x_mass_p [i * y_num_pole + j] = 0;
			velosity_y_mass_p [i * y_num_pole + j] = 0;
			energy_mass_p [i * y_num_pole + j] = 0;
			davlenie_mass_p [i * y_num_pole + j] = 0;
		}
	}

	for(int k = 0; k < t_steps; k++){
		//вычисление поля	
		if(pole){
			for(int i = 0; i < y_num_pole; i++){
				for(int j = 0; j < x_num_pole; j++){
					plotnost_mass_p [i * y_num_pole + j] = 0;
					velosity_x_mass_p [i * y_num_pole + j] = 0;
					velosity_y_mass_p [i * y_num_pole + j] = 0;
					energy_mass_p [i * y_num_pole + j] = 0;
					davlenie_mass_p [i * y_num_pole + j] = 0;
				}
			}
			for(int i = 0; i < y_num_pole; i++){
				for(int j = 0; j < x_num_pole; j++){
					for(int a = 0; a < num_all; a++){
						if(rast(koord_x_mass [a], koord_x_mass_p [i * y_num_pole + j], koord_y_mass [a], koord_y_mass_p [i * y_num_pole + j]) < (h * 2)){
							double r = rast(koord_x_mass [a], koord_x_mass_p [i * y_num_pole + j], koord_y_mass [a], koord_y_mass_p [i * y_num_pole + j]);
							plotnost_mass_p [i * y_num_pole + j] += mass_mass [a] * W(r, h);
							velosity_x_mass_p [i * y_num_pole + j] += mass_mass [a] * (velosity_x_mass [a] / plotnost_mass [a]) * W(r, h);
							velosity_y_mass_p [i * y_num_pole + j] += mass_mass [a] * (velosity_y_mass [a] / plotnost_mass [a]) * W(r, h);
							energy_mass_p [i * y_num_pole + j] += (mass_mass [a] / plotnost_mass [a]) * energy_mass [a] * W(r, h);
						}
						davlenie_mass_p [i * y_num_pole + j] = davlenie (plotnost_mass_p [i * y_num_pole + j], energy_mass_p [i * y_num_pole + j]);
					}
					for(int a = 0; a < num_all_f; a++){
						if(rast(koord_x_mass_f [a], koord_x_mass_p [i * y_num_pole + j], koord_y_mass_f [a], koord_y_mass_p [i * y_num_pole + j]) < (h * 2)){
							double r = rast(koord_x_mass_f [a], koord_x_mass_p [i * y_num_pole + j], koord_y_mass_f [a], koord_y_mass_p [i * y_num_pole + j]);
							plotnost_mass_p [i * y_num_pole + j] += mass_mass_f [a] * W(r, h);
							velosity_x_mass_p [i * y_num_pole + j] += mass_mass_f [a] * (velosity_x_mass_f [a] / plotnost_mass_f [a]) * W(r, h);
							velosity_y_mass_p [i * y_num_pole + j] += mass_mass_f [a] * (velosity_y_mass_f [a] / plotnost_mass_f [a]) * W(r, h);
							energy_mass_p [i * y_num_pole + j] += (mass_mass_f [a] / plotnost_mass_f [a]) * energy_mass_f [a] * W(r, h);
						}
						davlenie_mass_p [i * y_num_pole + j] = davlenie (plotnost_mass_p [i * y_num_pole + j], energy_mass_p [i * y_num_pole + j]);
					}
				}
			}
		}
		if(pole == 1){
			files_pole("k_x_", k, x_num_pole * y_num_pole, koord_x_mass_p);
			files_pole("k_y_", k, x_num_pole * y_num_pole, koord_y_mass_p);
			files_pole("e_", k, x_num_pole * y_num_pole, energy_mass_p);
			files_pole("v_x_", k, x_num_pole * y_num_pole, velosity_x_mass_p);
			files_pole("v_y_", k, x_num_pole * y_num_pole, velosity_y_mass_p);
			files_pole("d_", k, x_num_pole * y_num_pole, davlenie_mass_p);
			files_pole("p_", k, x_num_pole * y_num_pole, plotnost_mass_p);
		}
		else{
			files("k_x_", k, num_all, num_all_f, koord_x_mass, koord_x_mass_f);
			files("k_y_", k, num_all, num_all_f, koord_y_mass, koord_y_mass_f);
			files("e_", k, num_all, num_all_f, energy_mass, energy_mass_f);
			files("v_x_", k, num_all, num_all_f, velosity_x_mass, velosity_x_mass_f);
			files("v_y_", k, num_all, num_all_f, velosity_y_mass, velosity_y_mass_f);
			files("d_", k, num_all, num_all_f, davlenie_mass, davlenie_mass_f);
			files("p_", k, num_all, num_all_f, plotnost_mass, plotnost_mass_f);
		}
		

		//вычисление по шагам
		for(int i = 0; i < num_all; i++){
			plotnost_mass_next [i] = plotnost (mass_mass, plotnost_mass, velosity_x_mass, velosity_y_mass, koord_x_mass, koord_y_mass, velosity_x_mass_f, velosity_y_mass_f, koord_x_mass_f, koord_y_mass_f, mass_mass_f, plotnost_mass_f, i, t_step, h, num_all, num_all_f, sreda);
			energy_mass_next [i] = energy (energy_mass, davlenie_mass, davlenie_mass_f, mass_mass, plotnost_mass, velosity_x_mass, velosity_y_mass, koord_x_mass, koord_y_mass,velosity_x_mass_f, velosity_y_mass_f, koord_x_mass_f, koord_y_mass_f, mass_mass_f, plotnost_mass_f, i, t_step, h, num_all, num_all_f, sreda);
			velosity_x_mass_next [i] = skotost_x (davlenie_mass, mass_mass, plotnost_mass, velosity_x_mass, velosity_y_mass, koord_x_mass, koord_y_mass, davlenie_mass_f, velosity_x_mass_f, velosity_y_mass_f, koord_x_mass_f, koord_y_mass_f, mass_mass_f, plotnost_mass_f, i, t_step, h, num_all, num_all_f, sreda);
			velosity_y_mass_next [i] = skotost_y (davlenie_mass, mass_mass, plotnost_mass, velosity_x_mass, velosity_y_mass, koord_x_mass, koord_y_mass, davlenie_mass_f, velosity_x_mass_f, velosity_y_mass_f, koord_x_mass_f, koord_y_mass_f, mass_mass_f, plotnost_mass_f, i, t_step, h, num_all, num_all_f, sreda);
			davlenie_mass_next [i] = davlenie (plotnost_mass [i], energy_mass [i]);

		}

		//смена координат
		if(change_koordin){
			for(int a = 0; a < num_all; a++){
					koord_x_mass_next [a] = velosity_x_mass [a] * t_step + koord_x_mass [a];
					koord_y_mass_next [a] = velosity_y_mass [a] * t_step + koord_y_mass [a];
			}
			for(int a = 0; a < num_all_f; a++){
				koord_x_mass_f [a] = velosity_x_mass_f [a] * t_step + koord_x_mass_f [a];
				koord_y_mass_f [a] = velosity_y_mass_f [a] * t_step + koord_y_mass_f [a];
			}
		}

		//перенос в основные массивы
		for(int i = 0; i < num_all; i++){
			plotnost_mass [i] = plotnost_mass_next [i];
			energy_mass [i] = energy_mass_next [i];
			davlenie_mass [i] = davlenie_mass_next [i];
			velosity_x_mass [i] = velosity_x_mass_next [i];
			velosity_y_mass [i] = velosity_y_mass_next [i];
			if(change_koordin){ 
				koord_x_mass [i] = koord_x_mass_next [i];
				koord_y_mass [i] = koord_y_mass_next [i];
			}
		}
	}

	if (py){
		if(d_3 == 1){
			system("python graph_2d.py");
		}
		else{
			system("python graph.py");
		}
	}
 
	return 0;
}
