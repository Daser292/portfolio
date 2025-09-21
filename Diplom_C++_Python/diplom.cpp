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

void files (string f, int k, int x_num, int y_num, int x_num_f, int y_num_f, double mass[], double mass_f[]){
	ofstream fl;
	fl.open ("./out/" + f + to_string (k) +".txt", fstream::out);
	for(int a = 0; a < (y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num; a++){
		fl << to_string (mass_f [a]) << '\n';
	}
	for (int a = 0; a < x_num * y_num; a++){
		fl << to_string (mass [a]) << '\n';
	}
	fl.close ();
}

void files_pole (string f, int k, int x_num, int y_num, double mass[]){
	ofstream fl;
	fl.open ("./out/" + f + to_string (k) +".txt", fstream::out);
	for (int a = 0; a < x_num * y_num; a++){
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

double energy (double en [], double d [], double d_f [], double m [], double p [], double v_x [], double v_y [], double k_x [], double k_y [], double v_x_f [], double v_y_f [], double k_x_f [], double k_y_f [], double m_f [], double p_f [], int i, double t, double h, int y_num, int x_num, int x_num_f, int y_num_f, int s){
	double e = 0;
	for(int j = 0; j < y_num * x_num; j++){
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
	for(int j = 0; j < (y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num; j++){
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

double skotost_x (double d[], double m[], double p[], double v_x[], double v_y[], double k_x[], double k_y[], double d_f[], double v_x_f [], double v_y_f [], double k_x_f [], double k_y_f [], double m_f [], double p_f [], int i, double t, double h, int y_num, int x_num, int x_num_f, int y_num_f, int s){
	double v = 0;
	for(int j = 0; j < y_num * x_num; j++){
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
	for(int j = 0; j < (y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num; j++){
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
	return v * t + v_x [i];
}

double skotost_y (double d[], double m[], double p[], double v_x[], double v_y[], double k_x[], double k_y[], double d_f[], double v_x_f [], double v_y_f [], double k_x_f [], double k_y_f [], double m_f [], double p_f [], int i, double t, double h, int y_num, int x_num, int x_num_f, int y_num_f, int s){
	double v = 0;
	for(int j = 0; j < y_num * x_num; j++){
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
	
	for(int j = 0; j < (y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num; j++){
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

double plotnost (double m [], double pl [], double v_x [], double v_y [], double k_x [], double k_y [], double v_x_f [], double v_y_f [], double k_x_f [], double k_y_f [], double m_f [], double p_f [], int i, double t, double h, int y_num, int x_num, int x_num_f, int y_num_f, int s){ 
	double p=0;
	for(int j = 0; j < y_num * x_num; j++){
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
	for(int j = 0; j < (y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num; j++){
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

double tempreture (double m [], double p [], double k_x [], double k_y [], double k_x_f [], double k_y_f [], double m_f [], double p_f [], double tp [], double tp_f [], double d_f [], double d [], int i, double t, double h, int y_num, int x_num, int x_num_f, int y_num_f){ 
	double tmp = 0;
	double c = 1; //???
	for(int j = 0; j < y_num * x_num; j++){
		if(rast (k_x [j], k_x [i], k_y [j], k_y [i]) < (h * 2) && 
		   rast (k_x [j], k_x [i], k_y [j], k_y [i]) > 0){
			tmp = tmp + 1 / (p [i] * c) * (m [j] / p [j]) * 
						(4 * d [i] * d [j]) / (d [i] + d [j]) * (tp [i] - tp [j]) * 
				        gradW (rast (k_x [i], k_x [j], k_y [i], k_y [j]), h) *
						(1.0 / rast (k_x [i], k_x [j], k_y [i], k_y [j]));
		}
	}
	for(int j = 0; j < (y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num; j++){
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

	double y, x; //координаты желаемого прямоугольника
	int y_num, x_num; //частицы по осям
	double P_L, P_R, D_L, D_R, V_L_Y, V_L_X, V_R_Y, V_R_X; //начения на разрыве
	int scatter; //по точкам?
	int py; //запустить питон?
	double h; //расстояние взаимодействия
	double t_step; //шаг по времени
	int t_steps; //количество шагов по времени
	int x_num_f, y_num_f; //количество фиктивных по x, y
	int lines; //вывод линий
	int x_num_pole; //x точек поля
	int y_num_pole; //y точек поля
	int pole; //вывод точек поля
	int change_koordin; //изменение координат
	int sreda; //вязкость
	int f; // фиктивные
	int kol_cir, kol_cir_num, cir; //разлет шара
	int d_3; //3х мерный график
	int soty; //распределение фиктивных как соты
	int num_r = 0; //количество частиц справа от разрыва
	int num_l = 0; //количество частиц слева от разрыва
	int num_r_f = 0; //количество фиктивных частиц справа от разрыва
	int num_l_f = 0; //количество фиктивных частиц справа от разрыва

	//тесты
	// Тест - 0, Торо - 1, Дуэль - 2, Распад шара - 3, Разрыв на плоскости - 4

	/*

	ПРОВЕРЯЙ КОЭФФИЦИЕНТ N ПРИ ОДНОМЕРНЫХ ТЕСТАХ

	x y лучше не менять

	*/

	int test = 0;
	t_step = 0.001;
	h = 0.04;
	t_steps = 300;
	x = 1;
	y = 1;

	//Тестовый набор
	if (test == 0){
		P_L = 1; 
		P_R = 0.125; 
		D_L = 1; 
		D_R = 0.1;  
		V_L_Y = 0;
		V_L_X = 0; 
		V_R_Y = 0; 
		V_R_X = 0; 
		x_num = 300;
		y_num = 30;
		x_num_pole = 50;
		y_num_pole = 50;
		lines = 0;
		change_koordin = 1;
		pole = 1;
		sreda = 1;
		f = 1;
		soty = 0;
		py = 1;
		d_3 = 1;
		scatter = 1;
		cir = 0;

		y_num_f = int((h * 2) * double(y_num) / y);
		x_num_f = int((h * 2) * double(x_num) / x);

		//количество частиц справа и слева
		if (P_L > P_R){
			double kof = double(P_L) / P_R + 1;
			num_r = int((double(x_num) / kof));
			num_l = int((double(x_num) / kof) * (double(P_L) / P_R));
		}
		else{
			double kof = double(P_R) / P_L + 1;
			num_r = int((double(x_num) / kof) * (double(P_R) / P_L));
			num_l = int((double(x_num) / kof));
		}

		num_r_f = int((h * 2) / (0.5 / num_r));
		num_l_f = int((h * 2) / (0.5 / num_l));

		x_num = num_l + num_r;
		x_num_f = num_l_f + num_r_f;

		cout << num_l << ' ' << num_r << ' ' << num_l_f << ' ' << num_r_f << endl;;
	}
	//первый тест Торо
	if (test == 1){
		P_L = 1; 
		P_R = 0.125; 
		D_L = 1; 
		D_R = 0.1;  
		V_L_Y = 0;
		V_L_X = 0; 
		V_R_Y = 0; 
		V_R_X = 0; 
		x_num = 400;
		y_num = 1;
		x_num_pole = 10;
		y_num_pole = 10;
		lines = 0;
		change_koordin = 1;
		pole = 0;
		sreda = 1;
		f = 1;
		soty = 0;
		py = 1;
		d_3 = 1;
		scatter = 1;
		cir = 0;
		y_num_f = int((h * 2) * double(y_num) / y);
		x_num_f = int((h * 2) * double(x_num) / x);
	}

	//дуэль
	if (test == 2){
		P_L = 1; 
		P_R = 1; 
		D_L = 1; 
		D_R = 1;  
		V_L_Y = 0;
		V_L_X = 0; 
		V_R_Y = 0; 
		V_R_X = -3; 		
		x_num = 2;
		y_num = 1;
		x_num_f = 0;
		y_num_f = 0;
		lines = 0;
		change_koordin = 1;
		pole = 0;
		sreda = 0;
		f = 0;
		soty = 0;
		py = 1;
		d_3 = 0;
		scatter = 1;
		cir = 0;
	}

	//распад газового шара
	if (test == 3){
		P_L = 1; 
		P_R = 1; 
		D_L = 1; 
		D_R = 1;  
		V_L_Y = 0;
		V_L_X = 0; 
		V_R_Y = 0; 
		V_R_X = 0; 
		lines = 0;
		change_koordin = 1;
		pole = 0;
		sreda = 0;
		f = 0;
		soty = 0;
		py = 1;
		d_3 = 1;
		scatter = 1;
		cir = 1;

		if(cir == 1){
			kol_cir = 15;
			kol_cir_num = 40;
			x_num = (kol_cir - 1) * kol_cir_num + 1;
			y_num = 1;
		}
	}

	//разрыв на плоскости
	if (test == 4){
		P_L = 1; 
		P_R = 0.1; 
		D_L = 1; 
		D_R = 0.125;  
		V_L_Y = 0;
		V_L_X = 0; 
		V_R_Y = 0; 
		V_R_X = 0; 
		x_num = 10;
		y_num = 10;
		x_num_pole = 10;
		y_num_pole = 10;
		lines = 0;
		change_koordin = 1;
		pole = 0;
		sreda = 0;
		f = 1;
		soty = 0;
		py = 1;
		d_3 = 1;
		scatter = 1;
		cir = 0;
		y_num_f = int((h * 2) * double(y_num) / y);
		x_num_f = int((h * 2) * double(x_num) / x);
	}

	if (f == 0){
		x_num_f = 0;
		y_num_f = 0;
	}

	if (pole == 0){
		x_num_pole = 0;
		y_num_pole = 0;
	}

	//перенос данных в питон
	x_y_num.open("./out/x_y_num.txt", fstream::out);
	x_y_num << to_string(x_num) << '\n';
	x_y_num << to_string(y_num) << '\n';
	x_y_num << to_string(x_num_f) << '\n';
	x_y_num << to_string(y_num_f) << '\n';
	x_y_num << to_string(scatter) << '\n';
	x_y_num << to_string(lines) << '\n';
	x_y_num << to_string(x_num_pole) << '\n';
	x_y_num << to_string(y_num_pole) << '\n';
	x_y_num << to_string(pole) << '\n';
	x_y_num.close();

	//основные
	double koord_x_mass [x_num * y_num], koord_y_mass [x_num * y_num], velosity_x_mass [x_num * y_num], velosity_y_mass [x_num * y_num];
	double davlenie_mass [x_num * y_num], plotnost_mass [x_num * y_num], energy_mass [x_num * y_num], mass_mass [x_num * y_num], temp_mass [x_num * y_num];

	//фиктивные
	double koord_x_mass_f [(y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num], koord_y_mass_f [(y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num];
	double energy_mass_f [(y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num], mass_mass_f [(y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num];
	double velosity_x_mass_f [(y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num], velosity_y_mass_f [(y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num];
	double davlenie_mass_f [(y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num], plotnost_mass_f [(y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num];
	double temp_mass_f [(y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num];

	//на следущем шаге 
	double koord_x_mass_next [x_num * y_num], koord_y_mass_next [x_num * y_num], velosity_x_mass_next [x_num * y_num], velosity_y_mass_next [x_num * y_num];
	double energy_mass_next [x_num * y_num], davlenie_mass_next [x_num * y_num], plotnost_mass_next [x_num * y_num], tempreture_mass_next [x_num * y_num];

	//поле
	double koord_x_mass_p [x_num_pole * y_num_pole], koord_y_mass_p [x_num_pole * y_num_pole];
	double temp_mass_p [x_num_pole * y_num_pole];
	double energy_mass_p [x_num_pole * y_num_pole];
	double velosity_x_mass_p [x_num_pole * y_num_pole], velosity_y_mass_p [x_num_pole * y_num_pole];
	double davlenie_mass_p [x_num_pole * y_num_pole], plotnost_mass_p [x_num_pole * y_num_pole];

	int osn = 0;
	int fict = 0;

	//определение начальных данных
	if(cir == 1){
		int num = 0;
		for(int i = 0; i < kol_cir; i++){
			if(i == 0){
				koord_x_mass[0] = 0;
				koord_y_mass[0] = 0;
				plotnost_mass [num] = 1 / (pow(i * 0.3,2) + 0.2);
				davlenie_mass [num] = 0.1; 
				mass_mass [num] = 0.1; 
				energy_mass [num] = en(P_L, D_L);
				temp_mass [num] = 0;
				velosity_x_mass [num] = 0;
				velosity_y_mass [num] = 0;
				num += 1;
			}
			else{
				for(int j = 0; j < kol_cir_num; j++){
					if(i % 2 == 0){
						koord_x_mass[num] = (i * 0.5) * cos(j * (360 / kol_cir_num) * PI / 180);
						koord_y_mass[num] = (i * 0.5) * sin(j * (360 / kol_cir_num) * PI / 180);
					}
					else{
						koord_x_mass[num] = (i * 0.5) * cos((j + 0.5) * (360 / kol_cir_num) * PI / 180);
						koord_y_mass[num] = (i * 0.5) * sin((j + 0.5) * (360 / kol_cir_num) * PI / 180);
					}
					plotnost_mass [num] = 1 / (pow(i * 0.5,2) + 0.2);
					davlenie_mass [num] = 0.1;
					mass_mass [num] = 0.1;
					energy_mass [num] = en(P_L, D_L);
					temp_mass [num] = 0;
					velosity_x_mass [num] = 0;
					velosity_y_mass [num] = 0;
					num += 1;
				}
			}
		}
	}
	else if (test == 0){
		//заполнение фиктивных частиц слева
		for(int i = 0; i < y_num + y_num_f * 2; i++){
			for(int j = 0; j < num_l_f; j++){
				plotnost_mass_f [fict] = P_L;
				energy_mass_f [fict] = en(P_L, D_L);
				davlenie_mass_f [fict] = D_L;
				mass_mass_f [fict] = plotnost_mass_f [fict] / (double(num_l_f) * 2 * double(y_num));
				velosity_x_mass_f [fict] = V_L_X;
				velosity_y_mass_f [fict] = V_L_Y;
				temp_mass_f [fict] = 0;
				//temp_mass_f [fict] = tmp (P_L, D_L);
				koord_x_mass_f [fict] = double (j - num_l_f) * (0.5 / num_l);
				koord_y_mass_f [fict] = double (i - y_num_f) / y_num;
				fict++;
			}
		}

		//заполнение фиктивных частиц справа
		for(int i = 0; i < y_num + y_num_f * 2; i++){
			for(int j = 0; j < num_r_f; j++){
				plotnost_mass_f [fict] = P_R;
				energy_mass_f [fict] = en(P_R, D_R);
				davlenie_mass_f [fict] = D_R;
				mass_mass_f [fict] = plotnost_mass_f [fict] / (double(num_r_f) * 2 * double(y_num));
				velosity_x_mass_f [fict] = V_R_X;
				velosity_y_mass_f [fict] = V_R_Y;
				temp_mass_f [fict] = 0;
				//temp_mass_f [fict] = tmp (P_R, D_R);
				koord_x_mass_f [fict] = 1 + double (j + 1) * (0.5 / num_r);
				koord_y_mass_f [fict] = double (i - y_num_f) / double(y_num);
				fict++;
			}
		}

		//заполение частиц x > 0 && x < 0.5
		for(int i = 0; i < y_num + y_num_f * 2; i++){
			for(int j = 0; j < num_l; j++){
				double tmp_y = double (i - y_num_f) / double(y_num);
				double tmp_x = double (j) * (0.5 / num_l);
				if (tmp_y < 0 || tmp_y >= 1){
					plotnost_mass_f [fict] = P_L;
					energy_mass_f [fict] = en(P_L, D_L);
					davlenie_mass_f [fict] = D_L;
					mass_mass_f [fict] = plotnost_mass_f [fict] / (double(num_l_f) * 2 * double(y_num));
					velosity_x_mass_f [fict] = V_L_X;
					velosity_y_mass_f [fict] = V_L_Y;
					temp_mass_f [fict] = 0;
					//temp_mass_f [fict] = tmp (P_L, D_L);
					koord_x_mass_f [fict] = tmp_x;
					koord_y_mass_f [fict] = tmp_y;
					fict++;
				}
				else{
					plotnost_mass [osn] = P_L;
					energy_mass [osn] = en(P_L, D_L);
					davlenie_mass [osn] = D_L;
					mass_mass [osn] = plotnost_mass [osn] / (double(num_l_f) * 2 * double(y_num));
					velosity_x_mass [osn] = V_L_X;
					velosity_y_mass [osn] = V_L_Y;
					temp_mass [osn] = 0;
					//temp_mass [osn] = tmp (P_L, D_L);
					koord_x_mass [osn] = tmp_x;
					koord_y_mass [osn] = tmp_y;
					osn++;
				}
			}
		}
		
		//заполение частиц x > 0.5 && x < 1
		for(int i = 0; i < y_num + y_num_f * 2; i++){
			for(int j = 0; j < num_r; j++){
				double tmp_y = double (i - y_num_f) / double(y_num);
				double tmp_x = 0.5 + double (j + 1) * (0.5 / num_r);
				if (tmp_y < 0 || tmp_y >= 1){
					plotnost_mass_f [fict] = P_R;
					energy_mass_f [fict] = en(P_R, D_R);
					davlenie_mass_f [fict] = D_R;
					mass_mass_f [fict] = plotnost_mass_f [fict] / (double(num_r_f) * 2 * double(y_num));
					velosity_x_mass_f [fict] = V_R_X;
					velosity_y_mass_f [fict] = V_R_Y;
					temp_mass_f [fict] = 0;
					//temp_mass_f [fict] = tmp (P_R, D_R);
					koord_x_mass_f [fict] = tmp_x;
					koord_y_mass_f [fict] = tmp_y;
					fict++;
				}
				else{
					plotnost_mass [osn] = P_R;
					energy_mass [osn] = en(P_R, D_R);
					davlenie_mass [osn] = D_R;
					mass_mass [osn] = plotnost_mass [osn] / (double(num_r_f) * 2 * double(y_num));
					velosity_x_mass [osn] = V_R_X;
					velosity_y_mass [osn] = V_R_Y;
					temp_mass [osn] = 0;
					//temp_mass [osn] = tmp (P_R, D_R);
					koord_x_mass [osn] = tmp_x;
					koord_y_mass [osn] = tmp_y;
					osn++;
				}
			}
		}
	}
	else{
		for(int i = 0; i < y_num + y_num_f * 2; i++){
			for(int j = 0; j < x_num + x_num_f; j++){
				if((((j - x_num_f) / double (x_num) < 0) || ((j - x_num_f) / double (x_num) >= 1)) ||
			       (((i - y_num_f) / double (y_num) < 0) || ((i - y_num_f) / double (y_num) >= 1))){
				   	if(j < ((x_num + x_num_f) / 2)){
						plotnost_mass_f [fict] = P_L;
						energy_mass_f [fict] = en(P_L, D_L);
						davlenie_mass_f [fict] = D_L;
						mass_mass_f [fict] = plotnost_mass_f [fict] / (double(x_num) * double(y_num));
						velosity_x_mass_f [fict] = V_L_X;
						velosity_y_mass_f [fict] = V_L_Y;
						temp_mass_f [fict] = tmp (P_L, D_L);
					}
					else{
						plotnost_mass_f [fict] = P_R;
						energy_mass_f [fict] = en(P_R, D_R);
						davlenie_mass_f [fict] = D_R;
						mass_mass_f [fict] = plotnost_mass_f [fict] / (double(x_num) * double(y_num));
						velosity_x_mass_f [fict] = V_R_X;
						velosity_y_mass_f [fict] = V_R_Y;
						temp_mass_f [fict] = tmp (P_R, D_R);
					}
					if(soty == 1){
						if(i % 2 == 0){
							koord_x_mass_f [fict] = (j + 0.5 - x_num_f) * x / double(x_num);
							koord_y_mass_f [fict] = (i - y_num_f) * y / double(y_num);
						}
						else{
							koord_x_mass_f [fict] = (j - x_num_f) * x / double(x_num);
							koord_y_mass_f [fict] = (i - y_num_f) * y / double(y_num);
						}
					}
					else{
						koord_x_mass_f [fict] = (j - x_num_f) * x / double(x_num);
						koord_y_mass_f [fict] = (i - y_num_f) * y / double(y_num);
					}
					fict++;
				}
				else{
					if(j < ((x_num + x_num_f) / 2)){
						plotnost_mass [osn] = P_L;
						energy_mass [osn] = en(P_L, D_L);
						davlenie_mass [osn] = D_L;
						mass_mass [osn] = plotnost_mass [osn] / (double(x_num) * double(y_num));
						velosity_x_mass [osn] = V_L_X;
						velosity_y_mass [osn] = V_L_Y;
						temp_mass [osn] = tmp (P_L, D_L);
					}
					else{
						plotnost_mass [osn] = P_R;
						energy_mass [osn] = en(P_R, D_R);
						davlenie_mass [osn] = D_R;
						mass_mass [osn] = plotnost_mass [osn] / (double(x_num) * double(y_num));
						velosity_x_mass [osn] = V_R_X;
						velosity_y_mass [osn] = V_R_Y;
						temp_mass [osn] = tmp (P_R, D_R);
					}
					koord_x_mass [osn] = (j - x_num_f) * x / double(x_num);
					koord_y_mass [osn] = (i - y_num_f) * y / double(y_num);
					osn++;
				}
			}
		}
	}

	//определение координат поля
	for(int i = 0; i < y_num_pole; i++){
		for(int j = 0; j < x_num_pole; j++){
			koord_x_mass_p [i * y_num_pole + j] = j / double(x_num_pole);
			koord_y_mass_p [i * y_num_pole + j] = i / double(y_num_pole);
			plotnost_mass_p [i * y_num_pole + j] = 0;
			temp_mass_p [i * y_num_pole + j] = 0;
			velosity_x_mass_p [i * y_num_pole + j] = 0;
			velosity_y_mass_p [i * y_num_pole + j] = 0;
			energy_mass_p [i * y_num_pole + j] = 0;
			davlenie_mass_p [i * y_num_pole + j] = 0;
		}
	}

	//список для копирования в фиктивные
	double ind_mirror [y_num_f * 2 * (num_r + num_l)][2];
	int reek = 0;

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
					temp_mass_p [i * y_num_pole + j] = 0;
				}
			}
			for(int i = 0; i < y_num_pole; i++){
				for(int j = 0; j < x_num_pole; j++){
					for(int a = 0; a < y_num * x_num; a++){
						if(rast(koord_x_mass [a], koord_x_mass_p [i * y_num_pole + j], koord_y_mass [a], koord_y_mass_p [i * y_num_pole + j]) < (h * 2)){
							double r = rast(koord_x_mass [a], koord_x_mass_p [i * y_num_pole + j], koord_y_mass [a], koord_y_mass_p [i * y_num_pole + j]);
							plotnost_mass_p [i * y_num_pole + j] += mass_mass [a] * W(r, h);
							velosity_x_mass_p [i * y_num_pole + j] += mass_mass [a] * (velosity_x_mass [a] / plotnost_mass [a]) * W(r, h);
							velosity_y_mass_p [i * y_num_pole + j] += mass_mass [a] * (velosity_y_mass [a] / plotnost_mass [a]) * W(r, h);
							energy_mass_p [i * y_num_pole + j] += (mass_mass [a] / plotnost_mass [a]) * energy_mass [a] * W(r, h);
							temp_mass_p [i * y_num_pole + j] += (mass_mass [a] / plotnost_mass [a]) * temp_mass [a] * W(r, h);
						}
						davlenie_mass_p [i * y_num_pole + j] = davlenie (plotnost_mass_p [i * y_num_pole + j], energy_mass_p [i * y_num_pole + j]);
					}
					for(int a = 0; a < (y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num; a++){
						if(rast(koord_x_mass_f [a], koord_x_mass_p [i * y_num_pole + j], koord_y_mass_f [a], koord_y_mass_p [i * y_num_pole + j]) < (h * 2)){
							double r = rast(koord_x_mass_f [a], koord_x_mass_p [i * y_num_pole + j], koord_y_mass_f [a], koord_y_mass_p [i * y_num_pole + j]);
							plotnost_mass_p [i * y_num_pole + j] += mass_mass_f [a] * W(r, h);
							velosity_x_mass_p [i * y_num_pole + j] += mass_mass_f [a] * (velosity_x_mass_f [a] / plotnost_mass_f [a]) * W(r, h);
							velosity_y_mass_p [i * y_num_pole + j] += mass_mass_f [a] * (velosity_y_mass_f [a] / plotnost_mass_f [a]) * W(r, h);
							energy_mass_p [i * y_num_pole + j] += (mass_mass_f [a] / plotnost_mass_f [a]) * energy_mass_f [a] * W(r, h);
							temp_mass_p [i * y_num_pole + j] += (mass_mass_f [a] / plotnost_mass_f [a]) * temp_mass_f [a] * W(r, h);
						}
						davlenie_mass_p [i * y_num_pole + j] = davlenie (plotnost_mass_p [i * y_num_pole + j], energy_mass_p [i * y_num_pole + j]);
					}
				}
			}
		}
		if(pole == 1 && k == 100){
			files_pole("k_x_", k, x_num_pole, y_num_pole, koord_x_mass_p);
			files_pole("k_y_", k, x_num_pole, y_num_pole, koord_y_mass_p);
			files_pole("e_", k, x_num_pole, y_num_pole, energy_mass_p);
			files_pole("v_x_", k, x_num_pole, y_num_pole, velosity_x_mass_p);
			files_pole("v_y_", k, x_num_pole, y_num_pole, velosity_y_mass_p);
			files_pole("d_", k, x_num_pole, y_num_pole, davlenie_mass_p);
			files_pole("p_", k, x_num_pole, y_num_pole, plotnost_mass_p);
			files_pole("t_", k, x_num_pole, y_num_pole, temp_mass_p);
		}
		/*else{
			files("k_x_", k, x_num, y_num, x_num_f, y_num_f, koord_x_mass, koord_x_mass_f);
			files("k_y_", k, x_num, y_num, x_num_f, y_num_f, koord_y_mass, koord_y_mass_f);
			files("e_", k, x_num, y_num, x_num_f, y_num_f, energy_mass, energy_mass_f);
			files("v_x_", k, x_num, y_num, x_num_f, y_num_f, velosity_x_mass, velosity_x_mass_f);
			files("v_y_", k, x_num, y_num, x_num_f, y_num_f, velosity_y_mass, velosity_y_mass_f);
			files("d_", k, x_num, y_num, x_num_f, y_num_f, davlenie_mass, davlenie_mass_f);
			files("p_", k, x_num, y_num, x_num_f, y_num_f, plotnost_mass, plotnost_mass_f);
			files("t_", k, x_num, y_num, x_num_f, y_num_f, temp_mass, temp_mass_f);
		}*/
		

		//вычисление по шагам
		for(int i = 0; i < y_num * x_num; i++){
			plotnost_mass_next [i] = plotnost (mass_mass, plotnost_mass, velosity_x_mass, velosity_y_mass, koord_x_mass, koord_y_mass, velosity_x_mass_f, velosity_y_mass_f, koord_x_mass_f, koord_y_mass_f, mass_mass_f, plotnost_mass_f, i, t_step, h, y_num, x_num, x_num_f, y_num_f, sreda);
			energy_mass_next [i] = energy (energy_mass, davlenie_mass, davlenie_mass_f, mass_mass, plotnost_mass, velosity_x_mass, velosity_y_mass, koord_x_mass, koord_y_mass,velosity_x_mass_f, velosity_y_mass_f, koord_x_mass_f, koord_y_mass_f, mass_mass_f, plotnost_mass_f, i, t_step, h, y_num, x_num, x_num_f, y_num_f, sreda);
			velosity_x_mass_next [i] = skotost_x (davlenie_mass, mass_mass, plotnost_mass, velosity_x_mass, velosity_y_mass, koord_x_mass, koord_y_mass, davlenie_mass_f, velosity_x_mass_f, velosity_y_mass_f, koord_x_mass_f, koord_y_mass_f, mass_mass_f, plotnost_mass_f, i, t_step, h, y_num, x_num, x_num_f, y_num_f, sreda);
			//velosity_y_mass_next [i] = skotost_y (davlenie_mass, mass_mass, plotnost_mass, velosity_x_mass, velosity_y_mass, koord_x_mass, koord_y_mass, davlenie_mass_f, velosity_x_mass_f, velosity_y_mass_f, koord_x_mass_f, koord_y_mass_f, mass_mass_f, plotnost_mass_f, i, t_step, h, y_num, x_num, x_num_f, y_num_f, sreda);
			davlenie_mass_next [i] = davlenie (plotnost_mass [i], energy_mass [i]);
			//tempreture_mass_next [i] = tempreture (mass_mass, plotnost_mass, koord_x_mass, koord_y_mass, koord_x_mass_f, koord_y_mass_f, mass_mass_f, plotnost_mass_f, temp_mass, temp_mass_f, davlenie_mass_f, davlenie_mass, i, t_step, h, y_num, x_num, x_num_f, y_num_f);
			tempreture_mass_next [i] = 0;
			velosity_y_mass_next [i] = 0;

		}

		//смена координат
		if(change_koordin){
			for(int a = 0; a < y_num * x_num; a++){
					koord_x_mass_next [a] = velosity_x_mass [a] * t_step + koord_x_mass [a];
					koord_y_mass_next [a] = velosity_y_mass [a] * t_step + koord_y_mass [a];
			}
			for(int a = 0; a < (y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num; a++){
				koord_x_mass_f [a] = velosity_x_mass_f [a] * t_step + koord_x_mass_f [a];
				koord_y_mass_f [a] = velosity_y_mass_f [a] * t_step + koord_y_mass_f [a];
			}
		}

		
		//фиктивные сверху и снизу копируют
		if (test == 0){
			if(k == 0){
				for(int a = 0; a < (y_num + y_num_f * 2) * (x_num + x_num_f) - x_num * y_num; a++){
					if ((koord_y_mass_f [a] < 0 || koord_y_mass_f [a] >= 1) && (koord_x_mass_f [a] >= 0 && koord_x_mass_f [a] < 1)){
						for(int j = 0; j < x_num * y_num; j++){
							if(koord_x_mass [j] == koord_x_mass_f [a]){
								ind_mirror [reek][0] = a;
								ind_mirror [reek][1] = j;
								reek++;
								break;
							}
						}
					}
				}
				//cout << reek << endl;
			}
			for(int a = 0; a < reek; a++){
				int index_1 = ind_mirror [a][0];
				int index_2 = ind_mirror [a][1];
				plotnost_mass_f [index_1] = plotnost_mass_next [index_2];
				energy_mass_f [index_1] = energy_mass_next [index_2];
				davlenie_mass_f [index_1] = davlenie_mass_next [index_2];
				velosity_x_mass_f [index_1] = velosity_x_mass_next [index_2];
				velosity_y_mass_f [index_1] = velosity_y_mass_next [index_2];
				temp_mass_f [index_1] = tempreture_mass_next [index_2];
			}
		}

		//перенос в основные массивы
		for(int i = 0; i < y_num * x_num; i++){
			plotnost_mass [i] = plotnost_mass_next [i];
			energy_mass [i] = energy_mass_next [i];
			davlenie_mass [i] = davlenie_mass_next [i];
			velosity_x_mass [i] = velosity_x_mass_next [i];
			velosity_y_mass [i] = velosity_y_mass_next [i];
			temp_mass [i] = tempreture_mass_next [i];
			if(change_koordin){ 
				koord_x_mass [i] = koord_x_mass_next [i];
				koord_y_mass [i] = koord_y_mass_next [i];
			}
		}
	}

	if (py){
		if(d_3 == 1){
			system("python graph3d.py");
		}
		else{
			system("python graph.py");
		}
	}
 
	return 0;
}
