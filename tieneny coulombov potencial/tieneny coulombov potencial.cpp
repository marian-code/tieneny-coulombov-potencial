// tieneny coulombov potencial.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <string>
#include <sstream>
#include <windows.h>

using namespace std;

const double pi = 3.14159265358979323;

double function_f1(double x, double y1, double y2)
{
	return y2;
}

double function_f2(double x, double y1, double y2, const double k_s, const double energy)
{
	return ((-2.0 / (x))*exp(-k_s*x) - 2 * energy)*y1;
}

vector<double> modified_midpoint(int i, double x, int n_steps, double h_small, double H_big, double y1, double y2, const double k_s, const double energy)
{
	int m;
	vector<double> z1, z2;
	vector<double> y;

	//priradit hodnoty na kraji intervalu y(i) -> z(0)
	z1.push_back(y1);
	z2.push_back(y2);

	//cout << "z1(" << 0 << ") " << z1[0] << "\t" << "z2(" << 0 << ") " << z2[0] << endl;

	z1.push_back(z1[0] + h_small*function_f1(x, z1[0], z2[0]));
	z2.push_back(z2[0] + h_small*function_f2(x, z1[0], z2[0], k_s, energy));

	//cout << "z1(" << 1 << ") " << z1[1] << "\t" << "z2(" << 1 << ") " << z2[1] << endl;

	//itegrovat so zjemnenym krokom na intervale od z(1) -> z(n)
	for (m = 1; m < n_steps; m++)
	{
		z1.push_back(z1[m - 1] + 2 * h_small*function_f1(x * m, z1[m], z2[m]));
		z2.push_back(z2[m - 1] + 2 * h_small*function_f2(x * m, z1[m], z2[m], k_s, energy));

		//cout << "z1(" << m + 1 << ") " << z1[m + 1] << "\t" << "z2(" << m + 1 << ") " << z2[m + 1] << endl;
	}

	//cout << "z1(" << m + 1 << ") " << (z1[n_steps] + z1[n_steps - 1] + h_small*function_f1(x + H_big, z1[n_steps], z2[n_steps])) / 2 << "\t" 
	//	 << "z2(" << m + 1 << ") " << (z2[n_steps] + z2[n_steps - 1] + h_small*function_f2(x + H_big, z1[n_steps], z2[n_steps], k_s, energy)) / 2 << endl;

	//spocitat dalsi krok -> y(i + H) pomocou z(n) a z(n - 1)
	y.push_back((z1[n_steps] + z1[n_steps - 1] + h_small*function_f1(x + H_big, z1[n_steps], z2[n_steps])) / 2);
	y.push_back((z2[n_steps] + z2[n_steps - 1] + h_small*function_f2(x + H_big, z1[n_steps], z2[n_steps], k_s, energy)) / 2);

	return y;
}

vector<vector<double>> Gragg_Bulirsch_Stoer(const double epsilon, const unsigned int max_n_steps, double H_max, const double H_min, const double x_initial, const double v_initial,
											const bool x_write, const bool H_write, const unsigned int output_points, const double k_s, const double energy, const double cut_off)
{
	bool repeat = true;
	bool reduce_step = false;
	unsigned int i, j, k, m, counter = 1;
	unsigned int n_steps;
	unsigned int i_max, increment;
	unsigned int divergence;
	double H_big, h_small;
	double local_error;
	double s;
	vector<double> scale;
	vector<double> H_evolution;
	vector<double> x, temp;
	vector<double> delta;
	vector<double> y1, y2;
	vector<vector<double>> Richardson_1, Richardson_2;
	vector<vector<double>> output;

	const double absolute_tolerance = epsilon;
	const double relative_tolerance = epsilon;

	output.resize(3);

	Richardson_1.resize(max_n_steps);
	for (i = 0; i < max_n_steps; i++)
		Richardson_1[i].resize(i + 1);

	Richardson_2.resize(max_n_steps);
	for (i = 0; i < max_n_steps; i++)
		Richardson_2[i].resize(i + 1);

	delta.resize(2);
	scale.resize(2);

	cout.precision(15);
	cout.setf(std::ios::fixed, std::ios::floatfield);

	H_big = H_min;//velkost kroku

	i = 0;
	x.push_back(H_min/1000);//okrajova podmienka

	y1.push_back(x_initial);//pocitocna podmienka x(0)
	y2.push_back(v_initial);//pocitocna podmienka x´(0)

	cout << "running Gragg Bulirsch Stoer... " << endl;

	while (x[i] < cut_off)
	{
		//cout << "----------------------- step " << counter << "-----------------------------------" << endl;
		//cout << "space step: " << i << endl;
		//cout << "stepsize H: " << H_big << endl << endl;

		if (x[i] >= cut_off)
		{
			cout << "pri zadanom cutoffe a energii nebol najdeny viazany stav" << endl;
			break;
		}

		/*if (reduce_step == false)
		{
			if (y1[i] < 1E-10)
			{
				cout << "reducing stepsize" << endl;
				cout << "x: " << x[i] << endl;
				cout << "stepsize H: " << H_big << endl;

				x.pop_back();
				y1.pop_back();
				y2.pop_back();
				H_max /= 10000;
				H_big = H_max;
				reduce_step = true;
				i--;

				cout << "after delete: " << endl;
				cout << "x: " << x[i] << "\t" << "y1: " << y1[i] << endl;
				cout << "stepsize H: " << H_big << endl << endl;

				continue;
			}
		}*/

		//podmienka na ukoncenie cyklu
		if (x[i] > 1)
			if (y1[i] < 1E-10) break;

		//max_n_steps urcuje max pocet podintervalov na intervale H
		for (k = 0; k < max_n_steps; k++)
		{
			n_steps = 2 * (k + 1); //pocet podintevalov v aktualnom kroku extrapolacie
			h_small = H_big / n_steps;

			//spocitat dalsi krok pomocou modified midpoint rule -> y(i + H)
			temp = modified_midpoint(i, x[i], n_steps, h_small, H_big, y1[i], y2[i], k_s, energy);//ak by som priradoval postupne tak by sa funkcia musela vykonat 2x
			Richardson_1[k][0] = temp[0];
			Richardson_2[k][0] = temp[1];

			//ak niesom v prvom kroku tak pomocou Richardsona dopocitat aproximacie vyssich radov
			if (k > 0)
			{
				for (j = 0; j < k; j++)
				{
					Richardson_1[k][j + 1] = Richardson_1[k][j] + (Richardson_1[k][j] - Richardson_1[k - 1][j]) / (pow(double(n_steps) / (2 * (k - j)), 2) - 1);
					Richardson_2[k][j + 1] = Richardson_2[k][j] + (Richardson_2[k][j] - Richardson_2[k - 1][j]) / (pow(double(n_steps) / (2 * (k - j)), 2) - 1);
				}

#pragma region kontrolne vypisy

				//kontrolny vypis
				/*cout << "Richardson 1" << endl;
				for (j = 0; j < k+1; j++)
				{
					for (m = 0; m < j + 1; m++)
					{
						if (Richardson_1[j][m] != 0)
							cout << Richardson_1[j][m] << "\t";
					}
					cout << endl;
				}
				cout << endl;
				cout << "Richardson 2" << endl;
				for (j = 0; j < k+1; j++)
				{
					for (m = 0; m < j + 1; m++)
					{
						if (Richardson_2[j][m] != 0)
							cout << Richardson_2[j][m] << "\t";
					}
					cout << endl;
				}
				cout << endl;*/

#pragma endregion

				//spocitat rozdiely aproximacii
				delta[0] = abs(Richardson_1[k][k] - Richardson_1[k][k - 1]);
				delta[1] = abs(Richardson_2[k][k] - Richardson_2[k][k - 1]);

				//num recepies odhad chyby
				/*scale[0] = absolute_tolerance + abs(Richardson_1[k][k])*relative_tolerance;
				scale[1] = absolute_tolerance + abs(Richardson_2[k][k])*relative_tolerance;

				local_error = 0;
				for (m = 0; m < 2; m++)
				{
					local_error += pow(delta[m] / scale[m], 2);

					cout << "delta/scale = " << delta[m] << " / " << scale[m] << " --> " << delta[m] / scale[m] << endl;
				}

				local_error = sqrt(local_error / 2);*/

				sort(delta.begin(), delta.end());
				local_error = delta[1];
				//cout << "lokalna chyba: " << local_error << endl;
				//cout << "delta: " << delta[0] << "\t" << delta[1] << endl << endl;

				if (local_error <epsilon)
				{
					repeat = false;
					//cout << "**********************************  break  *************************************" << endl;
					break;
				}

				//num recepies odhad chyby
				////ked sa dosiahne pozadovana lokalna presnost -> break
				//if (local_error <= 1)
				//{
				//	repeat = false;
				//	cout << "**********************************  break  *************************************" << endl;
				//	break;
				//}
			}
		}
		k--; //k treba o 1 zmensit lebo po skonceni cyklu sa este pricita 1 naviac

		// repeat = true ked sa v danom kroku nepodarilo pomocou richardsonovej extrapolacie
		//skonvergovat na pozadovanu presnost a je potrebne opakovat krok s inym H
		if (repeat == false)
		{
			//priradit novy krok do trajektorie
			y1.push_back(Richardson_1[k][k]);
			y2.push_back(Richardson_2[k][k]);

			//posunut sa do dalsieho casoveho kroku
			x.push_back(x[i] + H_big);

			i++;

			if (i % 1000 == 0)
			{
				cout << x[i] << "\t" << y1[i] << "\t" << H_big << endl;
			}
		}

		//spocitat velkost noveho kroku H - adaptive stepsize control
		//s = 0.98*pow(0.98 / local_error, 1.0 / (2 * k + 1));
		s = 0.9*pow((10E6*epsilon) / (10E6*local_error), 1.0 / (2 * k + 1));

		// empiricke obmedzenia, prevzate z num.recepies, h by sa nemalo menit velmi rychlo
		if (s > 10) s = 10; //h sa nikdy nezväèší o faktor viac ako 10
		else if (s < 0.2) s = 0.2; //h sa nikdy nezmenší o faktor viac ako 5

		H_big = H_big*s;

		if (H_big > H_max) H_big = H_max;
		if (H_big < H_min) H_big = H_min;

		if (repeat == false) H_evolution.push_back(H_big);
		repeat = true;

		if (y1[i] > y1[i - 1] && s < 1) divergence++;
		if (divergence > 10000) break;

		//cout << "<<<<<<<<<<<<<<<<<<<< repeat: " << (repeat ? "true" : "false") <<" >>>>>>>>>>>>>>>>>>"<< endl;
		//cout << "----------------------- end of step " << counter << "------------------------------" << endl << endl << endl;
		counter++;
	}

	//odstranit zaporne hodnoty ak sa nejake nachadzaju na konci vlnovej funkcie
	/*while (y1[y1.size() - 1] < 0)
	{
		y1.pop_back();
		x.pop_back();
	}*/

	cout << "Job done " << endl << endl;

	//ulozenie dat pre vystup s funkcie
	output.resize(1);
	output[0].resize(1);
	output[0][0] = i;
	output.push_back(x);
	output.push_back(y1);

#pragma region zápis na disk

	increment = 1;

	if (y1.size() > output_points)
	{
		i_max = output_points;
		while (y1.size() % i_max != 0)
		{
			i_max++;
		}
		increment = y1.size() / i_max;
	}

	if (x_write == true)
	{
		ostringstream fileNameStream1("subor1");
		fileNameStream1 << "../output/k_s_" << k_s << "/Phi_" << energy << ".dat";
		string fileName1 = fileNameStream1.str();

		ofstream outfile1; //zapis do suboru
		outfile1.precision(13);
		outfile1.setf(std::ios::fixed, std::ios::floatfield);

		outfile1.open(fileName1);
		for (i = 0; i < y1.size(); i+=increment)
		{
			outfile1 << x[i] << "\t";
			outfile1 << y1[i] << endl;
		}
		outfile1.close();
	}

	if (H_write == true)
	{
		ostringstream fileNameStream2("subor2");
		fileNameStream2 << "../output/k_s_" << k_s << "/H_" << energy << ".dat";
		string fileName2 = fileNameStream2.str();

		ofstream outfile2; //zapis do suboru
		outfile2.precision(13);
		outfile2.setf(std::ios::fixed, std::ios::floatfield);

		outfile2.open(fileName2);
		for (i = 0; i < y1.size() - 1; i += increment)
		{

			outfile2 << x[i] << "\t";
			outfile2 << H_evolution[i] << endl;
		}
		outfile2.close();
	}

	cout << "trajectories succesfully writen to disc!" << endl << endl;

#pragma endregion

	return output;
}

int main()
{
	unsigned int i, j;
	double energy, energy_low, energy_high;
	double y1_old;
	vector<double> x;
	vector<double> y1;
	vector<vector<double>> Richardson_1, Richardson_2;
	vector<vector<double>> input;

	const bool x_write = true;
	const bool H_write = true;
	const unsigned int output_points = 1000;
	const unsigned int max_n_steps = 10;
	const double epsilon = 1E-4;
	const double x_initial = 0;//pociatocna podmienka 
	const double v_initial = 1;//pociatocna podmienka 1. derivacia
	const double H_max = 0.5;
	const double H_min = 0.0000005;
	const double cut_off = pow(10,2);


	double k_s = 0.0;

	cout.precision(13);
	cout.setf(std::ios::fixed, std::ios::floatfield);

	//vytvorit priecinok pre vystup
	system("mkdir ..\\output");

	auto time_start = chrono::high_resolution_clock::now();

	//vytvorit podpriecinok pre prislusne k_s
	ostringstream fileNameStream1("subor1");
	fileNameStream1 << "mkdir ..\\output\\k_s_" << k_s;
	string fileName1 = fileNameStream1.str();
	const char * c = fileName1.c_str();
	system(c);
	cout << endl;

	energy_low = -0.8;
	energy_high = 0.0;
	y1_old = 10;

	for (j = 0; j < 100; j++)
	{
		energy = (energy_high + energy_low) / 2;
		//energy = -0.50015;
		cout << "----------------------- bisection cycle: " << j + 1 << " ------------------------" << endl;
		cout << "energy is: " << energy << endl << endl;

		input.swap(Gragg_Bulirsch_Stoer(epsilon, max_n_steps, H_max, H_min, x_initial, v_initial, x_write, H_write, output_points, k_s, energy, cut_off));
		i = input[0][0];
		x.swap(input[1]);
		y1.swap(input[2]);

		cout << "phi at cutoff -->  new:" << y1[y1.size() - 1] << " old: " << y1_old << endl;

		if (y1[y1.size() - 1] > y1_old)
		{
			y1_old = y1[y1.size() - 1];
			energy_low = energy;
		}
		else
		{
			y1_old = y1[y1.size() - 1];
			energy_high = energy;
		}

		cout << "new energy interval: " << "( " << energy_low << " , " << energy_high << " )" << endl;

		

		if (abs(y1[y1.size() - 1]) < 1E-5)
		{
			cout << "convergence archived!!!" << endl;
			break;
		}
		if (abs(energy_high - energy_low) < 1E-12)
		{
			cout << "convergence archived!!!" << endl;
			break;
		}

		cout << "-------------------------- end of cycle: " << j + 1 << " ------------------------" << endl << endl;
	}


	auto time_end = chrono::high_resolution_clock::now();
	cout << endl << "total CPU time: " << chrono::duration_cast<chrono::nanoseconds>(time_end - time_start).count()*10E-10 << "s" << endl << endl;


	system("PAUSE");

	return 0;
}



