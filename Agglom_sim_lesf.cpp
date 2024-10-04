#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

#include "settings_lesf.h"

using namespace std;

// half interval starting from 1/2
vector<double> F;
// full interval starting from 0
vector<double> s;
// full interval starting from 0
vector<double> Pc;

vector<double> C_cu;
vector<double> F_cu;

vector<double> C_grain;
vector<double> C_precip;

double dr, t;

#ifndef dt
double dt;
#endif

void setup_value(void)
{
	F.resize(r_max);
	s.resize(r_max);
	Pc.resize(r_max);

	C_cu.resize(r_max);
	F_cu.resize(r_max);

	C_grain.resize(r_max);
	C_precip.resize(r_max);

	for (int i = 0; i < r_max; i++)
	{
		s[i] = sat_0;
		// C_cu[i] = 1.0;
		C_cu[i] = 0.0;

		C_grain[i] = C_grain_0;
		C_precip[i] = 0.0;
	}

	dr = radius / (r_max - 1.0);

#ifndef dt
	double dt_water = cfl_num * dr * dr / (k_max * c_cap / mu);
	double dt_dif = cfl_num * dr * dr / (2.0 * D);
	double dt_evap = cfl_num * dr / v_evap;

	dt = min(dt_water, dt_dif);
	dt = min(dt, dt_evap);
#endif

	cout << "dr = " << dr << "\tdt = " << dt << endl;
}

double P(double s)
{
	if (s < s_min)
		s = s_min * 1.0000001;

	double s_w = (s - s_min) / (1.0 - s_min);

	/*van Geruchten equation*/
	return c_cap * pow(pow(s_w, -n_pc / (1.0 - n_pc)), 1. / n_pc);

	// return -c_cap * pow(s_w, -1.0 / 3.0);
}

double k(double s)
{
	if (s < s_min)
		s = s_min * 1.0000001;

	double s_w = (s - s_min) / (1.0 - s_min);

	/*Derivated from Brooks-Corey model (?)*/
	return k_max * pow(s_w, n_k);
}

double m_react(double C_cu, double C_grain, double s)
{
	if (C_cu < C_sat)
		return s * porosity * Beta * k_react * pow(C_grain / C_grain_0, 2.0 / 3.0) * (C_sat - C_cu);
	else
		return 0.0;
}

double m_precip(double C_cu, double s)
{
	if (C_cu > C_sat)
		return s * porosity * k_precip * (C_cu - C_sat);
	else
		return 0;
}

void calc_flux(void)
{
	for (int i = 0; i < r_max; i++)
		Pc[i] = P(s[i]);

	for (int i = 0; i < r_max - 1; i++)
	{
		double s_ave = (s[i] + s[i + 1]) / 2.0;
		F[i] = -(s[i] * porosity * k(s_ave) / mu) * (Pc[i + 1] - Pc[i]) / dr;

		double C_av = (C_cu[i] + C_cu[i + 1]) / 2.0;

		F_cu[i] = F[i] * C_av - s[i] * porosity * D / tort * (C_cu[i + 1] - C_cu[i]) / dr;

		// cout << -(s[i] * porosity * k(s_ave) / mu) * (Pc[i + 1] - Pc[i]) / dr <<endl;
	}

	F[r_max - 1] = std::max(s[r_max - 1] - s_min, 0.0) * porosity * v_evap;
	F_cu[r_max - 1] = 0.0;

	// cout << F[r_max - 1] << " : " << (1.0 / (radius * radius * porosity)) * ( F[r_max - 1] * radius * radius) / dr << endl;
}

void update_values(void)
{
	for (int i = 1; i < r_max; i++)
	{
		double r_ave = dr * i, r_down = r_ave - 0.5 * dr, r_up = r_ave + 0.5 * dr;

		double ds_dt = -(1.0 / (r_ave * r_ave * porosity)) * (F[i] * r_up * r_up - F[i - 1] * r_down * r_down) / dr;
		double dC_cu_dt = -(1.0 / (r_ave * r_ave * porosity * s[i])) * (F_cu[i] * r_up * r_up - F_cu[i - 1] * r_down * r_down) / dr;

		double dm_grain_dt = m_react(C_cu[i], C_grain[i], s[i]);
		double dm_precip_dt = m_precip(C_cu[i], s[i]);

		dC_cu_dt += ((dm_grain_dt - dm_precip_dt) / (porosity * s[i]) - ds_dt * C_cu[i] / s[i]);

		s[i] += ds_dt * dt;
		C_cu[i] += dC_cu_dt * dt;

		C_grain[i] -= dm_grain_dt * dt;
		C_precip[i] += dm_precip_dt * dt;

		if (C_grain[i] < 0.0)
			C_grain[i] = 0.0;
	}

	s[0] = s[1];
	C_cu[0] = C_cu[1];
	C_grain[0] = C_grain[1];
	C_precip[0] = C_precip[1];
}

int main(void)
{
	double t_write = dt_out;
	fstream f1, f2, f3, f4;

	setup_value();

	t = 0;
	f1.open("sat.txt", ios_base::out);
	f2.open("C_cu.txt", ios_base::out);
	f3.open("C_grain.txt", ios_base::out);
	f4.open("C_precip.txt", ios_base::out);

	f1 << "r\t";
	for (int i = 0; i < r_max - 1; i++)
		f1 << dr * i << "\t";
	f1 << endl;

	f2 << "r\t";
	for (int i = 0; i < r_max - 1; i++)
		f2 << dr * i << "\t";
	f2 << endl;

	f3 << "r\t";
	for (int i = 0; i < r_max - 1; i++)
		f3 << dr * i << "\t";
	f3 << endl;

	f4 << "r\t";
	for (int i = 0; i < r_max - 1; i++)
		f4 << dr * i << "\t";
	f4 << endl;

	while (t < t_max)
	{
		calc_flux();
		update_values();

		if (t > t_write)
		{
			cout << t << "\ts_e: " << s[r_max - 1] << "\ts_m: " << s[(r_max / 2)] << "\tC_cu_e: " << C_cu[r_max - 1] << "\tC_cu_m: " << C_cu[(r_max / 2)] << "\tC_grain_e: " << C_grain[r_max - 1] << "\tC_grain_m: " << C_grain[r_max / 2] << endl;

			f1 << t_write << "\t";
			f2 << t_write << "\t";
			f3 << t_write << "\t";
			f4 << t_write << "\t";
			for (int i = 0; i < r_max; i++)
			{
				f1 << s[i] << "\t";
				f2 << C_cu[i] << "\t";
				f3 << C_grain[i] << "\t";
				f4 << C_precip[i] << "\t";
			}
			f1 << endl;
			f2 << endl;
			f3 << endl;
			f4 << endl;

			t_write += dt_out;
		}

		t = t + dt;
	}

	f1.close();
	return 0;
}