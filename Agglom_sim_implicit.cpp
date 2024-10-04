#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>

#include "settings.h"

using namespace std;

// half interval starting from 1/2
vector<double> F;
vector<double> F_prev;
// full interval starting from 0
vector<double> s;
vector<double> s_prev;
// full interval starting from 0
vector<double> Pc;

vector<double> C_cu;
vector<double> F_cu;

vector<vector<double>> C_grain;
vector<double> C_precip;

vector<double> k_react;

vector<double> a, b, c, f, gam;

double dr, t;

#ifndef dt
double dt;
#endif

void Tri_Diag(double *a, double *b, double *c, double *f, double *u, int n)
{
	int i;
	double bet;

	u[0] = f[0] / c[0];
	bet = c[0];
	for (i = 1; i < n; i++)
	{
		gam[i] = a[i - 1] / bet;
		bet = c[i] - b[i] * gam[i];
		u[i] = (f[i] - b[i] * u[i - 1]) / bet;
	}
	for (i = (n - 2); i >= 0; i--)
		u[i] -= gam[i + 1] * u[i + 1];
}

void setup_value(void)
{
	F.resize(r_max);
	F_prev.resize(r_max);
	s.resize(r_max);
	s_prev.resize(r_max);
	Pc.resize(r_max);

	C_cu.resize(r_max);
	F_cu.resize(r_max);

	C_precip.resize(r_max);

	C_grain.resize(num_grain_classes, vector<double>(r_max));

	a.resize(r_max);
	b.resize(r_max);
	c.resize(r_max);
	f.resize(r_max);

	gam.resize(r_max);

	for (int i = 0; i < r_max; i++)
	{
		s[i] = sat_0;
		// C_cu[i] = 1.0;
		C_cu[i] = 0.0;
		F[i] = 0.0;

		for (int n = 0; n < num_grain_classes; n++)
			C_grain[n][i] = C_grain_0;
		C_precip[i] = 0.0;
	}

	dr = radius / (r_max - 1.0);

#ifndef dt
	double dt_water = cfl_num * dr * dr / (k_max * c_cap / mu);
	double dt_dif = cfl_num * dr * dr / (2.0 * D);
	double dt_evap = cfl_num * dr / v_evap;

	// dt =  min(dt_water, dt_dif);
	dt = min(dt_water, dt_evap);
#endif

	k_react.resize(num_grain_classes);

	for (int n = 0; n < num_grain_classes; n++)
	{
		k_react[n] = //....... function of which calss it is
	}

	//...or
	k_react[0] = ....;
	k_react[1] = ....;

	cout << "dt_water = " << dt_water << " dt_dif = " << dt_dif << " dt_evap = " << dt_evap << endl;
	cout << "dr = " << dr << "\tdt = " << dt << endl;
}

double P(double s)
{
	/*if (s < s_min)
		s = s_min * 1.0000001;*/

	double s_w = (s - s_min) / (1.0 - s_min);

	/*van Geruchten equation*/
	return c_cap * pow(pow(s_w, -n_pc / (1.0 - n_pc)), 1. / n_pc);

	// return -c_cap * pow(s_w, -1.0 / 3.0);
}

double k(double s)
{
	/*if (s < s_min)
		s = s_min * 1.0000001;*/

	double s_w = (s - s_min) / (1.0 - s_min);

	/*Derivated from Brooks-Corey model (?)*/
	return k_max * pow(s_w, n_k);
}

double K_react(double C_cu, double C_grain, double s, int grain_class)
{
	// k_react should be a vector that depends on initial size and surface kinetics

	if (C_cu < C_sat)
		return s * porosity * Beta * k_react[grain_class] * pow(C_grain / C_grain_0, 2.0 / 3.0);
	else
		return 0.0;
}

double K_precip(double C_cu, double s)
{
	if (C_cu > C_sat)
		return s * porosity * k_precip;
	else
		return 0.0;
}

void calc_flux_water(void)
{
	F.swap(F_prev);
	s.swap(s_prev);

	for (int i = 0; i < r_max; i++)
		Pc[i] = P(s[i]);

	for (int i = 0; i < r_max - 1; i++)
	{
		double s_ave = (s[i] + s[i + 1]) / 2.0;
		F[i] = -(s_ave * porosity * k(s_ave) / mu) * (Pc[i + 1] - Pc[i]) / dr;
	}

	F[r_max - 1] = std::max(s[r_max - 1] - s_min, 0.0) * porosity * v_evap;
}

void update_values_saturation(void)
{
	for (int i = 1; i < r_max; i++)
	{
		double r_ave = dr * i, r_down = r_ave - 0.5 * dr, r_up = r_ave + 0.5 * dr;

		double ds_dt = -(1.0 / (r_ave * r_ave * porosity)) * (F[i] * r_up * r_up - F[i - 1] * r_down * r_down) / dr;

		s[i] += ds_dt * dt;
	}

	s[0] = s[1];
}

void setup_conc_matrix(void)
{
	for (int i = 1; i < r_max - 1; i++)
	{
		double r_ave = dr * i, r_down = r_ave - 0.5 * dr, r_up = r_ave + 0.5 * dr;
		double s_up = (s[i] + s[i + 1]) / 2.0, s_down = (s[i] + s[i - 1]) / 2.0;
		double ds_dt = -(1.0 / (r_ave * r_ave * porosity)) * (F[i] * r_up * r_up - F[i - 1] * r_down * r_down) / dr;

		double K_r_sum = 0.0;
		for (int n = 0; n < num_grain_classes; n++)
			K_r_sum += K_react(C_cu[i], C_grain[n][i], s[i], n);
		double K_p = K_precip(C_cu[i], s[i]);

		a[i] = -0.5 * (r_up * r_up) / (r_ave * r_ave * porosity * s[i] * dr) * (0.5 * F[i] - s_up * porosity * D / tort / dr);
		b[i] = -0.5 * (r_down * r_down) / (r_ave * r_ave * porosity * s[i] * dr) * (-0.5 * F[i - 1] - s_down * porosity * D / tort / dr);
		c[i] = -0.5 * ((r_up * r_up) / (r_ave * r_ave * porosity * s[i] * dr) * (0.5 * F[i] + s_up * porosity * D / tort / dr) - (r_down * r_down) / (r_ave * r_ave * porosity * s[i] * dr) * (-0.5 * F[i - 1] + s_down * porosity * D / tort / dr)) - 1 / dt - (K_r_sum + K_p) / (porosity * s[i]) - ds_dt / s[i];

		f[i] = -C_cu[i] / dt - (K_r_sum + K_p) * C_sat / (porosity * s[i]) + 0.5 * (r_up * r_up) / (r_ave * r_ave * porosity * s[i] * dr) * (0.5 * F_prev[i] - s_up * porosity * D / tort / dr) * C_cu[i + 1];
	}

	a[0] = 1.0;
	b[0] = 0.0;
	c[0] = -1.0;
	f[0] = 0.0;

	int i = r_max - 1;
	double r_ave = dr * i, r_down = r_ave - 0.5 * dr, r_up = r_ave + 0.5 * dr;
	double s_down = (s[i] + s[i - 1]) / 2.0;
	double ds_dt = -(1.0 / (r_ave * r_ave * porosity)) * (F[i] * r_up * r_up - F[i - 1] * r_down * r_down) / dr;

	double K_r_sum = 0.0;
	for (int n = 0; n < num_grain_classes; n++)
		K_r_sum += K_react(C_cu[i], C_grain[n][i], s[i], n);
	double K_p = K_precip(C_cu[i], s[i]);

	a[i] = 0.0;
	b[i] = -(r_down * r_down) / (r_ave * r_ave * porosity * s[i] * dr) * (-0.5 * F[i - 1] - s_down * porosity * D / tort / dr);
	c[i] = -(r_down * r_down) / (r_ave * r_ave * porosity * s[i] * dr) * (-0.5 * F[i - 1] + s_down * porosity * D / tort / dr) - 1 / dt - (K_r_sum + K_p) / (porosity * s[i]) - ds_dt / s[i];

	f[i] = -C_cu[i] / dt - (K_r_sum + K_p) * C_sat / (porosity * s[i]);
}

void update_values_grain(void)
{
	for (int n = 0; n < num_grain_classes; n++)
	{
		for (int i = 1; i < r_max; i++)
		{

			double K_r = K_react(C_cu[i], C_grain[n][i], s[i], n);
			double K_p = K_precip(C_cu[i], s[i]);

			C_grain[n][i] -= K_r * dt * (C_sat - C_cu[i]);
			C_precip[i] += K_p * dt * (C_cu[i] - C_sat);

			if (C_grain[n][i] < 0.0)
				C_grain[n][i] = 0.0;
		}

		C_grain[n][0] = C_grain[n][1];
	}
	C_precip[0] = C_precip[1];
}

void main(void)
{
	double t_write = dt_out;
	fstream f1, f2, f3, f4;

	setup_value();

	// will need to change file format for 9 grain classes

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
		calc_flux_water();
		update_values_saturation();

		setup_conc_matrix();
		Tri_Diag(a.data(), b.data(), c.data(), f.data(), C_cu.data(), r_max);
		update_values_grain();

		if (t > t_write)
		{
			cout << t << "\ts_e: " << s[r_max - 1] << "\ts_m: " << s[(r_max / 2)] << "\tC_cu_e: " << C_cu[r_max - 1] << "\tC_cu_m: " << C_cu[(r_max / 2)] << "\tC_grain_e: " << C_grain[r_max - 1] << "\tC_grain_m: " << C_grain[(r_max / 2)] << endl;

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
}

/*void calc_flux_water(void)
{
	for (int i = 0; i < r_max; i++)
		Pc[i] = P(s[i]);

	for (int i = 0; i < r_max - 1; i++)
	{
		double s_ave = (s[i] + s[i + 1]) / 2.0;
		F[i] = -(s_ave * porosity * k(s_ave) / mu) * (Pc[i + 1] - Pc[i]) / dr;

		double C_av = (C_cu[i] + C_cu[i + 1]) / 2.0;

		F_cu[i] = F[i] * C_av - s_ave * porosity * D / tort * (C_cu[i + 1] - C_cu[i]) / dr;

		//cout << -(s[i] * porosity * k(s_ave) / mu) * (Pc[i + 1] - Pc[i]) / dr <<endl;
	}

	F[r_max - 1] = std::max(s[r_max - 1] - s_min, 0.0) * porosity * v_evap;
	F_cu[r_max - 1] = 0.0;

	//cout << F[r_max - 1] << " : " << (1.0 / (radius * radius * porosity)) * ( F[r_max - 1] * radius * radius) / dr << endl;
}

void update_values_saturation(void)
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
}*/