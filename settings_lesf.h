#pragma once

//water

//surface tension divided by pore size
#define c_cap	2500
#define k_max	1.0e-11
#define mu		0.01
#define v_evap	1.0e-8
#define radius	0.005
#define s_min	0.05
#define n_pc	2.0
#define n_k		2.0

#define sat_0	1

//agglomorate
#define porosity	0.05
#define r_max		20

//mass of grains per vol ore
#define C_grain_0	1.0

#define C_sat		50.0

#define Beta		1.0

//#define k_react		1e-3

#define k_precip	1.0

#define k_react	1.0e-6

//reactant
#define D		1.0e-9
#define tort	7.0

//system
#define t_max	2600000.0
//#define t_max 3000.0
#define dt_out	300.0

#define cfl_num		0.3
//#define dt		1e-3


