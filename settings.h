#pragma once

//water
#define num_grain_classes	9


//surface tension divided by pore size
#define c_cap	2500
#define k_max	1.0e-11
#define mu		0.001
#define v_evap	1.0e-8
#define radius	0.005
#define s_min	0.05
#define n_pc	2.0
#define n_k		2.0

#define sat_0	1

//agglomerate
#define porosity	0.05
#define r_max		10.0

//mass of grains per vol ore
#define C_grain_0	2.0

#define C_sat		150.0

#define Beta		1.0

//#define k_react		1e-3

#define k_precip	1.0

#define k_react_0	1.0e-6

//reactant
#define D		1.0e-4
//#define D		10e-5
#define tort	7.0

//system
//#define t_max 3000000.0
#define t_max 2000000.0
#define dt_out	1800.0

#define cfl_num		0.3
//#define dt		1e-3


