#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED

using namespace std;

const int MAX_SIZE_DIM4 = 10;


const int N1 = 3;
const int N2 = 3;
const int N3 = 3;
const int N4 = 8; // N4 == 4*(int)
const int N_loops = 2;
const float tau = 10.0;
const float T_preequilibration = 10*tau;
const float T_equilibration = 30*tau;
const float T_measurement = 100*tau;
const int threads_per_beta = 8;
const unsigned int multihit_number = 8;

const float level_I_tau = 5.0;
const int level_I_multihit_number = 8;
const int level_I_number_of_measurements = 10;

const float level_II_tau = 5.0;
const int level_II_multihit_number = 8;
const int level_II_number_of_measurements = 10;

const int loop_number_of_measurements = 8;

const int iBeta_low = 580;
const int iBeta_step = 50;
const int iBeta_high = 620;
const bool first_file_initialization = false;



#endif // MAIN_H_INCLUDED
