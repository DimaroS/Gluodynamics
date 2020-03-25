#include <iostream>
#include <chrono>
#include <random>
#include <thread>

#include "main.h"
#include "MatrixSU2.h"
#include "MatrixSU3.h"
#include "CycledArrayDim4.h"
#include "GluodynamicsDim4_SU.h"

using namespace std;

//1.73458666666667 for 3x3x3x8 pSU3

//new comment

typedef PeriodicGluodynamicsDim4_SU<MatrixSU2, DynamicUnsafeArrayDim4,
                    float, ranlux24> Gluodynamics;
const int matr_dim = 2;
const char bc_code = 'p';
const char system_preconfiguration_file_format[500] =
    "/media/user/LinuxDATA/GluodynamicsDATA/%d%d%d%d%cSU%dSystem_%d.save";
const char system_log_file_for_beta_format[500] =
    "/media/user/LinuxDATA/GluodynamicsDATA/%d%d%d%doutput_%cSU%d_%u_%d_%d.csv";
const char system_log_file_for_beta_preeq_format[500] =
    "/media/user/LinuxDATA/GluodynamicsDATA/%d%d%d%doutput_%cSU%d_%u_%d.csv";
const char system_configuration_file_format[500] =
    "/media/user/LinuxDATA/GluodynamicsDATA/%d_%d%d%d%dSystem_%cSU%d_%d_%d.save";
const char system_measurements_file_format[500] =
    "/media/user/LinuxDATA/GluodynamicsDATA/Measured_%d%d%d%doutput_%cSU%d_%u_%d.csv";
const char system_single_measurement_file_format[500] =
    "/media/user/LinuxDATA/GluodynamicsDATA/Measured%d_%d%d%d%dSystem_%cSU%d_%d_%d.save";
//   /media/user/LinuxDATA/
//   /home/itep/sychev/

void EquilibrateForBeta(unsigned int key, int iBeta, int seed, int thread_id) {
    float beta = 0.01*iBeta;
    float g0 = sqrt(2*(matr_dim)/beta);
    Gluodynamics System(N1, N2, N3, N4, 0, g0, key);

    char file_name[500];
    sprintf(file_name, system_log_file_for_beta_format, N1, N2, N3, N4, bc_code, matr_dim, key, iBeta, thread_id);
    ofstream out_loc(file_name, ios::out | ios::app);


    char file_name_in[500];
    sprintf(file_name_in, system_preconfiguration_file_format, N1, N2, N3, N4, bc_code, matr_dim, iBeta);
    ifstream in_sys(file_name_in, ios::in | ios::binary);
    in_sys >> System;
    in_sys.close();
    System.Seed(seed);



    for (int t = 0; t < int(T_equilibration/tau); t++) {
        chrono::time_point<chrono::high_resolution_clock>
                        start_time = chrono::high_resolution_clock::now();

        float hit_ratio = System.MonteCarloStep(tau, multihit_number);


        chrono::time_point<chrono::high_resolution_clock>
                        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> diff_time = end_time - start_time;


        double s = 1.0 + System.Action()*g0*g0/2/(matr_dim)/N1/N2/N3/N4/6;

        cout << "eq " << beta << '\t' << t + int(T_preequilibration/tau) << '\t' << s << '\t'
                << hit_ratio << "\t\t"
                << (int) 4*tau*N1*N2*N3*N4/hit_ratio/diff_time.count()
                << " hits/sec"<< endl;

        out_loc << "eq," << beta << ',' << t + int(T_preequilibration/tau) << ',' << s << ',' << hit_ratio << ','
            << (int) 4*tau*N1*N2*N3*N4/hit_ratio/diff_time.count() << '\n';
        out_loc.close();
        out_loc.open(file_name, ios::out | ios::app);
    }



    for (int t = 0; t < int(T_measurement/tau); t++) {
        chrono::time_point<chrono::high_resolution_clock>
                        start_time = chrono::high_resolution_clock::now();

        float hit_ratio = System.MonteCarloStep(tau, multihit_number);


        chrono::time_point<chrono::high_resolution_clock>
                        end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> diff_time = end_time - start_time;


        double s = 1.0 + System.Action()*g0*g0/2/(matr_dim)/N1/N2/N3/N4/6;


        char file_name_out[500];
        sprintf(file_name_out, system_configuration_file_format, iBeta, N1, N2, N3, N4, bc_code, matr_dim, thread_id, t);
        ofstream out_sys(file_name_out, ios::out | ios::binary | ios::trunc);
        out_sys << System;
        out_sys.close();



        cout << beta << '\t' << t << '\t' << s << '\t'
                << hit_ratio << "\t\t"
                << (int) 4*tau*N1*N2*N3*N4/hit_ratio/diff_time.count()
                << " hits/sec"<< endl;

        out_loc << beta << ',' << t << ',' << s << ',' << hit_ratio << ','
            << (int) 4*tau*N1*N2*N3*N4/hit_ratio/diff_time.count() << '\n';
        out_loc.close();
        out_loc.open(file_name, ios::out | ios::app);
    }

    out_loc.close();
}

void Equilibration(unsigned int key) {
    for (int iBeta = iBeta_low; iBeta <= iBeta_high; iBeta += iBeta_step) {
        float beta = 0.01*iBeta;
        float g0 = sqrt(2*(matr_dim)/beta);
        Gluodynamics System(N1, N2, N3, N4, 0, g0, key);

        char file_name[500];
        sprintf(file_name, system_log_file_for_beta_preeq_format, N1, N2, N3, N4,
                                                                    bc_code, matr_dim, key, iBeta);
        ofstream out_loc(file_name, ios::out | ios::app);


        if (iBeta != iBeta_low || first_file_initialization) {
            char file_name_in[500];
            sprintf(file_name_in, system_preconfiguration_file_format, N1, N2, N3, N4,
                                                            bc_code, matr_dim, iBeta - iBeta_step);
            ifstream in_sys(file_name_in, ios::in | ios::binary);
            in_sys >> System;
            in_sys.close();
        }

        System.Set_beta(beta);
        System.Set_t(0.0);

        for (int t = 0; t < int(T_preequilibration/tau); t++) {
            chrono::time_point<chrono::high_resolution_clock>
                            start_time = chrono::high_resolution_clock::now();

            float hit_ratio = System.MonteCarloStep(tau, multihit_number);


            chrono::time_point<chrono::high_resolution_clock>
                            end_time = chrono::high_resolution_clock::now();
            chrono::duration<double> diff_time = end_time - start_time;


            double s = 1.0 + System.Action()*g0*g0/2/(matr_dim)/N1/N2/N3/N4/6;



            cout << "pre " << beta << '\t' << t << '\t' << s << '\t'
                    << hit_ratio << "\t\t"
                    << (int) 4*tau*N1*N2*N3*N4/hit_ratio/diff_time.count()
                    << " hits/sec"<< endl;

            out_loc << "pre," << beta << ',' << t << ',' << s << ',' << hit_ratio << ','
                << (int) 4*tau*N1*N2*N3*N4/hit_ratio/diff_time.count() << '\n';
            out_loc.close();
            out_loc.open(file_name, ios::out | ios::app);
        }

        char file_name_out[500];
        sprintf(file_name_out, system_preconfiguration_file_format, N1, N2, N3, N4, bc_code, matr_dim, iBeta);
        ofstream out_sys(file_name_out, ios::out | ios::binary | ios::trunc);
        out_sys << System;
        out_sys.close();

        out_loc.close();
    }

    thread beta_threads[(iBeta_high - iBeta_low)/iBeta_step + 1][threads_per_beta];
    mt19937 seed_gen(key);

    for (int t = 0; t <= (iBeta_high - iBeta_low)/iBeta_step; t++) {
        for (int thread_for_beta_id = 0; thread_for_beta_id < threads_per_beta; thread_for_beta_id++) {
            beta_threads[t][thread_for_beta_id] = thread(EquilibrateForBeta,
                key, iBeta_low + t*iBeta_step, seed_gen(), thread_for_beta_id);
        }
    }

    for (int t = 0; t <= (iBeta_high - iBeta_low)/iBeta_step; t++) {
        for (int thread_for_beta_id = 0; thread_for_beta_id < threads_per_beta; thread_for_beta_id++) {
            beta_threads[t][thread_for_beta_id].join();
        }
    }
}







void MeasureCreutzRatioForBeta_ReadingThreadFunction (unsigned int key,
            int iBeta, int thread_id, double *s, double ***W) {
    float beta = 0.01*iBeta;
    float g0 = sqrt(2*(matr_dim)/beta);
    Gluodynamics System(N1, N2, N3, N4, 0, g0, key);

    char file_name[500];
    sprintf(file_name, system_log_file_for_beta_format, N1, N2, N3, N4,
                                    bc_code, matr_dim, key, iBeta, thread_id);
    ofstream out_loc(file_name, ios::out | ios::app);


    for (int t = 0; t < int(T_measurement/tau); t++) {
        char file_name_in[500];

        sprintf(file_name_in,
        system_configuration_file_format,
                iBeta, N1, N2, N3, N4, bc_code, matr_dim, thread_id, t);

        ifstream in_sys(file_name_in, ios::in | ios::binary);
        in_sys >> System;
        in_sys.close();

        s[t] = 1.0 + System.Action()*g0*g0/2/(matr_dim)/N1/N2/N3/N4/6;

        for (int i = 0; i < N_loops; i++) {
            for (int j = 0; j <= i; j++) {
                W[i][j][t] = System.AverageWilsonLoop(i+1, j+1);
            }
        }



        cout << beta << '\t' << thread_id << '\t' << t  << '\t'  << s[t] << '\t'
                << endl;

        out_loc << beta << ',' << t << ',' << s[t]
                    << '\n';
        out_loc.close();
        out_loc.open(file_name, ios::out | ios::app);
    }

    out_loc.close();
}





double Average (double **s, int i_max, int j_max) {
    double s_av = 0.0;

    for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
            s_av += s[i][j];
        }
    }
    s_av /= i_max*j_max;

    return s_av;
}



double Variance (double **s, int i_max, int j_max) {
    double s_av = Average(s, i_max, j_max);
    double s_var = 0.0;

    for (int i = 0; i < i_max; i++) {
        for (int j = 0; j < j_max; j++) {
            s_var += (s[i][j] - s_av)*(s[i][j] - s_av);
        }
    }
    s_var /= i_max*j_max - 1;

    return s_var;
}






void ActionAutocorrelation(double &s_av, double &s_d, double &tau_corr_av, double **s) {
    s_av = Average(s, threads_per_beta, int(T_measurement/tau));


    s_d = Variance(s, threads_per_beta, int(T_measurement/tau));
    s_d = sqrt(s_d);



    double **s_timecorr_normed = new double*[int(T_measurement/3/tau) + 1];
    for (int t = 0; t < int(T_measurement/3/tau) + 1; t++) {
        s_timecorr_normed[t] = new double[threads_per_beta];
    }
    int t_corr[threads_per_beta];
    double tau_corr[threads_per_beta];

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {

        for (int t = 0; t < int(T_measurement/3/tau); t++) {
            s_timecorr_normed[t][thread_id] = 0.0;
        }

        for (int t = 0; t < int(T_measurement/3/tau); t++) {
            for (int _tau = 0; _tau < int(T_measurement/tau) - t; _tau++) {
                s_timecorr_normed[t][thread_id] += (s[thread_id][_tau] - s_av)*(s[thread_id][_tau + t] - s_av);
            }
            s_timecorr_normed[t][thread_id] /= T_measurement/tau - t;
            s_timecorr_normed[t][thread_id] /= s_av*s_av;
        }

        for (int t = 1; t < int(T_measurement/3/tau); t++) {
            s_timecorr_normed[t][thread_id] /= s_timecorr_normed[0][thread_id];
        }
        s_timecorr_normed[0][thread_id] = 1.0;


        for (t_corr[thread_id] = 0;
                s_timecorr_normed[t_corr[thread_id]][thread_id] > 0.368
                && t_corr[thread_id] < int(T_measurement/3/tau) - 1;
                                                        t_corr[thread_id]++);

        tau_corr[thread_id] = 0.0;
        if (t_corr[thread_id] <= 1) {
            tau_corr[thread_id] = 1.0;
        } else {
            tau_corr[thread_id] = 0.5;
            for (int t = 1; t < t_corr[thread_id]-1; t++) {
                tau_corr[thread_id] += s_timecorr_normed[t][thread_id];
            }
            tau_corr[thread_id] += s_timecorr_normed[t_corr[thread_id]-1][thread_id]/2;
            tau_corr[thread_id] += (s_timecorr_normed[t_corr[thread_id]-1][thread_id]*s_timecorr_normed[t_corr[thread_id]-1][thread_id] - 0.135)
                    /(s_timecorr_normed[t_corr[thread_id]-1][thread_id] - s_timecorr_normed[t_corr[thread_id]][thread_id])/2;
            tau_corr[thread_id] /= 0.632;
        }
    }

    tau_corr_av = 0.0;
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        tau_corr_av += tau_corr[thread_id];
    }
    tau_corr_av /= threads_per_beta;

    for (int t = 0; t < int(T_measurement/3/tau) + 1; t++) {
        delete [] s_timecorr_normed[t];
    }
    delete [] s_timecorr_normed;

}




void MeasureCreutzRatioForBeta(unsigned int key, int iBeta) {
    float beta = 0.01*iBeta;


    char output_file_name[500];
    sprintf(output_file_name, system_measurements_file_format, N1, N2, N3, N4, bc_code, matr_dim, key, iBeta);
    ofstream out(output_file_name, ios::out | ios::app);



    double **s = new double*[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        s[thread_id] = new double[int(T_measurement/tau) + 1];
    }

    double ****W = new double***[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        W[thread_id] = new double**[N_loops];
        for (int i = 0; i < N_loops; i++) {
            W[thread_id][i] = new double*[N_loops];
            for (int j = 0; j < N_loops; j++) {
                W[thread_id][i][j] = new double[int(T_measurement/tau) + 1];
            }
        }
    }



    thread thread_array[threads_per_beta];

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        thread_array[thread_id] =
            thread(MeasureCreutzRatioForBeta_ReadingThreadFunction, key,
                iBeta, thread_id, s[thread_id], W[thread_id]);
    }

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        thread_array[thread_id].join();
    }






    double s_av = 0.0;
    double s_d = 0.0;
    double tau_corr_av = 0.0;
    ActionAutocorrelation(s_av, s_d, tau_corr_av, s);




    double W_jack[threads_per_beta][N_loops][N_loops][int(T_measurement/tau) + 1];
    double chi_jack[threads_per_beta][N_loops][N_loops][int(T_measurement/tau) + 1];

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int t_jack = 0; t_jack < int(T_measurement/tau); t_jack++) {
            for (int i = 0; i < N_loops; i++) {
                for (int j = 0; j < N_loops; j++) {
                    W_jack[thread_id][i][j][t_jack] = 0.0;
                }
            }

            for (int i = 0; i < N_loops; i++) {
                for (int j = 0; j <= i; j++) {
                    for (int t = 0; t < int(T_measurement/tau); t++) {
                        if (t == t_jack) {
                            continue;
                        }

                        W_jack[thread_id][i][j][t_jack] += W[thread_id][i][j][t];
                    }
                }
            }
            for (int i = 0; i < N_loops; i++) {
                for (int j = 0; j <= i; j++) {
                    W_jack[thread_id][i][j][t_jack] /= T_measurement/tau - 1;
                    W_jack[thread_id][j][i][t_jack] = W_jack[thread_id][i][j][t_jack];
                }
            }


            for (int i = 1; i < N_loops; i++) {
                for (int j = 1; j <= i; j++) {
                    chi_jack[thread_id][i][j][t_jack] =
                    W_jack[thread_id][i][j][t_jack]*W_jack[thread_id][i-1][j-1][t_jack]
                     /W_jack[thread_id][i-1][j][t_jack]/W_jack[thread_id][i][j-1][t_jack] > 0 ?
                     -log(W_jack[thread_id][i][j][t_jack]*W_jack[thread_id][i-1][j-1][t_jack]
                     /W_jack[thread_id][i-1][j][t_jack]/W_jack[thread_id][i][j-1][t_jack]) :
                      0.0;
                }
            }
        }
    }


    double W_av[N_loops][N_loops];
    double chi_av[N_loops][N_loops];
    double W_d[N_loops][N_loops];
    double chi_d[N_loops][N_loops];

    for (int i = 0; i < N_loops; i++) {
        for (int j = 0; j < N_loops; j++) {
            W_av[i][j] = chi_av[i][j] = W_d[i][j] = chi_d[i][j] = 0.0;
        }
    }

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int t = 0; t < int(T_measurement/tau); t++) {
            for (int i = 0; i < N_loops; i++) {
                for (int j = 0; j < N_loops; j++) {
                    W_av[i][j] += W_jack[thread_id][i][j][t];
                    chi_av[i][j] += chi_jack[thread_id][i][j][t];
                }
            }
        }
    }
    for (int i = 0; i < N_loops; i++) {
        for (int j = 0; j < N_loops; j++) {
            W_av[i][j] /= threads_per_beta*T_measurement/tau;
            chi_av[i][j] /= threads_per_beta*T_measurement/tau;
        }
    }

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int t = 0; t < int(T_measurement/tau); t++) {
            for (int i = 0; i < N_loops; i++) {
                for (int j = 0; j < N_loops; j++) {
                    W_d[i][j] += (W_jack[thread_id][i][j][t] - W_av[i][j]) * (W_jack[thread_id][i][j][t] - W_av[i][j]);
                    chi_d[i][j] += (chi_jack[thread_id][i][j][t] - chi_av[i][j]) * (chi_jack[thread_id][i][j][t] - chi_av[i][j]);
                }
            }
        }
    }
    for (int i = 0; i < N_loops; i++) {
        for (int j = 0; j < N_loops; j++) {
            W_d[i][j] /= threads_per_beta*(T_measurement/tau)/(T_measurement/tau - 1);
            W_d[i][j] = sqrt(W_d[i][j]);
            chi_d[i][j] /= threads_per_beta*(T_measurement/tau)/(T_measurement/tau - 1);
            chi_d[i][j] = sqrt(chi_d[i][j]);
        }
    }









    out << beta << ',' << s_av << ',' << s_d << ',' << tau_corr_av << '\n';
    out.close();
    out.open(output_file_name, ios::out | ios::app);


    for (int i = 0; i < N_loops; i++) {
        for (int j = 0; j <= i; j++) {
            out << beta << ',' << tau_corr_av*tau << ',' << i+1 << ',' << j+1
                << ',' << W_av[i][j] << ',' << W_d[i][j] << '\n';
        }
    }
    out.close();
    out.open(output_file_name, ios::out | ios::app);


    for (int i = 1; i < N_loops; i++) {
        for (int j = 1; j <= i; j++) {
            out << beta << ',' << tau_corr_av*tau << ',' << i+1 << ',' << j+1
                << ',' << chi_av[i][j] << ','
                << chi_d[i][j] << ',' << W_av[i][j] << ',' << W_d[i][j] << ','
                <<  W_av[i-1][j-1]<< ',' << W_d[i-1][j-1] << ',' << W_av[i-1][j]
                << ',' << W_d[i-1][j] << ','
                << W_av[i][j-1] << ',' << W_d[i][j-1] << '\n';
        }
    }

    out.close();


    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        delete [] s[thread_id];
    }
    delete [] s;

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int i = 0; i < N_loops; i++) {
            for (int j = 0; j < N_loops; j++) {
                delete [] W[thread_id][i][j];
            }
            delete [] W[thread_id][i];
        }
        delete [] W[thread_id];
    }
    delete [] W;
}







void MeasureAverageOrientedWilsonLoopForBeta_ReadingThreadFunction (unsigned int key,
            int iBeta, int thread_id, double *s, double ***W) {
    float beta = 0.01*iBeta;
    float g0 = sqrt(2*(matr_dim)/beta);
    Gluodynamics System(N1, N2, N3, N4, 0, g0, key);

    char file_name[500];
    sprintf(file_name, system_log_file_for_beta_format, N1, N2, N3, N4,
                                    bc_code, matr_dim, key, iBeta, thread_id);
    ofstream out_loc(file_name, ios::out | ios::app);


    for (int t = 0; t < int(T_measurement/tau); t++) {
        char file_name_in[500];

        sprintf(file_name_in,
        system_configuration_file_format,
                iBeta, N1, N2, N3, N4, bc_code, matr_dim, thread_id, t);

        ifstream in_sys(file_name_in, ios::in | ios::binary);
        in_sys >> System;
        in_sys.close();



        char file_name_out_measured[500];

        sprintf(file_name_out_measured,
        system_single_measurement_file_format,
                iBeta, N1, N2, N3, N4, bc_code, matr_dim, thread_id, t);

        ofstream out_mes;

        out_mes.open(file_name_out_measured, ios::out | ios::binary);







        s[t] = 1.0 + System.Action()*g0*g0/2/(matr_dim)/N1/N2/N3/N4/6;

        for (int i_direction = 0; i_direction <= 3; i_direction++) {
            for (int j_direction = 0; j_direction <= 3; j_direction++) {
                if (i_direction == j_direction) {
                    continue;
                }

                W[i_direction][j_direction][t] = System.AverageOrientedWilsonLoop(1, 1,
                                                i_direction, j_direction);
            }
        }


        for (int i_direction = 0; i_direction <= 3; i_direction++) {
            for (int j_direction = 0; j_direction <= 3; j_direction++) {
                if (i_direction == j_direction) {
                    continue;
                }

                out_mes.write((const char *) &W[i_direction][j_direction][t],
                            sizeof(double));
            }
        }
        out_mes.write((const char *) &s[t], sizeof(double));
        out_mes.close();



        cout << beta << '\t' << thread_id << '\t' << t  << '\t'  << s[t] << '\t'
                << endl;

        out_loc << beta << ',' << t << ',' << s[t]
                    << '\n';
        out_loc.close();
        out_loc.open(file_name, ios::out | ios::app);
    }

    out_loc.close();
}





void MeasureAverageOrientedWilsonLoopForBeta(unsigned int key, int iBeta) {
    float beta = 0.01*iBeta;


    char output_file_name[500];
    sprintf(output_file_name, system_measurements_file_format, N1, N2, N3, N4, bc_code, matr_dim, key, iBeta);
    ofstream out(output_file_name, ios::out | ios::app);



    double **s = new double*[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        s[thread_id] = new double[int(T_measurement/tau) + 1];
    }

    double ****W = new double***[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        W[thread_id] = new double**[4];
        for (int i_direction = 0; i_direction < 4; i_direction++) {
            W[thread_id][i_direction] = new double*[4];
            for (int j_direction = 0; j_direction < 4; j_direction++) {
                W[thread_id][i_direction][j_direction]
                                    = new double[int(T_measurement/tau) + 1];
            }
        }
    }



    thread thread_array[threads_per_beta];

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        thread_array[thread_id] =
            thread(MeasureAverageOrientedWilsonLoopForBeta_ReadingThreadFunction, key,
                iBeta, thread_id, s[thread_id], W[thread_id]);
    }

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        thread_array[thread_id].join();
    }






    double s_av = 0.0;
    double s_d = 0.0;
    double tau_corr_av = 0.0;
    ActionAutocorrelation(s_av, s_d, tau_corr_av, s);

//    double tau_corr_d = 0.0;
//    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
//        tau_corr_d += (tau_corr[thread_id] - tau_corr_av)*
//                        (tau_corr[thread_id] - tau_corr_av);
//    }
//    tau_corr_d /= threads_per_beta - 1;
//    tau_corr_d = sqrt(tau_corr_d);

















    double W_av[4][4];
    double W_d[4][4];

    for (int i_direction = 0; i_direction < 4; i_direction++) {
        for (int j_direction = 0; j_direction < 4; j_direction++) {
            W_av[i_direction][j_direction] = W_d[i_direction][j_direction] = 0.0;
        }
    }

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int t = 0; t < int(T_measurement/tau); t++) {
            for (int i_direction = 0; i_direction < 4; i_direction++) {
                for (int j_direction = 0; j_direction < 4; j_direction++) {
                    if (i_direction == j_direction) {
                        continue;
                    }

                    W_av[i_direction][j_direction] += W[thread_id][i_direction][j_direction][t];
                }
            }
        }
    }
    for (int i_direction = 0; i_direction < 4; i_direction++) {
        for (int j_direction = 0; j_direction < 4; j_direction++) {
            if (i_direction == j_direction) {
                continue;
            }

            W_av[i_direction][j_direction] /= threads_per_beta*T_measurement/tau;
        }
    }

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int t = 0; t < int(T_measurement/tau); t++) {
            for (int i_direction = 0; i_direction < 4; i_direction++) {
                for (int j_direction = 0; j_direction < 4; j_direction++) {
                    if (i_direction == j_direction) {
                        continue;
                    }

                    W_d[i_direction][j_direction] += (W[thread_id][i_direction][j_direction][t] - W_av[i_direction][j_direction])
                                        *(W[thread_id][i_direction][j_direction][t] - W_av[i_direction][j_direction]);
                }
            }
        }
    }
    for (int i_direction = 0; i_direction < 4; i_direction++) {
        for (int j_direction = 0; j_direction < 4; j_direction++) {
            if (i_direction == j_direction) {
                continue;
            }

            W_d[i_direction][j_direction] /= threads_per_beta*T_measurement/tau*(threads_per_beta*T_measurement/tau - 1);
            W_d[i_direction][j_direction] = sqrt(W_d[i_direction][j_direction]);
        }
    }






    out << beta << ',' << tau_corr_av << /*',' << tau_corr_d << */'\n';

//    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
//        out << beta << ',' << thread_id << ',' << tau_corr[thread_id] << '\n';
//    }

    out << beta << ',' << s_av << ',' << s_d << ',' << tau_corr_av << '\n';
    out.close();
    out.open(output_file_name, ios::out | ios::app);


    for (int i_direction = 0; i_direction < 4; i_direction++) {
        for (int j_direction = 0; j_direction < 4; j_direction++) {
            if (i_direction == j_direction) {
                continue;
            }

            out << beta << ',' << tau_corr_av*tau << ',' << i_direction << ',' << j_direction
                << ',' << W_av[i_direction][j_direction] << ',' << W_d[i_direction][j_direction] << '\n';
        }
    }


//    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
//        for (int t = 0; t < int(T_measurement/3/tau); t++) {
//            out << t << ',' << s_timecorr_normed[t][thread_id] << '\n';
//        }
//    }

    out.close();



    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        delete [] s[thread_id];
    }
    delete [] s;

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int i_direction = 0; i_direction < 4; i_direction++) {
            for (int j_direction = 0; j_direction < 4; j_direction++) {
                delete [] W[thread_id][i_direction][j_direction];
            }
            delete [] W[thread_id][i_direction];
        }
        delete [] W[thread_id];
    }
    delete [] W;
}





//  odd layers updates then even layers updates
void MeasureScalarGlueballForBeta_ReadingThreadFunction_v1 (unsigned int key,
            int iBeta, int thread_id, double *s, double **timeslice_energy) {
    float beta = 0.01*iBeta;
    float g0 = sqrt(2*(matr_dim)/beta);
    Gluodynamics System(N1, N2, N3, N4, 0, g0, key);

    char file_name[500];
    sprintf(file_name, system_log_file_for_beta_format, N1, N2, N3, N4,
                                    bc_code, matr_dim, key, iBeta, thread_id);
    ofstream out_loc(file_name, ios::out | ios::app);


    for (int t = 0; t < int(T_measurement/tau); t++) {
        char file_name_in[500];

        sprintf(file_name_in,
        system_configuration_file_format,
                iBeta, N1, N2, N3, N4, bc_code, matr_dim, thread_id, t);

        ifstream in_sys;

        in_sys.open(file_name_in, ios::in | ios::binary);
        in_sys >> System;
        in_sys.close();





        char file_name_out_measured[500];

        sprintf(file_name_out_measured,
        system_single_measurement_file_format,
                iBeta, N1, N2, N3, N4, bc_code, matr_dim, thread_id, t);

        ofstream out_mes;

        out_mes.open(file_name_out_measured, ios::out | ios::binary);





        s[t] = 1.0 + System.Action()*g0*g0/2/(matr_dim)/N1/N2/N3/N4/6;

        for (int timeslice = 0; timeslice < N4; timeslice += 2) {
            timeslice_energy[timeslice][t] = 0.0;
            for (int step = 0; step < level_I_number_of_measurements; step++) {
                for (int i = 0; i < N1; i++) {
                    for (int j = 0; j < N2; j++) {
                        for (int k = 0; k < N3; k++) {
                            timeslice_energy[timeslice][t] +=
    //                                        + System.SingleWilsonLoop(3, 3, 0, 1, i, j, k, timeslice)
                                            + System.SingleWilsonLoop(1, 1, 0, 1, i, j, k, timeslice)
    //                                        + System.SingleWilsonLoop(3, 3, 0, 2, i, j, k, timeslice)
                                            + System.SingleWilsonLoop(1, 1, 0, 2, i, j, k, timeslice)
    //                                        + System.SingleWilsonLoop(3, 3, 1, 2, i, j, k, timeslice)
                                            + System.SingleWilsonLoop(1, 1, 1, 2, i, j, k, timeslice);
                        }
                    }
                }

                System.Timeslice_I_MonteCarloStep(level_I_tau, level_I_multihit_number, timeslice, 0);
            }
            timeslice_energy[timeslice][t] /= level_I_number_of_measurements;
        }


        in_sys.open(file_name_in, ios::in | ios::binary);
        in_sys >> System;
        in_sys.close();



        for (int timeslice = 1; timeslice < N4; timeslice += 2) {
            timeslice_energy[timeslice][t] = 0.0;
            for (int step = 0; step < level_I_number_of_measurements; step++) {
                for (int i = 0; i < N1; i++) {
                    for (int j = 0; j < N2; j++) {
                        for (int k = 0; k < N3; k++) {
                            timeslice_energy[timeslice][t] +=
    //                                        + System.SingleWilsonLoop(3, 3, 0, 1, i, j, k, timeslice)
                                            + System.SingleWilsonLoop(1, 1, 0, 1, i, j, k, timeslice)
    //                                        + System.SingleWilsonLoop(3, 3, 0, 2, i, j, k, timeslice)
                                            + System.SingleWilsonLoop(1, 1, 0, 2, i, j, k, timeslice)
    //                                        + System.SingleWilsonLoop(3, 3, 1, 2, i, j, k, timeslice)
                                            + System.SingleWilsonLoop(1, 1, 1, 2, i, j, k, timeslice);
                        }
                    }
                }

                System.Timeslice_I_MonteCarloStep(level_I_tau, level_I_multihit_number, timeslice, 0);
            }
            timeslice_energy[timeslice][t] /= level_I_number_of_measurements;
        }

//
//        in_sys.open(file_name_in, ios::in | ios::binary);
//        in_sys >> System;
//        in_sys.close();
//
//
//
//        for (int timeslice = 2; timeslice < N4; timeslice += 4) {
//            timeslice_energy[timeslice][t] = 0.0;
//            for (int step = 0; step < level_I_number_of_measurements; step++) {
//                for (int i = 0; i < N1; i++) {
//                    for (int j = 0; j < N2; j++) {
//                        for (int k = 0; k < N3; k++) {
//                            timeslice_energy[timeslice][t] +=
//    //                                        + System.SingleWilsonLoop(3, 3, 0, 1, i, j, k, timeslice)
//                                            + System.SingleWilsonLoop(1, 1, 0, 1, i, j, k, timeslice)
//    //                                        + System.SingleWilsonLoop(3, 3, 0, 2, i, j, k, timeslice)
//                                            + System.SingleWilsonLoop(1, 1, 0, 2, i, j, k, timeslice)
//    //                                        + System.SingleWilsonLoop(3, 3, 1, 2, i, j, k, timeslice)
//                                            + System.SingleWilsonLoop(1, 1, 1, 2, i, j, k, timeslice);
//                        }
//                    }
//                }
//
//                System.Timeslice_I_MonteCarloStep(level_I_tau, level_I_multihit_number, timeslice, 2);
//            }
//            timeslice_energy[timeslice][t] /= level_I_number_of_measurements;
//        }
//
//
//        in_sys.open(file_name_in, ios::in | ios::binary);
//        in_sys >> System;
//        in_sys.close();
//
//
//        for (int timeslice = 3; timeslice < N4; timeslice += 4) {
//            timeslice_energy[timeslice][t] = 0.0;
////            double se[level_I_number_of_measurements];
//            for (int step = 0; step < level_I_number_of_measurements; step++) {
////                se[step] = -timeslice_energy[timeslice][t];
//                for (int i = 0; i < N1; i++) {
//                    for (int j = 0; j < N2; j++) {
//                        for (int k = 0; k < N3; k++) {
//                            timeslice_energy[timeslice][t] +=
//    //                                        + System.SingleWilsonLoop(3, 3, 0, 1, i, j, k, timeslice)
//                                            + System.SingleWilsonLoop(1, 1, 0, 1, i, j, k, timeslice)
//    //                                        + System.SingleWilsonLoop(3, 3, 0, 2, i, j, k, timeslice)
//                                            + System.SingleWilsonLoop(1, 1, 0, 2, i, j, k, timeslice)
//    //                                        + System.SingleWilsonLoop(3, 3, 1, 2, i, j, k, timeslice)
//                                            + System.SingleWilsonLoop(1, 1, 1, 2, i, j, k, timeslice);
//                        }
//                    }
//                }
////                se[step] += timeslice_energy[timeslice][t];
//                System.Timeslice_I_MonteCarloStep(level_I_tau, level_I_multihit_number, timeslice, 0);
//            }
////
////            double se_av = 0.0;
////            for (int step = 0; step < level_I_number_of_measurements; step++) {
////                se_av += se[step];
////            }
////            se_av /= level_I_number_of_measurements;
////
////            double se_d = 0.0;
////            for (int step = 0; step < level_I_number_of_measurements; step++) {
////                se_d += (se[step] - se_av)*(se[step] - se_av);
////            }
////            se_d /= level_I_number_of_measurements - 1;
////            se_d = sqrt(se_d);
////
////
////
////            double se_timecorr_normed[level_I_number_of_measurements/3 + 1];
////            int t_corr;
////            double tau_corr;
////
////
////            for (int t = 0; t < level_I_number_of_measurements/3; t++) {
////                se_timecorr_normed[t] = 0.0;
////            }
////
////            for (int t = 0; t < level_I_number_of_measurements/3; t++) {
////                for (int _tau = 0; _tau < level_I_number_of_measurements/3 - t; _tau++) {
////                    se_timecorr_normed[t] += (se[_tau] - se_av)*(se[_tau + t] - se_av);
////                }
////                se_timecorr_normed[t] /= level_I_number_of_measurements/3 - t;
////                se_timecorr_normed[t] /= se_av*se_av;
////            }
////
////            for (int t = 1; t < level_I_number_of_measurements/3; t++) {
////                se_timecorr_normed[t] /= se_timecorr_normed[0];
////            }
////            se_timecorr_normed[0] = 1.0;
////
////
////            for (t_corr = 0;
////                    se_timecorr_normed[t_corr] > 0.368
////                    && t_corr < level_I_number_of_measurements/3 - 1;
////                                                            t_corr++);
////
////            tau_corr = 0.0;
////            if (t_corr <= 1) {
////                tau_corr = 1.0;
////            } else {
////                tau_corr = 0.5;
////                for (int t = 1; t < t_corr-1; t++) {
////                    tau_corr += se_timecorr_normed[t];
////                }
////                tau_corr += se_timecorr_normed[t_corr-1]/2;
////                tau_corr += (se_timecorr_normed[t_corr-1]*se_timecorr_normed[t_corr-1] - 0.135)
////                        /(se_timecorr_normed[t_corr-1] - se_timecorr_normed[t_corr])/2;
////                tau_corr /= 0.632;
////            }
////
////
////            for (int t = 0; t < level_I_number_of_measurements/3; t++) {
////                out_loc << t << ',' << se_timecorr_normed[t] << '\n';
////            }
////
////
//
//
//
//
////            cout << se_av << ' ' << se_d << '\n';
//
//
//            timeslice_energy[timeslice][t] /= level_I_number_of_measurements;
//        }

        for (int timeslice = 0; timeslice < N4; timeslice++) {
            out_mes.write((const char *) &timeslice_energy[timeslice][t],
                            sizeof(double));
        }
        out_mes.write((const char *) &s[t], sizeof(double));
        out_mes.close();



        cout << beta << '\t' << thread_id << '\t' << t  << '\t'  << s[t] << '\t'
                << endl;

        out_loc << beta << ',' << t << ',' << s[t]
                    << '\n';
        out_loc.close();
        out_loc.open(file_name, ios::out | ios::app);
    }

    out_loc.close();
}

void MeasureScalarGlueballForBeta_v1_processing (unsigned int key, int iBeta,
                double s_av, double s_d, double tau_corr_av,
                double ***timeslice_energy,
                double timeslice_energy_av, double timeslice_energy_d,
                double ***timeslice_energy_corr_product,
                double **timeslice_energy_corr_product_jack,
                double *timeslice_energy_jack,
                double **timeslice_energy_corr_normed_jack,
                double *timeslice_energy_corr_normed_av,
                double *timeslice_energy_corr_normed_d,
                double ***timeslice_energy_corr_normed_old,
                double *timeslice_energy_corr_normed_av_old,
                double *timeslice_energy_corr_normed_d_old
                                                                            ) {
    float beta = 0.01*iBeta;



    char output_file_name[500];
    sprintf(output_file_name, system_measurements_file_format, N1, N2, N3, N4, bc_code, matr_dim, key, iBeta);
    ofstream out(output_file_name, ios::out | ios::app);




    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int t = 0; t < int(T_measurement/tau); t++) {
            for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
            timeslice_energy_corr_product[thread_id][timeslice_diff][t] = 0.0;
                for (int timeslice_base = 0; timeslice_base < N4; timeslice_base++) {
                    timeslice_energy_corr_product[thread_id][timeslice_diff][t] +=
                        timeslice_energy[thread_id][timeslice_base][t]
                        *timeslice_energy[thread_id][(timeslice_base + timeslice_diff)%N4][t];
                }

                timeslice_energy_corr_product[thread_id][timeslice_diff][t] /= N4;
            }
        }
    }







    for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
        for (int thread_id_base = 0; thread_id_base < threads_per_beta; thread_id_base++) {
            timeslice_energy_corr_product_jack[timeslice_diff][thread_id_base] = 0.0;
            for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
                if (thread_id == thread_id_base) {
                    continue;
                }

                for (int t = 0; t < int(T_measurement/tau); t++) {
                    timeslice_energy_corr_product_jack[timeslice_diff][thread_id_base]
                                    += timeslice_energy_corr_product[thread_id][timeslice_diff][t];
                }
            }

            timeslice_energy_corr_product_jack[timeslice_diff][thread_id_base]
                                    /= (threads_per_beta - 1)*T_measurement/tau;
        }
    }


    for (int thread_id_base = 0; thread_id_base < threads_per_beta; thread_id_base++) {
        timeslice_energy_jack[thread_id_base] = 0.0;
        for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
            if (thread_id == thread_id_base) {
                continue;
            }

            for (int timeslice = 0; timeslice < N4; timeslice++) {
                for (int t = 0; t < int(T_measurement/tau); t++) {
                    timeslice_energy_jack[thread_id_base] += timeslice_energy[thread_id][timeslice][t];
                }
            }
        }

        timeslice_energy_jack[thread_id_base] /= (threads_per_beta - 1)*N4*T_measurement/tau;
    }



    for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
        for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
            timeslice_energy_corr_normed_jack[timeslice_diff][thread_id] =
                timeslice_energy_corr_product_jack[timeslice_diff][thread_id]
                /timeslice_energy_jack[thread_id]/timeslice_energy_jack[thread_id] - 1.0;
        }
    }


    for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
        timeslice_energy_corr_normed_av[timeslice_diff] = 0.0;
        for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
            timeslice_energy_corr_normed_av[timeslice_diff] +=
                timeslice_energy_corr_normed_jack[timeslice_diff][thread_id];
        }
        timeslice_energy_corr_normed_av[timeslice_diff] /= threads_per_beta;
    }



    for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
        timeslice_energy_corr_normed_d[timeslice_diff] = 0.0;
        for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
            timeslice_energy_corr_normed_d[timeslice_diff] +=
                (timeslice_energy_corr_normed_jack[timeslice_diff][thread_id] - timeslice_energy_corr_normed_av[timeslice_diff])
                *(timeslice_energy_corr_normed_jack[timeslice_diff][thread_id] - timeslice_energy_corr_normed_av[timeslice_diff]);
        }
        timeslice_energy_corr_normed_d[timeslice_diff] /= threads_per_beta/(threads_per_beta - 1);
        timeslice_energy_corr_normed_d[timeslice_diff] = sqrt(timeslice_energy_corr_normed_d[timeslice_diff]);
    }










    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int t = 0; t < int(T_measurement/tau); t++) {
            for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
                timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t] = 0.0;
                for (int timeslice_base = 0; timeslice_base < N4; timeslice_base++) {
                    timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t] +=
                        (timeslice_energy[thread_id][timeslice_base][t] - timeslice_energy_av)
                        *(timeslice_energy[thread_id][(timeslice_base + timeslice_diff)%N4][t] - timeslice_energy_av);
                }

                timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t] /= N4;
                timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t] /= timeslice_energy_av*timeslice_energy_av;
            }
        }
    }


    for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
        timeslice_energy_corr_normed_av_old[timeslice_diff] = 0.0;
        for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
            for (int t = 0; t < int(T_measurement/tau); t++) {
                timeslice_energy_corr_normed_av_old[timeslice_diff] +=
                    timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t];
            }
        }
        timeslice_energy_corr_normed_av_old[timeslice_diff] /= threads_per_beta*T_measurement/tau;
    }


    for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
        timeslice_energy_corr_normed_d_old[timeslice_diff] = 0.0;
        for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
            for (int t = 0; t < int(T_measurement/tau); t++) {
                timeslice_energy_corr_normed_d_old[timeslice_diff] +=
                    (timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t] - timeslice_energy_corr_normed_av[timeslice_diff])
                    *(timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t] - timeslice_energy_corr_normed_av[timeslice_diff]);
            }
        }
        timeslice_energy_corr_normed_d_old[timeslice_diff] /= (threads_per_beta*T_measurement/tau)*(threads_per_beta*T_measurement/tau - 1);
        timeslice_energy_corr_normed_d_old[timeslice_diff] = sqrt(timeslice_energy_corr_normed_d_old[timeslice_diff]);
    }



    out << beta << ',' << s_av << ',' << s_d << ',' << tau_corr_av*tau << '\n';
    out << beta << ',' << timeslice_energy_av << ',' << timeslice_energy_d
                    << ',' << timeslice_energy_d/timeslice_energy_av  << '\n';
    out.close();
    out.open(output_file_name, ios::out | ios::app);


    for (int timeslice_diff = 2; timeslice_diff <= N4/2; timeslice_diff++) {
            out << beta << ',' << tau_corr_av << ',' << timeslice_diff << ','
                << timeslice_energy_corr_normed_av[timeslice_diff] << ','
                << timeslice_energy_corr_normed_d[timeslice_diff] << '\n';
    }
    for (int timeslice_diff = 2; timeslice_diff <= N4/2; timeslice_diff++) {
            out << beta << ',' << tau_corr_av << ',' << timeslice_diff << ','
                << timeslice_energy_corr_normed_av_old[timeslice_diff] << ','
                << timeslice_energy_corr_normed_d_old[timeslice_diff] << '\n';
    }
    out.close();
    out.open(output_file_name, ios::out | ios::app);


    out.close();


}

//  odd layers updates then even layers updates
void MeasureScalarGlueballForBeta_v1_0(unsigned int key, int iBeta) {
//    float beta = 0.01*iBeta;



    char output_file_name[500];
    sprintf(output_file_name, system_measurements_file_format, N1, N2, N3, N4, bc_code, matr_dim, key, iBeta);
    ofstream out(output_file_name, ios::out | ios::app);


    double **s = new double*[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        s[thread_id] = new double[int(T_measurement/tau) + 1];
    }

    double ***timeslice_energy = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy[thread_id] = new double*[N4];
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            timeslice_energy[thread_id][timeslice] = new double[int(T_measurement/tau) + 1];
        }
    }




    thread thread_array[threads_per_beta];

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        thread_array[thread_id] =
            thread(MeasureScalarGlueballForBeta_ReadingThreadFunction_v1, key,
                iBeta, thread_id, s[thread_id], timeslice_energy[thread_id]);
    }

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        thread_array[thread_id].join();
    }





    double s_av = 0.0;
    double s_d = 0.0;
    double tau_corr_av = 0.0;
    ActionAutocorrelation(s_av, s_d, tau_corr_av, s);

//    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
//        for (int t = 0; t < int(T_measurement/3/tau); t++) {
//            out << t << ',' << s_timecorr_normed[t][thread_id] << '\n';
//        }
//    }



    double timeslice_energy_av = 0.0;
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            for (int t = 0; t < int(T_measurement/tau); t++) {
                timeslice_energy_av += timeslice_energy[thread_id][timeslice][t];
            }
        }
    }
    timeslice_energy_av /= threads_per_beta*N4*T_measurement/tau;

    double timeslice_energy_d = 0.0;
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            for (int t = 0; t < int(T_measurement/tau); t++) {
                timeslice_energy_d += (timeslice_energy[thread_id][timeslice][t] - timeslice_energy_av)*
                                        (timeslice_energy[thread_id][timeslice][t] - timeslice_energy_av);
            }
        }
    }
    timeslice_energy_d /= (threads_per_beta*N4*T_measurement/tau)*(threads_per_beta*N4*T_measurement/tau);
    timeslice_energy_d = sqrt(timeslice_energy_d);




    double ***timeslice_energy_corr_product = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy_corr_product[thread_id] = new double*[N4];
        for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
            timeslice_energy_corr_product[thread_id][timeslice_diff] = new double[int(T_measurement/tau) + 1];
        }
    }
    double **timeslice_energy_corr_product_jack = new double*[N4];
    for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
        timeslice_energy_corr_product_jack[timeslice_diff] = new double[threads_per_beta];
    }
    double *timeslice_energy_jack = new double[threads_per_beta];
    double **timeslice_energy_corr_normed_jack = new double*[N4];
    for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
        timeslice_energy_corr_normed_jack[timeslice_diff] = new double[threads_per_beta];
    }
    double *timeslice_energy_corr_normed_av = new double[N4];
    double *timeslice_energy_corr_normed_d = new double[N4];
    double ***timeslice_energy_corr_normed_old = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy_corr_normed_old[thread_id] = new double*[N4];
        for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
            timeslice_energy_corr_normed_old[thread_id][timeslice_diff] = new double[int(T_measurement/tau) + 1];
        }
    }
    double *timeslice_energy_corr_normed_av_old = new double[N4];
    double *timeslice_energy_corr_normed_d_old = new double[N4];


    MeasureScalarGlueballForBeta_v1_processing (key, iBeta,
                s_av, s_d, tau_corr_av,
                timeslice_energy,
                timeslice_energy_av, timeslice_energy_d,
                timeslice_energy_corr_product,
                timeslice_energy_corr_product_jack,
                timeslice_energy_jack,
                timeslice_energy_corr_normed_jack,
                timeslice_energy_corr_normed_av,
                timeslice_energy_corr_normed_d,
                timeslice_energy_corr_normed_old,
                timeslice_energy_corr_normed_av_old,
                timeslice_energy_corr_normed_d_old);


    delete [] timeslice_energy_jack;
    delete [] timeslice_energy_corr_normed_av;
    delete [] timeslice_energy_corr_normed_d;
    delete [] timeslice_energy_corr_normed_av_old;
    delete [] timeslice_energy_corr_normed_d_old;


    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        delete [] s[thread_id];
    }
    delete [] s;

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
            delete [] timeslice_energy_corr_product[thread_id][timeslice_diff];
        }
        delete [] timeslice_energy_corr_product[thread_id];
    }
    delete [] timeslice_energy_corr_product;

    for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
        delete [] timeslice_energy_corr_product_jack[timeslice_diff];
    }
    delete [] timeslice_energy_corr_product_jack;

    for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
        delete [] timeslice_energy_corr_normed_jack[timeslice_diff];
    }
    delete [] timeslice_energy_corr_normed_jack;

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
            delete [] timeslice_energy_corr_normed_old[thread_id][timeslice_diff];
        }
        delete [] timeslice_energy_corr_normed_old[thread_id];
    }
    delete [] timeslice_energy_corr_normed_old;



    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            delete [] timeslice_energy[thread_id][timeslice];
        }
        delete [] timeslice_energy[thread_id];
    }
    delete [] timeslice_energy;
}





void MeasureScalarGlueballForBeta_v1_0_precalculated(unsigned int key, int iBeta) {
//    float beta = 0.01*iBeta;




    char output_file_name[500];
    sprintf(output_file_name, system_measurements_file_format, N1, N2, N3, N4, bc_code, matr_dim, key, iBeta);
    ofstream out(output_file_name, ios::out | ios::app);


    double **s = new double*[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        s[thread_id] = new double[int(T_measurement/tau) + 1];
    }

    double ***timeslice_energy = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy[thread_id] = new double*[N4];
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            timeslice_energy[thread_id][timeslice] = new double[int(T_measurement/tau) + 1];
        }
    }



    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int t = 0; t < int(T_measurement/tau); t++) {
            char file_name_out_measured[500];

            sprintf(file_name_out_measured,
            system_single_measurement_file_format,
                    iBeta, N1, N2, N3, N4, bc_code, matr_dim, thread_id, t);

            ifstream out_mes;

            out_mes.open(file_name_out_measured, ios::in | ios::binary);

            for (int timeslice = 0; timeslice < N4; timeslice++) {
                out_mes.read((char *) &timeslice_energy[thread_id][timeslice][t],
                                sizeof(double));
            }
            out_mes.read((char *) &s[thread_id][t], sizeof(double));

            out_mes.close();
        }
    }





    double s_av = 0.0;
    double s_d = 0.0;
    double tau_corr_av = 0.0;
    ActionAutocorrelation(s_av, s_d, tau_corr_av, s);

//    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
//        for (int t = 0; t < int(T_measurement/3/tau); t++) {
//            out << t << ',' << s_timecorr_normed[t][thread_id] << '\n';
//        }
//    }

    double timeslice_energy_av = 0.0;
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            for (int t = 0; t < int(T_measurement/tau); t++) {
                timeslice_energy_av += timeslice_energy[thread_id][timeslice][t];
            }
        }
    }
    timeslice_energy_av /= threads_per_beta*N4*T_measurement/tau;

    double timeslice_energy_d = 0.0;
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            for (int t = 0; t < int(T_measurement/tau); t++) {
                timeslice_energy_d += (timeslice_energy[thread_id][timeslice][t] - timeslice_energy_av)*
                                        (timeslice_energy[thread_id][timeslice][t] - timeslice_energy_av);
            }
        }
    }
    timeslice_energy_d /= (threads_per_beta*N4*T_measurement/tau)*(threads_per_beta*N4*T_measurement/tau);
    timeslice_energy_d = sqrt(timeslice_energy_d);


    double ***timeslice_energy_corr_product = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy_corr_product[thread_id] = new double*[N4];
        for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
            timeslice_energy_corr_product[thread_id][timeslice_diff] = new double[int(T_measurement/tau) + 1];
        }
    }
    double **timeslice_energy_corr_product_jack = new double*[N4];
    for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
        timeslice_energy_corr_product_jack[timeslice_diff] = new double[threads_per_beta];
    }
    double *timeslice_energy_jack = new double[threads_per_beta];
    double **timeslice_energy_corr_normed_jack = new double*[N4];
    for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
        timeslice_energy_corr_normed_jack[timeslice_diff] = new double[threads_per_beta];
    }
    double *timeslice_energy_corr_normed_av = new double[N4];
    double *timeslice_energy_corr_normed_d = new double[N4];
    double ***timeslice_energy_corr_normed_old = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy_corr_normed_old[thread_id] = new double*[N4];
        for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
            timeslice_energy_corr_normed_old[thread_id][timeslice_diff] = new double[int(T_measurement/tau) + 1];
        }
    }
    double *timeslice_energy_corr_normed_av_old = new double[N4];
    double *timeslice_energy_corr_normed_d_old = new double[N4];


    MeasureScalarGlueballForBeta_v1_processing (key, iBeta,
                s_av, s_d, tau_corr_av,
                timeslice_energy,
                timeslice_energy_av, timeslice_energy_d,
                timeslice_energy_corr_product,
                timeslice_energy_corr_product_jack,
                timeslice_energy_jack,
                timeslice_energy_corr_normed_jack,
                timeslice_energy_corr_normed_av,
                timeslice_energy_corr_normed_d,
                timeslice_energy_corr_normed_old,
                timeslice_energy_corr_normed_av_old,
                timeslice_energy_corr_normed_d_old);


    delete [] timeslice_energy_jack;
    delete [] timeslice_energy_corr_normed_av;
    delete [] timeslice_energy_corr_normed_d;
    delete [] timeslice_energy_corr_normed_av_old;
    delete [] timeslice_energy_corr_normed_d_old;



//    for (int t_diff = 0; t_diff <= N4/2; t_diff++) {
//        for (int t = 0; t < int(T_measurement/3/tau); t++) {
//            delete [] timeslice_c_timecorr_normed[t_diff][t];
//        }
//        delete [] timeslice_c_timecorr_normed[t_diff];
//    }
//    delete [] timeslice_c_timecorr_normed;

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        delete [] s[thread_id];
    }
    delete [] s;


    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
            delete [] timeslice_energy_corr_product[thread_id][timeslice_diff];
        }
        delete [] timeslice_energy_corr_product[thread_id];
    }
    delete [] timeslice_energy_corr_product;

    for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
        delete [] timeslice_energy_corr_product_jack[timeslice_diff];
    }
    delete [] timeslice_energy_corr_product_jack;

    for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
        delete [] timeslice_energy_corr_normed_jack[timeslice_diff];
    }
    delete [] timeslice_energy_corr_normed_jack;

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice_diff = 0; timeslice_diff < N4; timeslice_diff++) {
            delete [] timeslice_energy_corr_normed_old[thread_id][timeslice_diff];
        }
        delete [] timeslice_energy_corr_normed_old[thread_id];
    }
    delete [] timeslice_energy_corr_normed_old;



    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            delete [] timeslice_energy[thread_id][timeslice];
        }
        delete [] timeslice_energy[thread_id];
    }
    delete [] timeslice_energy;
}





//  half-spaces updates
void MeasureScalarGlueballForBeta_ReadingThreadFunction_v2 (unsigned int key,
            int iBeta, int thread_id, double *s, double **timeslice_energy) {
    float beta = 0.01*iBeta;
    float g0 = sqrt(2*(matr_dim)/beta);
    Gluodynamics System(N1, N2, N3, N4, 0, g0, key);

    char file_name[500];
    sprintf(file_name, system_log_file_for_beta_format, N1, N2, N3, N4,
                                    bc_code, matr_dim, key, iBeta, thread_id);
    ofstream out_loc(file_name, ios::out | ios::app);


    for (int t = 0; t < int(T_measurement/tau); t++) {
        char file_name_in[500];

        sprintf(file_name_in,
        system_configuration_file_format,
                iBeta, N1, N2, N3, N4, bc_code, matr_dim, thread_id, t);

        ifstream in_sys;

        in_sys.open(file_name_in, ios::in | ios::binary);
        in_sys >> System;
        in_sys.close();




        char file_name_out_measured[500];

        sprintf(file_name_out_measured,
        system_single_measurement_file_format,
                iBeta, N1, N2, N3, N4, bc_code, matr_dim, thread_id, t);

        ofstream out_mes;

        out_mes.open(file_name_out_measured, ios::out | ios::binary);






        s[t] = 1.0 + System.Action()*g0*g0/2/(matr_dim)/N1/N2/N3/N4/6;

        for (int timeslice = 0; timeslice < N4; timeslice++) {
            timeslice_energy[timeslice][t] = 0.0;
        }

        for (int step = 0; step < level_II_number_of_measurements; step++) {
            for (int timeslice = 0; timeslice < N4; timeslice++) {
                for (int i = 0; i < N1; i++) {
                    for (int j = 0; j < N2; j++) {
                        for (int k = 0; k < N3; k++) {
                            timeslice_energy[timeslice][t] +=
    //                                        + System.SingleWilsonLoop(3, 3, 0, 1, i, j, k, timeslice)
                                            + System.SingleWilsonLoop(1, 1, 0, 1, i, j, k, timeslice)
    //                                        + System.SingleWilsonLoop(3, 3, 0, 2, i, j, k, timeslice)
                                            + System.SingleWilsonLoop(1, 1, 0, 2, i, j, k, timeslice)
    //                                        + System.SingleWilsonLoop(3, 3, 1, 2, i, j, k, timeslice)
                                            + System.SingleWilsonLoop(1, 1, 1, 2, i, j, k, timeslice);
                        }
                    }
                }
            }

            System.Timeslice_II_MonteCarloStep(level_II_tau, level_II_multihit_number, N4/4, N4/2);
            System.Timeslice_II_MonteCarloStep(level_II_tau, level_II_multihit_number, N4/4 + N4/2, N4/2);
        }


        for (int timeslice = 0; timeslice < N4; timeslice ++) {
            timeslice_energy[timeslice][t] /= level_II_number_of_measurements;
            out_mes.write((const char *) &timeslice_energy[timeslice][t],
                            sizeof(double));
        }
        out_mes.write((const char *) &s[t], sizeof(double));
        out_mes.close();




        cout << beta << '\t' << thread_id << '\t' << t  << '\t'  << s[t] << '\t'
                << endl;

        out_loc << beta << ',' << t << ',' << s[t]
                    << '\n';
        out_loc.close();
        out_loc.open(file_name, ios::out | ios::app);
    }

    out_loc.close();
}




void MeasureScalarGlueballForBeta_v2_processing(unsigned int key, int iBeta,
                double s_av, double s_d, double tau_corr_av,
                double ***timeslice_energy,
                double timeslice_energy_av, double timeslice_energy_d,
                double ***timeslice_energy_corr_product,
                double **timeslice_energy_corr_product_jack,
                double *timeslice_energy_jack,
                double **timeslice_energy_corr_normed_jack,
                double *timeslice_energy_corr_normed_av,
                double *timeslice_energy_corr_normed_d,
                double ***timeslice_energy_corr_normed_old,
                double *timeslice_energy_corr_normed_av_old,
                double *timeslice_energy_corr_normed_d_old
                                                                    ) {
    float beta = 0.01*iBeta;



    char output_file_name[500];
    sprintf(output_file_name, system_measurements_file_format, N1, N2, N3, N4, bc_code, matr_dim, key, iBeta);
    ofstream out(output_file_name, ios::out | ios::app);





    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int t = 0; t < int(T_measurement/tau); t++) {
            for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
                timeslice_energy_corr_product[thread_id][timeslice_diff][t] = 0.0;
            }
        }
    }
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int t = 0; t < int(T_measurement/tau); t++) {
            for (int timeslice_diff = 2; timeslice_diff <= N4/2; timeslice_diff++) {
                for (int timeslice_base = N4/2 + 1 - timeslice_diff; timeslice_base < N4/2; timeslice_base++) {
                    timeslice_energy_corr_product[thread_id][timeslice_diff][t] +=
                        timeslice_energy[thread_id][timeslice_base][t]
                        *timeslice_energy[thread_id][(timeslice_base + timeslice_diff)%N4][t];
                }
                for (int timeslice_base = N4 + 1 - timeslice_diff; timeslice_base < N4; timeslice_base++) {
                    timeslice_energy_corr_product[thread_id][timeslice_diff][t] +=
                        timeslice_energy[thread_id][timeslice_base][t]
                        *timeslice_energy[thread_id][(timeslice_base + timeslice_diff)%N4][t];
                }

                timeslice_energy_corr_product[thread_id][timeslice_diff][t] /= (timeslice_diff - 1)*2;
            }
        }
    }





    for (int timeslice_diff = 2; timeslice_diff <= N4/2; timeslice_diff++) {
        for (int thread_id_base = 0; thread_id_base < threads_per_beta; thread_id_base++) {
            timeslice_energy_corr_product_jack[timeslice_diff][thread_id_base] = 0.0;
            for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
                if (thread_id_base == thread_id) {
                    continue;
                }

                for (int t = 0; t < int(T_measurement/tau); t++) {
                    timeslice_energy_corr_product_jack[timeslice_diff][thread_id_base] +=
                                    timeslice_energy_corr_product[thread_id][timeslice_diff][t];
                }
            }
            timeslice_energy_corr_product_jack[timeslice_diff][thread_id_base] /= (threads_per_beta - 1)*T_measurement/tau;
        }
    }


    for (int thread_id_base = 0; thread_id_base < threads_per_beta; thread_id_base++) {
        timeslice_energy_jack[thread_id_base] = 0.0;

        for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
            if (thread_id_base == thread_id) {
                continue;
            }

            for (int timeslice = 0; timeslice < N4; timeslice++) {
                for (int t = 0; t < int(T_measurement/tau); t++) {
                    timeslice_energy_jack[thread_id_base] += timeslice_energy[thread_id][timeslice][t];
                }
            }
        }

        timeslice_energy_jack[thread_id_base] /= (threads_per_beta - 1)*N4*T_measurement/tau;
    }



    for (int timeslice_diff = 2; timeslice_diff <= N4/2; timeslice_diff++) {
        for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
            timeslice_energy_corr_normed_jack[timeslice_diff][thread_id] =
                timeslice_energy_corr_product_jack[timeslice_diff][thread_id]
                /timeslice_energy_jack[thread_id]/timeslice_energy_jack[thread_id] - 1.0;
        }
    }


    for (int timeslice_diff = 2; timeslice_diff <= N4/2; timeslice_diff++) {
        timeslice_energy_corr_normed_av[timeslice_diff] = 0.0;
        for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
            timeslice_energy_corr_normed_av[timeslice_diff] +=
                timeslice_energy_corr_normed_jack[timeslice_diff][thread_id];
        }
        timeslice_energy_corr_normed_av[timeslice_diff] /= threads_per_beta;
    }



    for (int timeslice_diff = 2; timeslice_diff <= N4/2; timeslice_diff++) {
        timeslice_energy_corr_normed_d[timeslice_diff] = 0.0;
        for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
            timeslice_energy_corr_normed_d[timeslice_diff] +=
                (timeslice_energy_corr_normed_jack[timeslice_diff][thread_id] - timeslice_energy_corr_normed_av[timeslice_diff])
                *(timeslice_energy_corr_normed_jack[timeslice_diff][thread_id] - timeslice_energy_corr_normed_av[timeslice_diff]);
        }

        timeslice_energy_corr_normed_d[timeslice_diff] /= threads_per_beta/(threads_per_beta - 1);
        timeslice_energy_corr_normed_d[timeslice_diff] = sqrt(timeslice_energy_corr_normed_d[timeslice_diff]);
    }









    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int t = 0; t < int(T_measurement/tau); t++) {
            for (int timeslice_diff = 2; timeslice_diff <= N4/2; timeslice_diff++) {
                timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t] = 0.0;
                for (int timeslice_base = N4/2 + 1 - timeslice_diff; timeslice_base < N4/2; timeslice_base++) {
                    timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t] +=
                        (timeslice_energy[thread_id][timeslice_base][t] - timeslice_energy_av)
                        *(timeslice_energy[thread_id][(timeslice_base + timeslice_diff)%N4][t] - timeslice_energy_av);
                }
                for (int timeslice_base = N4 + 1 - timeslice_diff; timeslice_base < N4; timeslice_base++) {
                    timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t] +=
                        (timeslice_energy[thread_id][timeslice_base][t] - timeslice_energy_av)
                        *(timeslice_energy[thread_id][(timeslice_base + timeslice_diff)%N4][t] - timeslice_energy_av);
                }

                timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t] /= (timeslice_diff - 1)*2;
                timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t] /= timeslice_energy_av*timeslice_energy_av;
            }
        }
    }

    for (int timeslice_diff = 2; timeslice_diff <= N4/2; timeslice_diff++) {
        timeslice_energy_corr_normed_av_old[timeslice_diff] = 0.0;
        for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
            for (int t = 0; t < int(T_measurement/tau); t++) {
                timeslice_energy_corr_normed_av_old[timeslice_diff] +=
                    timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t];
            }
        }
        timeslice_energy_corr_normed_av_old[timeslice_diff] /= threads_per_beta*T_measurement/tau;
    }

    for (int timeslice_diff = 2; timeslice_diff <= N4/2; timeslice_diff++) {
        timeslice_energy_corr_normed_d_old[timeslice_diff] = 0.0;
        for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
            for (int t = 0; t < int(T_measurement/tau); t++) {
                timeslice_energy_corr_normed_d_old[timeslice_diff] +=
                    (timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t] - timeslice_energy_corr_normed_av_old[timeslice_diff])
                    *(timeslice_energy_corr_normed_old[thread_id][timeslice_diff][t] - timeslice_energy_corr_normed_av_old[timeslice_diff]);
            }
        }
        timeslice_energy_corr_normed_d_old[timeslice_diff] /= (threads_per_beta*T_measurement/tau)*(threads_per_beta*T_measurement/tau - 1);
        timeslice_energy_corr_normed_d_old[timeslice_diff] = sqrt(timeslice_energy_corr_normed_d_old[timeslice_diff]);
    }



    out << beta << ',' << s_av << ',' << s_d << ',' << tau_corr_av*tau << '\n';
    out << beta << ',' << timeslice_energy_av << ',' << timeslice_energy_d
                    << ',' << timeslice_energy_d/timeslice_energy_av  << '\n';
    out.close();
    out.open(output_file_name, ios::out | ios::app);


    for (int timeslice_diff = 2; timeslice_diff <= N4/2; timeslice_diff++) {
            out << beta << ',' << tau_corr_av << ',' << timeslice_diff << ','
                << timeslice_energy_corr_normed_av_old[timeslice_diff] << ','
                << timeslice_energy_corr_normed_d_old[timeslice_diff] << '\n';
    }
    for (int timeslice_diff = 2; timeslice_diff <= N4/2; timeslice_diff++) {
            out << beta << ',' << tau_corr_av << ',' << timeslice_diff << ','
                << timeslice_energy_corr_normed_av[timeslice_diff] << ','
                << timeslice_energy_corr_normed_d[timeslice_diff] << '\n';
    }
    out.close();
    out.open(output_file_name, ios::out | ios::app);


    out.close();



}


//  half-spaces updates
void MeasureScalarGlueballForBeta_v2_0(unsigned int key, int iBeta) {
//    float beta = 0.01*iBeta;



    char output_file_name[500];
    sprintf(output_file_name, system_measurements_file_format, N1, N2, N3, N4, bc_code, matr_dim, key, iBeta);
    ofstream out(output_file_name, ios::out | ios::app);


    double **s = new double*[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        s[thread_id] = new double[int(T_measurement/tau) + 1];
    }

    double ***timeslice_energy = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy[thread_id] = new double*[N4];
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            timeslice_energy[thread_id][timeslice] = new double[int(T_measurement/tau) + 1];
        }
    }




    thread thread_array[threads_per_beta];

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        thread_array[thread_id] =
            thread(MeasureScalarGlueballForBeta_ReadingThreadFunction_v2, key,
                iBeta, thread_id, s[thread_id], timeslice_energy[thread_id]);
    }

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        thread_array[thread_id].join();
    }





    double s_av = 0.0;
    double s_d = 0.0;
    double tau_corr_av = 0.0;
    ActionAutocorrelation(s_av, s_d, tau_corr_av, s);

//    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
//        for (int t = 0; t < int(T_measurement/3/tau); t++) {
//            out << t << ',' << s_timecorr_normed[t][thread_id] << '\n';
//        }
//    }



    double timeslice_energy_av = 0.0;
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            for (int t = 0; t < int(T_measurement/tau); t++) {
                timeslice_energy_av += timeslice_energy[thread_id][timeslice][t];
            }
        }
    }
    timeslice_energy_av /= threads_per_beta*N4*T_measurement/tau;

    double timeslice_energy_d = 0.0;
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            for (int t = 0; t < int(T_measurement/tau); t++) {
                timeslice_energy_d += (timeslice_energy[thread_id][timeslice][t] - timeslice_energy_av)*
                                        (timeslice_energy[thread_id][timeslice][t] - timeslice_energy_av);
            }
        }
    }
    timeslice_energy_d /= (threads_per_beta*N4*T_measurement/tau)*(threads_per_beta*N4*T_measurement/tau);
    timeslice_energy_d = sqrt(timeslice_energy_d);




    double ***timeslice_energy_corr_product = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy_corr_product[thread_id] = new double*[N4/2 + 1];
        for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
            timeslice_energy_corr_product[thread_id][timeslice_diff] = new double[int(T_measurement/tau) + 1];
        }
    }
    double **timeslice_energy_corr_product_jack = new double*[N4/2 + 1];
    for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
        timeslice_energy_corr_product_jack[timeslice_diff] = new double[threads_per_beta];
    }
    double *timeslice_energy_jack = new double[threads_per_beta];
    double **timeslice_energy_corr_normed_jack = new double*[N4/2 + 1];
    for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
        timeslice_energy_corr_normed_jack[timeslice_diff] = new double[threads_per_beta];
    }
    double *timeslice_energy_corr_normed_av = new double[N4/2 + 1];
    double *timeslice_energy_corr_normed_d = new double[N4/2 + 1];
    double ***timeslice_energy_corr_normed_old = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy_corr_normed_old[thread_id] = new double*[N4/2 + 1];
        for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
            timeslice_energy_corr_normed_old[thread_id][timeslice_diff] = new double[int(T_measurement/tau) + 1];
        }
    }
    double *timeslice_energy_corr_normed_av_old = new double[N4/2 + 1];
    double *timeslice_energy_corr_normed_d_old = new double[N4/2 + 1];


    MeasureScalarGlueballForBeta_v2_processing (key, iBeta,
                s_av, s_d, tau_corr_av,
                timeslice_energy,
                timeslice_energy_av, timeslice_energy_d,
                timeslice_energy_corr_product,
                timeslice_energy_corr_product_jack,
                timeslice_energy_jack,
                timeslice_energy_corr_normed_jack,
                timeslice_energy_corr_normed_av,
                timeslice_energy_corr_normed_d,
                timeslice_energy_corr_normed_old,
                timeslice_energy_corr_normed_av_old,
                timeslice_energy_corr_normed_d_old);


    delete [] timeslice_energy_jack;
    delete [] timeslice_energy_corr_normed_av;
    delete [] timeslice_energy_corr_normed_d;
    delete [] timeslice_energy_corr_normed_av_old;
    delete [] timeslice_energy_corr_normed_d_old;



    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        delete [] s[thread_id];
    }
    delete [] s;



    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
            delete [] timeslice_energy_corr_product[thread_id][timeslice_diff];
        }
        delete [] timeslice_energy_corr_product[thread_id];
    }
    delete [] timeslice_energy_corr_product;

    for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
        delete [] timeslice_energy_corr_product_jack[timeslice_diff];
    }
    delete [] timeslice_energy_corr_product_jack;

    for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
        delete [] timeslice_energy_corr_normed_jack[timeslice_diff];
    }
    delete [] timeslice_energy_corr_normed_jack;

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
            delete [] timeslice_energy_corr_normed_old[thread_id][timeslice_diff];
        }
        delete [] timeslice_energy_corr_normed_old[thread_id];
    }
    delete [] timeslice_energy_corr_normed_old;



    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice = 0; timeslice < N4/2 + 1; timeslice++) {
            delete [] timeslice_energy[thread_id][timeslice];
        }
        delete [] timeslice_energy[thread_id];
    }
    delete [] timeslice_energy;
}


// Precalculated average loop, v2 version
void MeasureScalarGlueballForBeta_v2_1(unsigned int key, int iBeta) {
//    float beta = 0.01*iBeta;



    char output_file_name[500];
    sprintf(output_file_name, system_measurements_file_format, N1, N2, N3, N4, bc_code, matr_dim, key, iBeta);
    ofstream out(output_file_name, ios::out | ios::app);


    double **s = new double*[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        s[thread_id] = new double[int(T_measurement/tau) + 1];
    }

    double ***timeslice_energy = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy[thread_id] = new double*[N4];
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            timeslice_energy[thread_id][timeslice] = new double[int(T_measurement/tau) + 1];
        }
    }




    thread thread_array[threads_per_beta];

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        thread_array[thread_id] =
            thread(MeasureScalarGlueballForBeta_ReadingThreadFunction_v2, key,
                iBeta, thread_id, s[thread_id], timeslice_energy[thread_id]);
    }

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        thread_array[thread_id].join();
    }





    double s_av = 0.0;
    double s_d = 0.0;
    double tau_corr_av = 0.0;
    ActionAutocorrelation(s_av, s_d, tau_corr_av, s);

//    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
//        for (int t = 0; t < int(T_measurement/3/tau); t++) {
//            out << t << ',' << s_timecorr_normed[t][thread_id] << '\n';
//        }
//    }



    const double timeslice_energy_av = 1.7345867*3*N1*N2*N3;
    const double timeslice_energy_d = 0.01;





    double ***timeslice_energy_corr_product = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy_corr_product[thread_id] = new double*[N4/2 + 1];
        for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
            timeslice_energy_corr_product[thread_id][timeslice_diff] = new double[int(T_measurement/tau) + 1];
        }
    }
    double **timeslice_energy_corr_product_jack = new double*[N4/2 + 1];
    for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
        timeslice_energy_corr_product_jack[timeslice_diff] = new double[threads_per_beta];
    }
    double *timeslice_energy_jack = new double[threads_per_beta];
    double **timeslice_energy_corr_normed_jack = new double*[N4/2 + 1];
    for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
        timeslice_energy_corr_normed_jack[timeslice_diff] = new double[threads_per_beta];
    }
    double *timeslice_energy_corr_normed_av = new double[N4/2 + 1];
    double *timeslice_energy_corr_normed_d = new double[N4/2 + 1];
    double ***timeslice_energy_corr_normed_old = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy_corr_normed_old[thread_id] = new double*[N4/2 + 1];
        for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
            timeslice_energy_corr_normed_old[thread_id][timeslice_diff] = new double[int(T_measurement/tau) + 1];
        }
    }
    double *timeslice_energy_corr_normed_av_old = new double[N4/2 + 1];
    double *timeslice_energy_corr_normed_d_old = new double[N4/2 + 1];


    MeasureScalarGlueballForBeta_v2_processing (key, iBeta,
                s_av, s_d, tau_corr_av,
                timeslice_energy,
                timeslice_energy_av, timeslice_energy_d,
                timeslice_energy_corr_product,
                timeslice_energy_corr_product_jack,
                timeslice_energy_jack,
                timeslice_energy_corr_normed_jack,
                timeslice_energy_corr_normed_av,
                timeslice_energy_corr_normed_d,
                timeslice_energy_corr_normed_old,
                timeslice_energy_corr_normed_av_old,
                timeslice_energy_corr_normed_d_old);


    delete [] timeslice_energy_jack;
    delete [] timeslice_energy_corr_normed_av;
    delete [] timeslice_energy_corr_normed_d;
    delete [] timeslice_energy_corr_normed_av_old;
    delete [] timeslice_energy_corr_normed_d_old;


    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        delete [] s[thread_id];
    }
    delete [] s;


    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
            delete [] timeslice_energy_corr_product[thread_id][timeslice_diff];
        }
        delete [] timeslice_energy_corr_product[thread_id];
    }
    delete [] timeslice_energy_corr_product;

    for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
        delete [] timeslice_energy_corr_product_jack[timeslice_diff];
    }
    delete [] timeslice_energy_corr_product_jack;

    for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
        delete [] timeslice_energy_corr_normed_jack[timeslice_diff];
    }
    delete [] timeslice_energy_corr_normed_jack;

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
            delete [] timeslice_energy_corr_normed_old[thread_id][timeslice_diff];
        }
        delete [] timeslice_energy_corr_normed_old[thread_id];
    }
    delete [] timeslice_energy_corr_normed_old;



    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice = 0; timeslice < N4/2 + 1; timeslice++) {
            delete [] timeslice_energy[thread_id][timeslice];
        }
        delete [] timeslice_energy[thread_id];
    }
    delete [] timeslice_energy;
}



void MeasureScalarGlueballForBeta_v2_0_precalculated(unsigned int key, int iBeta) {
//    float beta = 0.01*iBeta;



    char output_file_name[500];
    sprintf(output_file_name, system_measurements_file_format, N1, N2, N3, N4, bc_code, matr_dim, key, iBeta);
    ofstream out(output_file_name, ios::out | ios::app);


    double **s = new double*[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        s[thread_id] = new double[int(T_measurement/tau) + 1];
    }

    double ***timeslice_energy = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy[thread_id] = new double*[N4];
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            timeslice_energy[thread_id][timeslice] = new double[int(T_measurement/tau) + 1];
        }
    }



    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int t = 0; t < int(T_measurement/tau); t++) {
            char file_name_out_measured[500];

            sprintf(file_name_out_measured,
            system_single_measurement_file_format,
                    iBeta, N1, N2, N3, N4, bc_code, matr_dim, thread_id, t);

            ifstream out_mes;

            out_mes.open(file_name_out_measured, ios::in | ios::binary);

            for (int timeslice = 0; timeslice < N4; timeslice++) {
                out_mes.read((char *) &timeslice_energy[thread_id][timeslice][t],
                                sizeof(double));
            }
            out_mes.read((char *) &s[thread_id][t], sizeof(double));

            out_mes.close();
        }
    }





    double s_av = 0.0;
    double s_d = 0.0;
    double tau_corr_av = 0.0;
    ActionAutocorrelation(s_av, s_d, tau_corr_av, s);

//    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
//        for (int t = 0; t < int(T_measurement/3/tau); t++) {
//            out << t << ',' << s_timecorr_normed[t][thread_id] << '\n';
//        }
//    }



    double timeslice_energy_av = 0.0;
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            for (int t = 0; t < int(T_measurement/tau); t++) {
                timeslice_energy_av += timeslice_energy[thread_id][timeslice][t];
            }
        }
    }
    timeslice_energy_av /= threads_per_beta*N4*T_measurement/tau;

    double timeslice_energy_d = 0.0;
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice = 0; timeslice < N4; timeslice++) {
            for (int t = 0; t < int(T_measurement/tau); t++) {
                timeslice_energy_d += (timeslice_energy[thread_id][timeslice][t] - timeslice_energy_av)*
                                        (timeslice_energy[thread_id][timeslice][t] - timeslice_energy_av);
            }
        }
    }
    timeslice_energy_d /= (threads_per_beta*N4*T_measurement/tau)*(threads_per_beta*N4*T_measurement/tau);
    timeslice_energy_d = sqrt(timeslice_energy_d);




    double ***timeslice_energy_corr_product = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy_corr_product[thread_id] = new double*[N4/2 + 1];
        for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
            timeslice_energy_corr_product[thread_id][timeslice_diff] = new double[int(T_measurement/tau) + 1];
        }
    }
    double **timeslice_energy_corr_product_jack = new double*[N4/2 + 1];
    for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
        timeslice_energy_corr_product_jack[timeslice_diff] = new double[threads_per_beta];
    }
    double *timeslice_energy_jack = new double[threads_per_beta];
    double **timeslice_energy_corr_normed_jack = new double*[N4/2 + 1];
    for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
        timeslice_energy_corr_normed_jack[timeslice_diff] = new double[threads_per_beta];
    }
    double *timeslice_energy_corr_normed_av = new double[N4/2 + 1];
    double *timeslice_energy_corr_normed_d = new double[N4/2 + 1];
    double ***timeslice_energy_corr_normed_old = new double**[threads_per_beta];
    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        timeslice_energy_corr_normed_old[thread_id] = new double*[N4/2 + 1];
        for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
            timeslice_energy_corr_normed_old[thread_id][timeslice_diff] = new double[int(T_measurement/tau) + 1];
        }
    }
    double *timeslice_energy_corr_normed_av_old = new double[N4/2 + 1];
    double *timeslice_energy_corr_normed_d_old = new double[N4/2 + 1];


    MeasureScalarGlueballForBeta_v2_processing (key, iBeta,
                s_av, s_d, tau_corr_av,
                timeslice_energy,
                timeslice_energy_av, timeslice_energy_d,
                timeslice_energy_corr_product,
                timeslice_energy_corr_product_jack,
                timeslice_energy_jack,
                timeslice_energy_corr_normed_jack,
                timeslice_energy_corr_normed_av,
                timeslice_energy_corr_normed_d,
                timeslice_energy_corr_normed_old,
                timeslice_energy_corr_normed_av_old,
                timeslice_energy_corr_normed_d_old);


    delete [] timeslice_energy_jack;
    delete [] timeslice_energy_corr_normed_av;
    delete [] timeslice_energy_corr_normed_d;
    delete [] timeslice_energy_corr_normed_av_old;
    delete [] timeslice_energy_corr_normed_d_old;


    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        delete [] s[thread_id];
    }
    delete [] s;


    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
            delete [] timeslice_energy_corr_product[thread_id][timeslice_diff];
        }
        delete [] timeslice_energy_corr_product[thread_id];
    }
    delete [] timeslice_energy_corr_product;

    for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
        delete [] timeslice_energy_corr_product_jack[timeslice_diff];
    }
    delete [] timeslice_energy_corr_product_jack;

    for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
        delete [] timeslice_energy_corr_normed_jack[timeslice_diff];
    }
    delete [] timeslice_energy_corr_normed_jack;

    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice_diff = 0; timeslice_diff < N4/2 + 1; timeslice_diff++) {
            delete [] timeslice_energy_corr_normed_old[thread_id][timeslice_diff];
        }
        delete [] timeslice_energy_corr_normed_old[thread_id];
    }
    delete [] timeslice_energy_corr_normed_old;



    for (int thread_id = 0; thread_id < threads_per_beta; thread_id++) {
        for (int timeslice = 0; timeslice < N4/2 + 1; timeslice++) {
            delete [] timeslice_energy[thread_id][timeslice];
        }
        delete [] timeslice_energy[thread_id];
    }
    delete [] timeslice_energy;
}



void Measurement(unsigned int key) {
    thread beta_threads[(iBeta_high - iBeta_low)/iBeta_step + 1];

    for (int t = 0; t <= (iBeta_high - iBeta_low)/iBeta_step; t++) {
        beta_threads[t] = thread(MeasureCreutzRatioForBeta, key, iBeta_low + t*iBeta_step);
    }
    for (int t = 0; t <= (iBeta_high - iBeta_low)/iBeta_step; t++) {
        beta_threads[t].join();
    }
}





int main() {
    chrono::time_point<chrono::system_clock>
                        start_time = chrono::system_clock::now();

    unsigned int key = start_time.time_since_epoch().count();







    Equilibration(key);

//    Measurement(key);









    chrono::time_point<chrono::system_clock>
                    end_time = chrono::system_clock::now();
    chrono::duration<double> diff_time = end_time - start_time;
    cout << "Execution time: " << diff_time.count() << " sec" << endl;

    return EXIT_SUCCESS;
}
