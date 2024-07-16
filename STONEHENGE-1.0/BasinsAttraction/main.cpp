// -----------------------------------------------------------------
// main.cpp
// -----------------------------------------------------------------
// This is the main file for a program which trace the Basins of 
// attraction for the piezo-magneto-elastic beam
// -----------------------------------------------------------------
//  programmer: 
//       João Peterson (ligier.peterson@gmail.com)
//       Americo Cunha (americo.cunha@uerj.br)
//
//  last update: Jan, 2021
// -----------------------------------------------------------------

#include <iostream>
#include <thread>
#include <mutex>
#include <fstream>  // ofstream
#include <cmath>    // round
#include "header.h"

using namespace std;

void task_t2(double *d_params, int *i_params, double *param, double** AT, double** IC);
void task_t3(double *d_params, int *i_params, double *param, double** AT, double** IC);
void task_t4(double *d_params, int *i_params, double *param, double** AT, double** IC);

int main()
{
    cout << "// -------------------- //\n";
	cout << "// Basins of attraction\n";
	cout << "// -------------------- //\n" << endl;
	cout << "Please, wait...\n" << endl;

	// ------------------------------------------- //
    // Simulation parameters
    // ------------------------------------------- //
    double t0 = 0.0;            // initial time
    double t1 = 5e3;         // final time
    double h  = 5e-2;           // time step
    long int Ndt = round((t1-t0)/h); // number of time steps
    long int Nss = round(0.9*Ndt);   // Steady-state
    double tol = 0.01;          // tolerance for comparison test
    int Nat = 10;               // number of attractors

    double x0_min = -3.0;   // min. of initial displacement
    double x0_max = 3.0;    // max. of initial displacement
    double xdot0_min = -3.0; // min. of initial velocity
    double xdot0_max = 3.0;  // max. of initial velocity
    int N_x = 50;           // number of points
    int N_xdot = 50;        // number of points

    double *x0 = new double[N_x];       // vector of x0 values
    linspace(x0_min, x0_max, N_x, x0);

    double *xdot0 = new double[N_xdot]; // vector of xdot0 values
    linspace(xdot0_min, xdot0_max, N_xdot, xdot0);

    double v0 = 0.0;
    double ic[3];
    ic[2] = v0;

    double d_params [] = {t0, t1, h, tol, x0_min, x0_max, xdot0_min, xdot0_max, v0};
    int i_params[] = {Ndt, Nss, Nat, N_x, N_xdot};
    // ------------------------------------------- //
	// Physical parameters
	// ------------------------------------------- //
	double ksi, chi, f, omega, lambda, kappa;
	//double ksi, chi, f, lambda, kappa;
	ksi    = 0.01;      // mechanical damping ratio
	chi    = 0.05;      // dimensionless piezoelectric coupling term (mechanical)
	f      = 0.083;     // dimensionless excitation amplitude
	omega  = 0.5;       // dimensionless excitation frequency
	lambda = 0.05;      // dimensionless time constant reciprocal
	kappa  = 0.5;       // dimensionless piezoelectric coupling term (electrical)

	double param[] ={ksi, chi, f, omega, lambda, kappa};

	// ------------------------------------------- //
	// Allocating memory
	// ------------------------------------------- //
	//attractors output matrix
    /*
     This will create a matrix with dimensions [Ndt-Nss+1][3*Nat],
     the data is meant to be stored in that way:
     -- first attractor --   -- second attractor --
     x0[0]  xdot0[0]  v0[0]  x0[0]  xdot0[0]  v0[0] ...
     x0[1]  xdot0[1]  v0[1]  x0[1]  xdot0[1]  v0[1] ...
     ...     ...      ...    ...     ...      ...
    */
    double **AT = new double*[Ndt-Nss+1];
    for(int i=0; i<Ndt-Nss+1; i++){AT[i] = new double[3*Nat];}
    for(int i=0; i<Ndt-Nss+1; i++){
        for(int j=0; j<3*Nat; j++){AT[i][j] = 0;}
    }

    // IC output matrix
    /*
     This matrix is organized as follows:
     dimensions: [N_x][N_xdot]
     z_[0][0]  z_[0][1]  z_[0][2]
     z_[1][0]  z_[1][1]  z_[1][2]
     ...       ...       ...
    */
    double **IC = new double *[N_x];
    for(int i=0; i<N_x; i++){IC[i] = new double[N_xdot];}

    // response matrix (dimensions: [Ndt][3])
    double **y = new double*[Ndt];
    for(int i=0; i<Ndt; i++){y[i] = new double[3];}

    double *Qvolt = new double[(Ndt-Nss)/50];

    // ------------------------------------------- //
    // Integrating...
    // ------------------------------------------- //
    thread t2(task_t2, d_params, i_params, param, ref(AT), ref(IC));
    thread t3(task_t3, d_params, i_params, param, ref(AT), ref(IC));
    thread t4(task_t4, d_params, i_params, param, ref(AT), ref(IC));

    for(int nx=0; nx<N_x; nx=nx+4)
    {
        //percentage(nx, N_x);
        ic[0] = x0[nx];                                     // updating displacement
        for(int nxdot=0; nxdot<N_xdot; nxdot++)
        {
            ic[1] = xdot0[nxdot];                           // updating velocity

            RKF4(ic, 3, param, t0, t1, h, pvi, y);          // Integrate the dynamical system with RKF4

            double **yss = new double *[Ndt-Nss+1];         // steady-state response
            for(int i=0; i<Ndt-Nss+1; i++)
            {
                yss[i] = new double[3];
                yss[i][0] = y[Nss+i-1][0];
                yss[i][1] = y[Nss+i-1][1];
                yss[i][2] = y[Nss+i-1][2];
            }
            for(int i=0; i<(Ndt-Nss)/50; i++){
                Qvolt[i] = yss[(i*50)][2];
            }

            if(abs(z1test(Qvolt, (Ndt-Nss)/50)) < 0.2)          // if regular
            {
                for(int nat=0; nat<Nat; nat++)                  // cataloging attractors
                {
                    if(AT[0][nat*3] == 0)                       // if attractor is empty
                    {
                        fill_AT(yss, 0, Ndt-Nss, AT, nat*3, nat*3+2);   // fill attractor
                        IC[nx][nxdot] = nat+1;                          // fill IC matrix
                        break;
                    }
                    else if(test(yss, Ndt-Nss+1, AT, nat*3, tol))       // if same attractor
                    {
                        IC[nx][nxdot] = nat+1;                          // fill IC matrix
                        break;
                    }
                }
            }
            else                                                        // chaotic
            {
                IC[nx][nxdot] = 0;                                      // fill IC matrix
            }

            for(int i=0; i<Ndt-Nss+1; i++){delete [] yss[i];}
            delete [] yss;
            yss = 0;
        }
    }

    // Join threads
    t2.join();
    t3.join();
    t4.join();

    // save information
    ofstream ATfile("AT_[" + to_string(int(x0_min)) + "," + to_string(int(x0_max)) + "]_" +
                    "[" + to_string(int(xdot0_min)) + "," + to_string(int(xdot0_max)) + "]_" +
                    to_string(N_x) + "X" + to_string(N_xdot) + "_f" + to_string(int(f*1000)) +
                    "_O" + to_string(int(omega*100)) + ".dat", ios::trunc);
    for(int i=0; i<Ndt-Nss+1; i++){
        for(int j=0; j<3*Nat; j++){
            ATfile << AT[i][j] << " ";
        }
        ATfile << endl;
    }
    ATfile.close();

    cout << "ok";
    ofstream ICfile ("IC_[" + to_string(int(x0_min)) + "," + to_string(int(x0_max)) + "]_" +
                     "[" + to_string(int(xdot0_min)) + "," + to_string(int(xdot0_max)) + "]_" +
                    to_string(N_x) + "X" + to_string(N_xdot) + "_f" + to_string(int(f*1000)) +
                    "_O" + to_string(int(omega*100)) + ".dat", ios::trunc);
    for(int j=0; j<N_xdot; j++){
        for(int i=0; i<N_x; i++){
            ICfile << IC[i][N_xdot-1-j] << " ";
        }
        ICfile << endl;
    }
    ICfile.close();

    delete [] x0;
    x0 = 0;
    delete [] xdot0;
    xdot0 = 0;
    for(int i=0; i<Ndt; i++){delete [] y[i];}
    delete [] y;
    y = 0;

    for(int i=0; i<N_x; i++){delete [] IC[i];}
    delete [] IC;
    IC = 0;

    for(int i=0; i<Ndt-Nss+1; i++){delete [] AT[i];}
    delete [] AT;
    AT = 0;

    delete [] Qvolt;
    Qvolt = 0;

    return 0;
}

void task_t2(double *d_params, int *i_params, double *param, double** AT, double** IC){
    double t0 = d_params[0];
    double t1 = d_params[1];
    double h = d_params[2];
    double tol = d_params[3];
    double x0_min = d_params[4];
    double x0_max = d_params[5];
    double xdot0_min = d_params[6];
    double xdot0_max = d_params[7];
    double v0 = d_params[8];

    int Ndt = i_params[0];
    int Nss = i_params[1];
    int Nat = i_params[2];
    int N_x = i_params[3];
    int N_xdot = i_params[4];

    // vectors
    double *x0 = new double[N_x];       // vector of x0 values
    linspace(x0_min, x0_max, N_x, x0);

    double *xdot0 = new double[N_xdot]; // vector of xdot0 values
    linspace(xdot0_min, xdot0_max, N_xdot, xdot0);

    double ic[3];
    ic[2] = v0;

    // response matrix (dimensions: [Ndt][3])
    double **y = new double*[Ndt];
    for(int i=0; i<Ndt; i++){y[i] = new double[3];}

    double *Qvolt = new double[(Ndt-Nss)/50];


    // calculating
    for(int nx=1; nx<N_x; nx=nx+4)
    {
        //percentage(nx, N_x);
        ic[0] = x0[nx];                                     // updating displacement
        for(int nxdot=0; nxdot<N_xdot; nxdot++)
        {
            ic[1] = xdot0[nxdot];                           // updating velocity

            RKF4(ic, 3, param, t0, t1, h, pvi, y);     // Integrate the dynamical system with RKF4

            double **yss = new double *[Ndt-Nss+1];         // steady-state response
            for(int i=0; i<Ndt-Nss+1; i++)
            {
                yss[i] = new double[3];
                yss[i][0] = y[Nss+i-1][0];
                yss[i][1] = y[Nss+i-1][1];
                yss[i][2] = y[Nss+i-1][2];
            }
            for(int i=0; i<(Ndt-Nss)/50; i++){
                Qvolt[i] = yss[(i*50)][2];
            }

            if(abs(z1test(Qvolt, (Ndt-Nss)/50)) < 0.2)     // if regular
            {
                for(int nat=0; nat<Nat; nat++)                  // cataloging attractors
                {
                    if(AT[0][nat*3] == 0)                       // if attractor is empty
                    {
                        fill_AT(yss, 0, Ndt-Nss, AT, nat*3, nat*3+2);   // fill attractor
                        //IC[nx*N_xdot+nxdot][0] = ic[0];                 // fill IC matrix
                        //IC[nx*N_xdot+nxdot][1] = ic[1];
                        //IC[nx*N_xdot+nxdot][3] = nat+1;
                        IC[nx][nxdot] = nat+1;
                        break;
                    }
                    else if(test(yss, Ndt-Nss+1, AT, nat*3, tol))        // if same attractor
                    {
                        //IC[nx*N_xdot+nxdot][0] = ic[0];                 // fill IC matrix
                        //IC[nx*N_xdot+nxdot][1] = ic[1];
                        //IC[nx*N_xdot+nxdot][3] = nat+1;
                        IC[nx][nxdot] = nat+1;
                        break;
                    }
                }
            }
            else                                    // chaotic
            {
                //IC[nx*N_xdot+nxdot][0] = ic[0];     // fill IC matrix
                //IC[nx*N_xdot+nxdot][1] = ic[1];
                //IC[nx*N_xdot+nxdot][3] = 0;
                IC[nx][nxdot] = 0;
            }


            for(int i=0; i<Ndt-Nss+1; i++){delete [] yss[i];}
            delete [] yss;
            yss = 0;

        }
    }
    delete [] x0;
    x0 = 0;
    delete [] xdot0;
    xdot0 = 0;
    for(int i=0; i<Ndt; i++){delete [] y[i];}
    delete [] y;
    y = 0;
    delete [] Qvolt;
    Qvolt = 0;
}

void task_t3(double *d_params, int *i_params, double *param, double** AT, double** IC){
    double t0 = d_params[0];
    double t1 = d_params[1];
    double h = d_params[2];
    double tol = d_params[3];
    double x0_min = d_params[4];
    double x0_max = d_params[5];
    double xdot0_min = d_params[6];
    double xdot0_max = d_params[7];
    double v0 = d_params[8];

    int Ndt = i_params[0];
    int Nss = i_params[1];
    int Nat = i_params[2];
    int N_x = i_params[3];
    int N_xdot = i_params[4];

    // vectors
    double *x0 = new double[N_x];       // vector of x0 values
    linspace(x0_min, x0_max, N_x, x0);

    double *xdot0 = new double[N_xdot]; // vector of xdot0 values
    linspace(xdot0_min, xdot0_max, N_xdot, xdot0);

    double ic[3];
    ic[2] = v0;

    // response matrix (dimensions: [Ndt][3])
    double **y = new double*[Ndt];
    for(int i=0; i<Ndt; i++){y[i] = new double[3];}

    double *Qvolt = new double[(Ndt-Nss)/50];


    // calculating
    for(int nx=2; nx<N_x; nx=nx+4)
    {
        //percentage(nx, N_x);
        ic[0] = x0[nx];                                     // updating displacement
        for(int nxdot=0; nxdot<N_xdot; nxdot++)
        {
            ic[1] = xdot0[nxdot];                           // updating velocity

            RKF4(ic, 3, param, t0, t1, h, pvi, y);     // Integrate the dynamical system with RKF4

            double **yss = new double *[Ndt-Nss+1];         // steady-state response
            for(int i=0; i<Ndt-Nss+1; i++)
            {
                yss[i] = new double[3];
                yss[i][0] = y[Nss+i-1][0];
                yss[i][1] = y[Nss+i-1][1];
                yss[i][2] = y[Nss+i-1][2];
            }
            for(int i=0; i<(Ndt-Nss)/50; i++){
                Qvolt[i] = yss[(i*50)][2];
            }

            if(abs(z1test(Qvolt, (Ndt-Nss)/50)) < 0.2)     // if regular
            {
                for(int nat=0; nat<Nat; nat++)                  // cataloging attractors
                {
                    if(AT[0][nat*3] == 0)                       // if attractor is empty
                    {
                        fill_AT(yss, 0, Ndt-Nss, AT, nat*3, nat*3+2);   // fill attractor
                        //IC[nx*N_xdot+nxdot][0] = ic[0];                 // fill IC matrix
                        //IC[nx*N_xdot+nxdot][1] = ic[1];
                        //IC[nx*N_xdot+nxdot][3] = nat+1;
                        IC[nx][nxdot] = nat+1;
                        break;
                    }
                    else if(test(yss, Ndt-Nss+1, AT, nat*3, tol))        // if same attractor
                    {
                        //IC[nx*N_xdot+nxdot][0] = ic[0];                 // fill IC matrix
                        //IC[nx*N_xdot+nxdot][1] = ic[1];
                        //IC[nx*N_xdot+nxdot][3] = nat+1;
                        IC[nx][nxdot] = nat+1;
                        break;
                    }
                }
            }
            else                                    // chaotic
            {
                //IC[nx*N_xdot+nxdot][0] = ic[0];     // fill IC matrix
                //IC[nx*N_xdot+nxdot][1] = ic[1];
                //IC[nx*N_xdot+nxdot][3] = 0;
                IC[nx][nxdot] = 0;
            }


            for(int i=0; i<Ndt-Nss+1; i++){delete [] yss[i];}
            delete [] yss;
            yss = 0;

        }
    }
    delete [] x0;
    x0 = 0;
    delete [] xdot0;
    xdot0 = 0;
    for(int i=0; i<Ndt; i++){delete [] y[i];}
    delete [] y;
    y = 0;
    delete [] Qvolt;
    Qvolt = 0;
}

void task_t4(double *d_params, int *i_params, double *param, double** AT, double** IC){
    double t0 = d_params[0];
    double t1 = d_params[1];
    double h = d_params[2];
    double tol = d_params[3];
    double x0_min = d_params[4];
    double x0_max = d_params[5];
    double xdot0_min = d_params[6];
    double xdot0_max = d_params[7];
    double v0 = d_params[8];

    int Ndt = i_params[0];
    int Nss = i_params[1];
    int Nat = i_params[2];
    int N_x = i_params[3];
    int N_xdot = i_params[4];

    // vectors
    double *x0 = new double[N_x];       // vector of x0 values
    linspace(x0_min, x0_max, N_x, x0);

    double *xdot0 = new double[N_xdot]; // vector of xdot0 values
    linspace(xdot0_min, xdot0_max, N_xdot, xdot0);

    double ic[3];
    ic[2] = v0;

    // response matrix (dimensions: [Ndt][3])
    double **y = new double*[Ndt];
    for(int i=0; i<Ndt; i++){y[i] = new double[3];}

    double *Qvolt = new double[(Ndt-Nss)/50];


    // calculating
    for(int nx=3; nx<N_x; nx=nx+4)
    {
        //percentage(nx, N_x);
        ic[0] = x0[nx];                                     // updating displacement
        for(int nxdot=0; nxdot<N_xdot; nxdot++)
        {
            ic[1] = xdot0[nxdot];                           // updating velocity

            RKF4(ic, 3, param, t0, t1, h, pvi, y);     // Integrate the dynamical system with RKF4

            double **yss = new double *[Ndt-Nss+1];         // steady-state response
            for(int i=0; i<Ndt-Nss+1; i++)
            {
                yss[i] = new double[3];
                yss[i][0] = y[Nss+i-1][0];
                yss[i][1] = y[Nss+i-1][1];
                yss[i][2] = y[Nss+i-1][2];
            }
            for(int i=0; i<(Ndt-Nss)/50; i++){
                Qvolt[i] = yss[(i*50)][2];
            }

            if(abs(z1test(Qvolt, (Ndt-Nss)/50)) < 0.2)     // if regular
            {
                for(int nat=0; nat<Nat; nat++)                  // cataloging attractors
                {
                    if(AT[0][nat*3] == 0)                       // if attractor is empty
                    {
                        fill_AT(yss, 0, Ndt-Nss, AT, nat*3, nat*3+2);   // fill attractor
                        //IC[nx*N_xdot+nxdot][0] = ic[0];                 // fill IC matrix
                        //IC[nx*N_xdot+nxdot][1] = ic[1];
                        //IC[nx*N_xdot+nxdot][3] = nat+1;
                        IC[nx][nxdot] = nat+1;
                        break;
                    }
                    else if(test(yss, Ndt-Nss+1, AT, nat*3, tol))        // if same attractor
                    {
                        //IC[nx*N_xdot+nxdot][0] = ic[0];                 // fill IC matrix
                        //IC[nx*N_xdot+nxdot][1] = ic[1];
                        //IC[nx*N_xdot+nxdot][3] = nat+1;
                        IC[nx][nxdot] = nat+1;
                        break;
                    }
                }
            }
            else                                    // chaotic
            {
                //IC[nx*N_xdot+nxdot][0] = ic[0];     // fill IC matrix
                //IC[nx*N_xdot+nxdot][1] = ic[1];
                //IC[nx*N_xdot+nxdot][3] = 0;
                IC[nx][nxdot] = 0;
            }


            for(int i=0; i<Ndt-Nss+1; i++){delete [] yss[i];}
            delete [] yss;
            yss = 0;

        }
    }
    delete [] x0;
    x0 = 0;
    delete [] xdot0;
    xdot0 = 0;
    for(int i=0; i<Ndt; i++){delete [] y[i];}
    delete [] y;
    y = 0;
    delete [] Qvolt;
    Qvolt = 0;
}

