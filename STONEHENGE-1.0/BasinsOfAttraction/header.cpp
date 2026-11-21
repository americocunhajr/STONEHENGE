// -----------------------------------------------------------------
// header.cpp
// -----------------------------------------------------------------
// This is the function file to compute Basins of attraction and 
// others function support for the piezo-magneto-elastic beam
// -----------------------------------------------------------------
//  programmer: 
//       João Peterson (ligier.peterson@gmail.com)
//        Americo Cunha (americo.cunha@uerj.br)
//
//  last update: Jan, 2021
// -----------------------------------------------------------------

#include <iostream> // cout
#include <fstream>  // ofstream
#include <cmath>    // round
#include <time.h>
#include <vector>
#include <algorithm>
#include "header.h"

using namespace std;

void percentage(int a, int b){
    float n = (float)(a+1) / (float)(b+1);
    float percents[10] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    for(int i=0; i<10; i++){
        if(n == percents[i]){
            cout << "\r";
            cout << "[" << (n)*100 << "%]";
       }
    }
}

void linspace(double a, double b, int N, double output[]){
    for(int i=0; i<N; i++){output[i] = a + i*(b-a)/(N-1);}
}

bool test(double **yss, int ly, double **AT, int colat, double tol){
	bool result = true;
	double ATmax, ATmin, ymax, ymin;
	for(int i=0; i<3; i++) // changing dim 2
	{
        ATmax = AT[0][colat+i];
        ATmin = AT[0][colat+i];
        ymax  = yss[0][i];
        ymin  = yss[0][i];

        for(int j=1; j<ly; j++)
        {
            ATmax = max(ATmax, AT[j][i+colat]);
            ATmin = min(ATmin, AT[j][i+colat]);
            ymax  = max(ymax, yss[j][i]);
            ymin  = min(ymin, yss[j][i]);
        }

		if(abs(ATmin-ymin) > tol || abs(ATmax-ymax) > tol)
		{
			result = false;
		}
	}

	return result;
}

void RKF4(double IC[], int s, double params[], double t0, double t1, double h,
              void (*func)(double, double[], double[], double[]), double **y)
{
    int Ndt = round((t1-t0)/h);   // number of time steps
    double *k1, *k2, *k3, *k4;    // Runge-Kutta's k
    double *y1, *y2, *y3, *y4;
    k1 = new double[s];
    k2 = new double[s];
    k3 = new double[s];
    k4 = new double[s];
    y1 = new double[s];
    y2 = new double[s];
    y3 = new double[s];
    y4 = new double[s];

    double t;

    t = t0;  // time
    for(int i=0; i<s; i++){y[0][i] = IC[i];}   // initial conditions

    for(int i=0; i<(Ndt-1); i++)
    {
            for(int m=0; m<s; m++){y1[m] = y[i][m];}
        (*func)(t, y1, params, k1);
            for(int m=0; m<s; m++){y2[m] = y[i][m] + 0.5*h*k1[m];}
        (*func)(t+0.5*h, y2, params, k2);
            for(int m=0; m<s; m++){y3[m] = y[i][m] + 0.5*h*k2[m];}
        (*func)(t+0.5*h, y3, params, k3);
            for(int m=0; m<s; m++){y4[m] = y[i][m] + h*k3[m];}
        (*func)(t+h, y4, params, k4);

        t += h;
        for(int m=0; m<s; m++){
            y[i+1][m] = y[i][m] + h*(0.1666666666666667)*(k1[m] + 2.0*k2[m] + 2.0*k3[m] + k4[m]);
        }
    }

    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] y1;
    delete[] y2;
    delete[] y3;
    delete[] y4;
}

void pvi(double t, double y[], double param[], double ydot[])
{
    // parameters
    /*
    double ksi    = param[0];
    double chi    = param[1];
    double f      = param[2];
    double omega  = param[3];
    double lambda = param[4];
    double kappa  = param[5];
    */
    //
    // system of equations
    ydot[0] = y[1];
    ydot[1] = -2.0*param[0]*y[1] + 0.5*y[0]*(1.0 - y[0]*y[0]) + param[1]*y[2] + param[2]*cos(param[3]*t);
    //ydot[1] = -2.0*param[0]*y[1] - y[0] + param[1]*y[2] + param[2]*cos(param[3]*t);
    ydot[2] = -param[4]*y[2] - param[5]*y[1];
}
// ----------------------------------------------------------------------------
void fill_AT(double **y, int RowYStart, int RowYEnd, double **AT, int ColATStart, int ColATEnd)
{
    for(int j=ColATStart; j<ColATEnd+1; j++){
        for(int i=RowYStart; i<RowYEnd+1; i++){
            AT[i-RowYStart][j] = y[i][j-ColATStart];
        }
    }
}
// ----------------------------------------------------------------------------
double z1test(double x[], int size_x)
{
    int j[size_x];
    for(int i=0; i<size_x; i++){j[i] = i+1;}
    int roundN10 = round(size_x/10);
    double t[roundN10];
    for(int i=0; i<roundN10; i++){t[i] = i+1;}
    double M[roundN10];

    int ns = 300; // number of samples
    double c[ns];
    for(int i=0; i<ns; i++){c[i] = ((double) rand()/(RAND_MAX))*4*atan(1.0);}

    double kcorr[ns];

    double *p, *q;
    p = new double[size_x];
    q = new double[size_x];

    double xcos[size_x], xsin[size_x];

    for(int its=0; its<ns; its++){
        for(int i=0; i<size_x; i++){
            xcos[i] = x[i]*cos(j[i]*c[its]);
            xsin[i] = x[i]*sin(j[i]*c[its]);
        }

        cumsum(xcos, size_x, p);
        cumsum(xsin, size_x, q);

        for(int n=0; n<roundN10; n++){
            double aux[size_x-n+1];
            for(int i=0; i<(size_x-n+1); i++){
                aux[i] = pow((p[i+n+1]-p[i]), 2.0) + pow((q[i+n+1]-q[i]), 2.0);
            }
            M[n] = mean(aux, size_x-n+1) - pow(mean(x, size_x), 2.0)*(1.0-cos(n*c[its]))/(1.0-cos(c[its]));
        }

        kcorr[its] = corr(t, M, roundN10);
    }

    vector<double> kcorr_less;
    vector<double> kcorr_more;

    for(int i=0; i<ns; i++)
    {
        if(c[i] < mean(c, ns))
        {
            kcorr_less.push_back(kcorr[i]);
        }
        else
        {
            kcorr_more.push_back(kcorr[i]);
        }
    }

    double array_k_less[kcorr_less.size()];
    double array_k_more[kcorr_more.size()];

    // checking for oversampling

    /*
    if( ( (max(x) - min(x))/(mean_abs_diff(x, size_x)) ) > 10 ||
       median(array_k_less, int(kcorr_less.size())) - median(array_k_more, int(kcorr_more.size())) > 0.5)
    {
        cout << "data is probably oversampled." << endl;
        cout << "Use coarser sampling or reduce the maximum value of c." << endl;
    }
    */


    if( ( (*max_element(x, x+size_x) - *min_element(x, x+size_x))/(mean_abs_diff(x, size_x)) ) > 10 ||
       median(array_k_less, int(kcorr_less.size())) - median(array_k_more, int(kcorr_more.size())) > 0.5)
    {
        cout << "data is probably oversampled." << endl;
        cout << "Use coarser sampling or reduce the maximum value of c." << endl;
    }



    delete [] p;
    delete [] q;

    return median(kcorr, ns);
}
// ----------------------------------------------------------------------------
double mean_abs_diff(double x[], int size_x)
{
    double y[size_x-1];

    for(int i=0; i<size_x-1; i++)
    {
        y[i] = abs(x[i+1] - x[i]);
    }

    return mean(y, size_x-1);
}
// ----------------------------------------------------------------------------
void cumsum(double x[], int size_x, double result[])
{
    result[0] = x[0];
    for(int i=1; i<size_x; i++){
        result[i] = x[i] + result[i-1];
    }
}
// ----------------------------------------------------------------------------
double sum(double a[], int size_a)
{
	double s = 0;
	for (int i = 0; i < size_a; i++)
	{
		s += a[i];
	}
	return s;
}
// ----------------------------------------------------------------------------
double mean(double a[], int size_a)
{
	return sum(a, size_a) / size_a;
}
// ----------------------------------------------------------------------------
double corr(double x[], double y[], int size_n)
{
    /* pearson coefficient

                    _n
                    \   (x_i - mean(x))(y_i - mean(y))
                    /_1
      r =  ------------------------------------------------------------
            sqrt(sum((x_i - mean(x))^2)) * sqrt(sum((y_i - mean(y))^2))

    */

    double N = 0;
    double D1 = 0;
    double D2 = 0;
    double mean_x = mean(x, size_n);
    double mean_y = mean(y, size_n);

    for(int i=0; i<size_n; i++){
        N += (x[i] - mean_x)*(y[i] - mean_y);
        D1 += pow((x[i] - mean_x), 2.0);
        D2 += pow((y[i] - mean_y), 2.0);
    }
    D1 = sqrt(D1);
    D2 = sqrt(D2);

    return (N/(D1*D2));
}
// ----------------------------------------------------------------------------
double median(double x[], int size_x)
{
	double result;

	my_sort(x, size_x);

	if(size_x % 2 == 0)
		result = (x[(size_x-1)/2] + x[(size_x+1)/2])/2;
	else
		result = x[size_x/2];

	return result;
}
// ----------------------------------------------------------------------------
void my_sort(double x[], int size_x)
{
    double a;
    for(int i=0; i<size_x; i++){
        for(int j=i; j<size_x; j++){
            if(x[j]<x[i]){
                a = x[j];
                x[j] = x[i];
                x[i] = a;
            }
        }
    }
}
