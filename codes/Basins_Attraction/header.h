/* -----------------------------------------------------------------
 * header.h
 * -----------------------------------------------------------------
 * This is the function file to compute Basins of attraction and 
 * others function support for the piezo-magneto-elastic beam
 * -----------------------------------------------------------------
 *  programmer: 
 *       João Peterson (ligier.peterson@gmail.com)
 *       Americo Cunha (americo.cunha@uerj.br)
 *
 *  last update: Jan, 2021
*/ -----------------------------------------------------------------
#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED

#include <iostream> // cout
#include <fstream>  // ofstream
#include <cmath>    // round
#include <vector>
#include <time.h>

using namespace std;

void percentage(int a, int b);

void linspace(double a, double b, int N, double output[]);

bool test(double **yss, int ly, double **AT, int colat, double tol);

void RKF4(double IC[], int s, double params[], double t0, double t1, double h,
              void (*func)(double, double[], double[], double[]), double **y);

void pvi(double t, double y[], double param[], double ydot[]);

void fill_AT(double **y, int RowYStart, int RowYEnd, double **AT, int ColATStart, int ColATEnd);

double z1test(double x[], int size_x);

double mean_abs_diff(double x[], int size_x);

void cumsum(double x[], int size_x, double result[]);

double sum(double a[], int size_a);

double mean(double a[], int size_a);

double corr(double x[], double y[], int size_n);

double median(double x[], int size_x);

void my_sort(double x[], int size_x);
#endif // HEADER_H_INCLUDED
