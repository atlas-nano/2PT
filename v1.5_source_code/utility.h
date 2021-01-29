/* Routines from Numerical Recipes (slightly modified). */
#ifndef UTILITY_H
#define UTILITY_H
#include <math.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include "constant.h"
#include "structure.h"

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

void jacobi(double **, int , double d[], double **);
int jacobi0(double **, int , double d[], double **);
void eigsrt(double d[], double **, int );
void moment(double *, int, float *, float *, float *, float *, float *, float *);
int eigsrt0(double d[], double **, int );
void ROTATE(double **,int ,int ,int ,int ,double *,double *);
double scweighting(double);
double sqweighting(double);
double polylog(double , double );
double search2PT(double );
double HSDF(double *,int ,double ,int ,double , double);
double TTDF(double *,int ,double );
void twoPT(double *,double *,double ,double ,double ,int ,double, double);
void HSweighting(double *,double *,double *,double *,double *,double *,double ,double ,double ,double , double ,double ,double *,double *,double *,double);
int twoPTmf(double *,double ,int ,double ,int ,double ,double &);
double cal_AB(double &, double &, double , double , double );
double Sgmf(double ,double , double ,int ,double );
void copyright(ostream *);
ele * loadelements(); 
char* toLowerCase(char* );
#endif
