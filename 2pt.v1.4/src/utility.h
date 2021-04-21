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

void jacobi(double **a, int n, double d[], double **v);
int jacobi0(double **a, int n, double d[], double **v);
void eigsrt(double d[], double **v, int n);
int eigsrt0(double d[], double **v, int n);
void ROTATE(double **a,int i,int j,int k,int l,double *tau,double *s);
int solve_cubic_eq(double a,double b,double c,double *c3rts);
double scweighting(double upper);
double sqweighting(double upper);
double polylog(double n, double z);
double search2PT(double K);
double HSDF(double *pwr,int lenth,double fmin,int nmol,double fract_f);
double TTDF(double *pwr,int lenth,double fmin);
void twoPT(double *tmpg,double *tmps,double s0,double sv,double v,int nmol,double fract_f);
void HSweighting(double *wep,double *wap,double *wcvp,double *wsehs,double *wsp,double *wspd,double y,double mass,double nmol,double tranT, double rotT,double volume,double *wsr,double *war,double *rT,double rs);
void copyright(ostream *);
ele * loadelements(); 
char* toLowerCase(char* );
#endif
