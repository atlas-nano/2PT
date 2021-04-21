#ifndef PROPERTY_H
#define PROPERTY_H
#include <stdlib.h>
#include "constant.h"

using namespace std;

class PROPERTY
{
public:
   PROPERTY () {
     mass=chg=dipole=mu[0]=mu[1]=mu[2]=0;
     peroid=step=0;
     time=T=P=V=rho=Et=Ek=Ep=fvol=0;
     Ebond=Eangle=Etorsion=Einversion=Evdw=Eel=Ehb=0;
     strs[0]=strs[1]=strs[2]=strs[3]=strs[4]=strs[5]=0;
     cmpv[0]=cmpv[1]=cmpv[2]=0;
   }
   ~PROPERTY () {
   }
   double mass;
   double chg;
   double dipole;
   double mu[3];
   double cmpv[3]; //center of mass position

   int peroid; 
   int step;  //md step
   double time; //md time in ps
   double T;  //temperature in K
   double P;  //pressure in GPa
   double V;  //volume in A3
   double fvol; //free volume in A3
   double rho; //density in g/cc
   double Et; //total energy in kcal/mol
   double Ek; //total kinetic energy in kcal/mol
   double Ep; //total potential energy in kcal/mol
   double Eke; //total electron kinetic energy in kcal/mol (from CPMD)
   double Ete; //total energy (Et+Eke) in kcal/mol (from CPMD)
   double Ebond; //bond
   double Eangle; //angle
   double Etorsion; //torsion
   double Einversion; //inversion
   double Evdw; //van der Waals
   double Eel;  //Coulomb
   double Ehb;  //hydrogen bond
   double strs[6]; //stress tensor  xx,yy,zz,xy,xz,yz

};

#endif
