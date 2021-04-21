#ifndef STRUCTURE_H
#define STRUCTURE_H
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "constant.h"
#include "statistics.h"

using namespace std;

struct ele
{
	string name;
	string symbol;
	double mass;
};

class CELL
{
public:
   CELL();
   ~CELL();

   double la;    //length of a
   double lb;    //length of b
   double lc;    //length of c
   double alpha;
   double beta;
   double gamma;
   double volume;
   double o[3];  //origin
   double cb[3]; //center of box
   double a[3];  //vector of a
   double b[3];  //vector of b
   double c[3];  //vector of c
   double H[3][3];
   double Hinv[3][3];

   int s2r2H(double *);
   int H2others();
   int labc2H();
   int cal_Hinv();
   int cal_volume();
};

class ATOM
{
public:
   ATOM () {
      mass=chg=eng_tmp=0.0;
      ffid=mol=id=-1;
      ncnt=0;
      int i;
      for(i=0;i<maxcnt;i++) { bod[i]=1; cnt[i]=-1; }
      faratm=-1;
      fardst=0;
      strcpy(name,"X");
      chain = 88;
      strcpy(fftype,"X_");
      radius=-1.0;
      crad=-1.0;
      shift[0]=shift[1]=shift[2]=0;
      grp=0;
   }
   ~ATOM () { }
   double f[FORCETYPE][3];  //net forces (bond,angle,torsion,inversion,vdw,vdw corr,elec,ewaldsum,net)
   double pv[6];    //position vector in angstroms, pv[4-6] stores pv of previous step, used for Verlet
   double vel[3];   //velocity in A/ps
   double vt[3];    //translational velocity
   double vr[3];    //rotational velocity
   double vv[3];    //vibrational velocity
   double strs[6];  //atomic stress component for surface tension
   double fv[3];     //atomic force from lammps
   double img[3];     //atomic image flag from lammps
   double shift[3]; //shift to base unit cell, used in atom remapping
   double mass;     //mass in g/mol
   double chg;      //charge in electrons
   int    ffid;     //ATOMTYPE
   char   fftype[6];//forcefild type
   int    mol;      //molecule id
   int    grp;      //group id
   int    id;       //atom id
   int    ncnt;     //number of atoms connected to it
   int    cnt[maxcnt];   //ids of atoms connected to it
   int    bod[maxcnt];   //bond order
   char   name[3];  //atom name
   int    chain;    //atom chain
   int    faratm;   //id of the farthermost connect atom
   int    fardst;   //number of atoms away from faratm
   double radius;   //atomic radius
   double crad;     //covalent radius for covalent bond calculation
   double eng_tmp;

   STAT eng;      //atom energy from lammps
private:
};

class MOLECULE
{
public:
   MOLECULE();
   ~MOLECULE();

   int natom;
   int head;
   int tail;
   int linear;      //linear molecule
   int *atm;
   double mass;     //mass
   double chg;      //charge
   double **inertia;  //moment of inertia tensor
   double pv[3];    //position vector (center of mass) in angstroms
   double vel[3];   //center of mass velocity
   double mu[3];    //dipole
   double omega[3]; //angular velocity
   double angmom[3];//angular momentum
   double anguv[3]; //angular velocity
   double pI[3];    //principle moment of inertia (kg*m2)

   STAT eng;      //energy
   ATOM   **atom;
   class Memory * memory;

   int init_atom();
   int cal_mass();   //calculate center of mass 
   int cal_chgmu();  //calculate charge and dipole
   int cal_vc();    //calculate velocity components
   int find_cm();    //calculate center of mass pv[3]
};

class BOND
{
public:
   BOND();
   ~BOND();

   int    id;       //bond id
   int    ffid;
   int    lbd;      //bond when folded into one cell (for proton transfer)
   int atm[2];    //atom pairs
   double len0;     //equilibrium bond length angstroms
   double len;      //current bond length angstroms
   double *ffp;

   ATOM   *atom[2];

   int    cal_len();
};

class GROUP
{
public:
  GROUP();
  ~GROUP();

  int natom;
  int nmol;
  int constraint;  //number of constraints
  int rotsym;      //rotational symmetry
  int linear;      //linear molecule
  int color;
  int *atm;
  int *mol;
  double atomsize; 
  double bondsize;
  double labelsize;
  double chg;
  double mass;
  double mu[3];
  double vol;

    STAT eng;      //energy
  int init_atom();
};

#endif
