#ifndef MODEL_H
#define MODEL_H
#include <stdlib.h>
#include <string.h>
#include "constant.h"
#include "structure.h"
#include "property.h"
#include "trajectory.h"
#include "control.h"

using namespace std;

class MODEL
{
public:
   MODEL();
   ~MODEL();

   int nmol;    //number of molecules
   int natom;   //number of atoms
   int nbond;   //number of bonds
   int nangle;  //number of angles
   int ntor;    //number of torsions
   int ninv;    //number of inversions
   int nimp;    //number of improper torsions
   int ngrp;    //number of groups

   int periodic;  //periodicity
   CONTROL ctl;       //control file
   CELL cell;         //cell properties
   TRAJECTORY trj;    //trajectory
   PROPERTY prp;      //model properties

   char   name[1024]; //name of the model

   ATOM *atom;
   BOND *bond;
   MOLECULE *mol;
   GROUP *grp;

   int xyzfreq;  //coordinate dump frequency
   int velfreq;  //velocity dump frequency
   int engfreq;  //energy dump frequency
   int strfreq;  //stress dump frequency
   double timestep; //simulation time step in ps

   void  init_atom();
   void  init_bond();
   void  init_mol();
   void  init_grp();

   int rd_strt();
   int rd_lmpdata(char *);     //read lammps data file
   int rd_bgf(char *);         //read bgf file
   int rd_grp(char *);         //read grp file
   int rd_amberprmtop(char *); //read amber prmtop file

   int bond2mol();    //find molecules from connected bonds
   void cnt2valence(); //connect to valance
   void element2prp(int do_mass); 

   void cal_mass(); //calc total mass of the model
   void cal_grp_mass(); //calc mass of each group
   void ck_range(); //check atom range for nonperiodic system;
   void cal_vcomp(); //calc velocity components
   void cal_T();  //calculate temperature of model using atomic velocities
   void find_cm(); //find the center of mass position of the model
   void find_molcm(); //find the center of mass position of each molecule
   void find_molingrp(); //find molecules in each group

   //trajectory
   void init_trj();
};

#endif
