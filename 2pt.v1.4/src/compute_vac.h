#ifndef COMPUTE_VAC_H
#define COMPUTE_VAC_H

#include "statistics.h"
#include "compute.h"
#include <fftw3.h>

using namespace std;

class ComputeVAC : public Compute {
public:
    ComputeVAC();
    ~ComputeVAC();
    void init(MODEL *);
    void rd_frame(MODEL *,int);
    void report(ostream *, MODEL *);

private:
   void compute(ostream *,MODEL *);
   void calc_thermo_vals(ostream *,MODEL *);
   void allocate();
   void setup(MODEL *);
   void getDOF(MODEL *);
   void calc_vacv(MODEL *);
   void compute_vel_temp();
   void calc_pwr_spectrum();
   void calc_acf(int, int, MODEL *);
   void calc_acf_ang(int, int, MODEL *);
   void do_2pt(ostream *,MODEL *);
   void do_vib_analysis();
   void mem_usage(MODEL *);
   void dump_freq(MODEL *);
   void print_thermo(char *, int, int, ostream *, ostream *);
   void get_rot_temp(MODEL *,int,double *);
   double get_2pt_val(int,int,int);
   void fix_missing_grp_mol_entry(MODEL *);

   class Memory * memory;
   STAT ***Pvac;   //Pvac[ngrp][vtype][nstep]
   ATOM *vacatom,*tatom;
   CELL vaccell;
   PROPERTY vacprp;
    GROUP vacgrp;

   fftw_complex *in,*out;
   fftw_plan p1;
   fftw_plan p2;

   int natom; //number of atoms
   int trj_atom_eng;
    double **Cvc;	    //Cvc[ngrp][vtype] = constant volume heat capacity from fluc(H)/NkT^2
   double **vacT;  //average temperature  //vacT[ngrp][vtype]
   double **vacDF; //degrees of freedom   //vacDF[ngrp][vtype]
   double **vacE;  //average energy vacE[ngrp][vtype]
   double **pI;  //nmol*3*sizeof(double)
   float ****vacvv; //vacvv[vtype][step][atom][xyz]
   double ***thermo; //thermo[14][ngrp][vactype]
   double **f2pt, **K2pt, **hsdf;
   double *wep,*wsp,*wap,*wcvp,*wspd,*wsr,*war,*wer,*wcvr;
   int  *vacread; //atom read list
   int has_atom_eng;

   double vacdtime; //trajectory dump frequency
   double pwrfreq; //power spectrum frequency scaling factor
   double vacP;  //average pressure
   double *vacV;  //average volume for each group
   double trjT;  //average temperature from trj file
   double trdf; //translation and rotational degrees of freedom removed
   int  cmol; //consider molecules
   int  vactype; //types of velocities
   int  vacmaxf,vacnsteps; //correlation steps, number of origins
   int  c2pt; //2pt correction type
   int  vacnpass,vacipass; //number of passes for reading the whole trj
   int  vactatmpassed,vactmolpassed,vacatmpass,vacmolpass,vaclastmol; //number of atoms read each pass
   int  vacatmperpass,vacpassatms,vacatmend,vacatmpassed;
   int  nframe,tot_N,nused,ngrp,nmol;
};

#endif
