#ifndef CONTROL_H
#define CONTROL_H
#include <stdlib.h>
#include <unistd.h>
#include "constant.h"
#include <string.h>

using namespace std;

class CONTROL
{
public:
   CONTROL();
   ~CONTROL();

   //input flags
   int     in_strt_flag;  //structure flag
   int     in_trj_flag;   //input trj flag: c2trj or asctrj or lmptrj
   int     in_grp_flag;   //group file flag

   //general analysis flags
   int    ana_flag;       //analysis flag
   int    ana_vac_tstep_flag; //timestep flag

   //velocity autocorrelation for thermodynamics flags
   int    ana_vac_flag;   //velocity autocorrelation analysis flag
   int    ana_vac_vol_flag;   //velocity autocorrelation volume flag
   int    ana_vac_const_flag; 
   int    ana_vac_rotsym_flag; //rotational symmetry flag
   int    ana_vac_linear_flag; //linear atom flag
   int    ana_vac_press_flag;  //system avg. pressure flag
   int    ana_vac_temp_flag;  //system avg. temperature flag
   int    ana_vac_eng_flag;    //system avg. energy (KE+PE) flag
   int    ana_vac_dump_freq_flag; //trajectory dump frequency flag (steps)

   //trajectory input flag
   int    in_lmp_thermo_flag;   //lammps thermodynamics flag
   int    in_lmp_atomeng_flag;  //lammps thermodynamics flag
   int    in_coord_amber_flag;  //amber coordinates flag
   int    in_coord_charmm_flag; //charmm coordinates flag
   int    in_coord_xyz_flag;    //xyz coordinates flag;

   //values
   double ana_vac_tstep; //simulation timestep

   char   in_strt[1024];        //structure file
   trjItem * in_trj;            //trajectory file(s)
   char   in_grp[1024];         //group file
   char   in_lmp_thermo[1024];  //LAMMPS thermodynamics file
   char   in_lmp_atomeng[1024]; //LAMMPS thermodynamics file

   char   ana_out[1024]; //name of analysis output

   int    in_trj_n;       //number of input trj file
   int    in_coord_n;     //number of input trj file
   int    ana_iframe;
   int    ana_fframe;
   int    ana_sframe;


   int    ana_vac_2pt; //flag for 2pt calculations, 0:no, 1:yes, 2:2pt for molecules
   int	  ana_vac_show_2pt; // flag for showing 2pt breakup, 0:no (default), 1:yes
   int    ana_vac_dump_freq; //trajectory dump frequency (steps)
   double ana_vac_corlen;  //maximum correlation length (percentage of the total frame number) in VAC calc
   double ana_vac_mem; //memory allocation in Mega Bytes
   double ana_vac_vol;   //system volume
   double ana_vac_press; //system pressure
   double ana_vac_temp;  //temperature
   double ana_vac_void_vol;
   char   ana_vac_linear[512]; //flag for linear molecule
   char   ana_vac_const[512]; //constraint degree of freedom
   char   ana_vac_rotsym[512]; //rotational symmetry number
   double ana_vac_eng[2];   //MD energy

   int    rd_ctl(char *);
   void   init_trj(int);
   void   out_ctl(ostream *);
};

#endif
