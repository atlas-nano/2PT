/*********************************************************************************
*                     Two-Phase Thermodynamics (2PT) Program                     *
*                          Shiang-Tai Lin (stlin@ntw.edu.tw)                     *
* Department of Chemical Engineering, National Tiawan University, Taipei, Taiwan *
*                     Prabal K. Maiti (maiti@physics.iisc.ernet.in)              *
*  Department of Physics, Indian Institute of Science, Bangalore, India, 560012  *
*                          Tod A Pascal (tpascal@wag.caltech.edu)                *
*                                     and                                        *
*                        William A Goddard III (wag@wag.caltech.edu)             *
*         Materials and Process Simulation Center, Caltech, Pasadena, CA USA     *
*                                Copyright (c) 2010                              *
*                                All rights Reserved                             *
*                             *************************                          *
*                                      Cite:                                     *
* S.T. Lin, M. Blanco and W.A. Goddard, J. Chem. Phys., 2003, 119, 11792-11805   *
* S.T. Lin, P.K. Maiti and W.A. Goddard, J. Phys. Chem. B., 2010, 114, 8191-8198 *
* T.A. Pascal, S.T. Lin and W.A. Goddard, PCCP, 2011, 13(1), 169-181             *
***********************************************************************************/

//parse input file
#include <iostream>
#include <fstream>
#include <sstream>
#include "constant.h"
#include "control.h"
#include "memory.h"

/* ---------------------------------------------------------------------- */
CONTROL::~CONTROL () 
{
  if(in_trj_flag) delete [] in_trj;
}

/* ---------------------------------------------------------------------- */
CONTROL::CONTROL () 
{
  in_strt_flag=0;
  in_trj_flag=0;
  in_grp_flag=0;

  ana_flag=0;

  ana_iframe=1;
  ana_fframe=0;
  ana_sframe=1;

  strcpy(in_grp,"none");
  in_trj=NULL;
  in_trj_n=0;
  in_coord_n=0;

  ana_vac_flag=1;
  ana_vac_corlen=0.5;
  ana_vac_mem=500; //memory allocation in Mega Bytes
  ana_vac_2pt=1; // 2pt calculations 
				 // 0:1PT for monoatomic systems
				 // 1:2PT for monoatomic systems (original of Lin,Blanco,Goddard JCP 2003) 
				 //	2:2PT-MF for monotonic systems (Sandia group)
				 //	3:1PT for molecules 
				 //	4:2PT for molecules (Lin,Maiti and Goddard, JPC 2012)
				 //	5:2PT-MF for molecules
  ana_vac_method = 1;
  ana_vac_mol_flag = 0;
  ana_vac_show_2pt = 0; //show 2pt breakup

  strcpy(ana_vac_const,"0"); //group constratined degrees of freedom
  ana_vac_const_flag = 0;
  strcpy(ana_vac_rotsym,"2"); //moleculue rotational symmtery
  ana_vac_rotsym_flag = 0;
  strcpy(ana_vac_linear,"0"); //linear molecule flag
  ana_vac_linear_flag = 0;
  ana_vac_eng[0]=ana_vac_eng[1]=0.0;
  ana_vac_eng_flag = 0;
  ana_vac_tstep = 0.001; //simulation timestep in ps
  ana_vac_tstep_flag=0; 
  ana_vac_vol = 0; //system volume in A^3
  ana_vac_vol_flag = 0; 
  ana_vac_press = 0; //system pressure in GPa
  ana_vac_press_flag = 0;
  ana_vac_temp = 0.0; //system temperature in K
  ana_vac_vsfactor = 1.0; //velocity scaling factor
  ana_vac_temp_flag = 0;
  ana_vac_void_vol=0.0;
  ana_vac_dump_freq_flag = 0;
  ana_vac_dump_freq = 0;

  in_lmp_thermo_flag=0;
  in_lmp_atomeng_flag=0;
  in_coord_amber_flag=0;
  in_coord_charmm_flag=0;
  in_coord_xyz_flag=0;

  has_grp_file_eng = 0;
  ana_vac_check_grp_eng = 1;
  ana_vac_per_mol_flag = 0;

  strcpy(ana_out,"2pt.default");
}

/* ---------------------------------------------------------------------- */
int CONTROL::rd_ctl(char *ctlf)
{
  int i;
  char opt[1024], null[1024]; 
  istringstream ss;
  ifstream inctl(ctlf,ios::in);

  cout<<"Reading control file "<<ctlf<<"...";
  if(!inctl.is_open()) return filenotopen;

  inctl>>opt;
  while(!inctl.eof()) {
    if(opt[0]=='#'||opt[0]=='*'||opt[0]=='!'||opt[0]=='/') { inctl.getline(opt,1024); }

    else if (strcmp(opt,"IN_LMPDATA")==0) {inctl>>in_strt; in_strt_flag=strtlmp; }
    else if (strcmp(opt,"IN_AMBERPRMTOP")==0) {inctl>>in_strt; in_strt_flag=strtamber; }
    else if (strcmp(opt,"IN_BGF")==0) {inctl>>in_strt; in_strt_flag=strtbgf; }

    else if (strcmp(opt,"IN_GROUPFILE")==0) { inctl>>in_grp; in_grp_flag=1; }

    else if (strcmp(opt,"IN_C2TRJ")==0     ||
	     strcmp(opt,"IN_ASCTRJ")==0    ||
	     strcmp(opt,"IN_USRTRJ")==0    ||
	     strcmp(opt,"IN_CPMDTRJ")==0   ||
	     strcmp(opt,"IN_LMPTRJ")==0    ||
	     strcmp(opt,"IN_AMBERTRJ")==0  ||
	     strcmp(opt,"IN_CHARMMTRJ")==0 ||
	     strcmp(opt,"IN_XYZTRJ")==0)   { 
      if (strcmp(opt,"IN_C2TRJ")==0) in_trj_flag=c2trj; 
      else if (strcmp(opt,"IN_ASCTRJ")==0)    in_trj_flag=asctrj;
      else if (strcmp(opt,"IN_USRTRJ")==0)    in_trj_flag=usrtrj;
      else if (strcmp(opt,"IN_CPMDTRJ")==0)   in_trj_flag=cpmdtrj;
      else if (strcmp(opt,"IN_LMPTRJ")==0)    in_trj_flag=lmptrj;
      else if (strcmp(opt,"IN_AMBERTRJ")==0)  in_trj_flag=ambertrj;
      else if (strcmp(opt,"IN_CHARMMTRJ")==0) in_trj_flag=charmmtrj;
      else if (strcmp(opt,"IN_XYZTRJ")==0) in_trj_flag=xyztrj; 
      inctl.getline(opt,1024);
      //determine the number of trjs
      ss.clear();
      ss.str(opt);
      in_trj_n=0;
      while (ss>>null) in_trj_n++;
      if(in_trj == NULL) { init_trj(in_trj_n); in_coord_n = in_trj_n; }
      else if (in_coord_n < in_trj_n) in_trj_n = in_coord_n;
      //now read in each
      ss.clear();
      ss.str(opt);
      i=0;
      while (ss>>in_trj[i].vel) {
	in_trj[i].molopt = ana_vac_2pt;
	i++;
      }
    } 

    //LAMMPS Options
    else if (strcmp(opt,"IN_LMP_THERMO_FILE")==0) { inctl>>in_lmp_thermo; in_lmp_thermo_flag=1; }
    else if (strcmp(opt,"IN_LMP_ATOMENG_FILE")==0) { inctl>>in_lmp_atomeng; in_lmp_atomeng_flag=1; }

    //AMBER/CHARMM/XYZ Options
    else if (strcmp(opt,"IN_AMBER_COORD_FILE")==0   ||
	     strcmp(opt,"IN_CHARMM_COORD_FILE")==0  ||
	     strcmp(opt,"IN_XYZ_COORD_FILE")==0) {
      if(strcmp(opt,"IN_AMBER_COORD_FILE")==0) in_coord_amber_flag=1;
      else if(strcmp(opt,"IN_CHARMM_COORD_FILE")==0) in_coord_charmm_flag=1;
      else in_coord_xyz_flag=1; 

      inctl.getline(opt,1024);
      //determine the number of trjs
      ss.clear();
      ss.str(opt);
      in_coord_n=0;
      while (ss>>null) in_coord_n++;
      if(in_trj == NULL) init_trj(in_coord_n); 
      if (in_coord_n > in_trj_n) in_coord_n = in_trj_n;
      //now read in each
      ss.clear();
      ss.str(opt);
      i=0;
      while (ss>>in_trj[i].coord) i++;
    }
    else if (strcmp(opt,"OUT_PREFIX")==0) { inctl>>ana_out; }
    else if (strcmp(opt,"ANALYSIS_OUT")==0) {inctl>>ana_out; } //backward compatibility
    //Trajectory selection
    else if (strcmp(opt,"ANALYSIS_FRAME_INITIAL")==0) { inctl>>ana_iframe; }
    else if (strcmp(opt,"ANALYSIS_FRAME_FINAL")==0) { inctl>>ana_fframe;  }
    else if (strcmp(opt,"ANALYSIS_FRAME_STEP")==0) { inctl>>ana_sframe;  }

    //Some settings
    else if (strcmp(opt,"ANALYSIS_PER_MOLECULE")==0) { inctl>>ana_vac_per_mol_flag; ana_flag=ana_vac_flag=1; }
    else if (strcmp(opt,"ANALYSIS_CHECK_GRP_ENG")==0) { inctl>>ana_vac_check_grp_eng; ana_flag=ana_vac_flag=1; }
    else if (strcmp(opt,"ANALYSIS_VAC_CORLENGTH")==0) { inctl>>ana_vac_corlen; ana_flag=ana_vac_flag=1; }
    else if (strcmp(opt,"ANALYSIS_VAC_MEMORYMB")==0) { inctl>>ana_vac_mem; }
    else if (strcmp(opt,"ANALYSIS_SHOW2PT")==0) { inctl>>ana_vac_show_2pt; }
    else if (strcmp(opt,"ANALYSIS_MOLECULE_FLAG")==0) { inctl>>ana_vac_mol_flag; if(ana_vac_mol_flag>1||ana_vac_mol_flag<0) ana_vac_mol_flag=0; }
    else if (strcmp(opt,"ANALYSIS_2PT_METHOD")==0) { inctl>>ana_vac_method; if(ana_vac_method>2||ana_vac_method<0) ana_vac_method=0; }
    else if (strcmp(opt,"ANALYSIS_VAC_2PT")==0) { 
		inctl>>ana_vac_2pt; //backward compatibility 
		if(ana_vac_2pt>5||ana_vac_2pt<0) ana_vac_2pt = 0;
		ana_vac_method = ana_vac_2pt;
		if(ana_vac_2pt>2) {ana_vac_mol_flag = 1; ana_vac_method -= 3; }
	}
    else if (strcmp(opt,"MD_FIXED_DF")==0) { ana_vac_const_flag = 1; inctl>>ana_vac_const; } 
    else if (strcmp(opt,"ANALYSIS_VAC_FIXED_DF")==0) { ana_vac_const_flag = 1; inctl>>ana_vac_const; } //backward compatalibity
    else if (strcmp(opt,"MOL_ROTN_SYMMETRY")==0) { ana_vac_rotsym_flag=ana_vac_flag = 1; inctl>>ana_vac_rotsym; }
    else if (strcmp(opt,"ANALYSIS_VAC_ROTN_SYMMETRY")==0) { ana_vac_rotsym_flag=ana_vac_flag = 1; inctl>>ana_vac_rotsym; } //backward compatibility
    else if (strcmp(opt,"MOL_LINEAR_FLAG")==0) { ana_vac_linear_flag=ana_vac_flag=1; inctl>>ana_vac_linear; }
    else if (strcmp(opt,"ANALYSIS_VAC_LINEAR_MOL")==0) { ana_vac_linear_flag=ana_vac_flag=1; inctl>>ana_vac_linear; } //backwards compatibility

    //Thermodynamics input
    else if (strcmp(opt,"MD_AVGVOLUME")==0) { inctl>>ana_vac_vol; ana_vac_flag=ana_vac_vol_flag= 1; }
    else if (strcmp(opt,"VOID_VOLUME")==0) { ana_vac_flag=1; inctl>>ana_vac_void_vol; }
    else if (strcmp(opt,"MD_TSTEP")==0) { inctl>>ana_vac_tstep; ana_vac_tstep_flag=1;}
    else if (strcmp(opt,"ANALYSIS_LMP_TSTEP")==0) { inctl>>ana_vac_tstep; ana_vac_tstep_flag=1;} //backward compatibility
    else if (strcmp(opt,"MD_AVGENERGY")==0) { 
	ana_vac_eng_flag=2;
	inctl.getline(opt,1024);
	ss.clear();
	ss.str(opt);
	i=0;
	while (ss>>null) {
	    if(i==2) {
		cout<<" Error in ana_roient "<<opt<<endl;
		return unregkey;
	    }
	    if(strcmp(null,"G")==0||strcmp(null,"g")==0) ana_vac_eng_flag=1;
	    else ana_vac_eng[i] = atof(null);
	    i++;
	}
    }
    else if (strcmp(opt,"MD_AVGPRESSURE")==0) { inctl>>ana_vac_press; ana_vac_flag=ana_vac_press_flag=1; }
    else if (strcmp(opt,"MD_AVGTEMPERATURE")==0) { inctl>>ana_vac_temp; ana_vac_flag=ana_vac_temp_flag=1; }
    else if (strcmp(opt,"TRAJ_DUMPFREQ")==0) { inctl>>ana_vac_dump_freq; ana_vac_flag=ana_vac_dump_freq_flag=1; }
	else if (strcmp(opt,"TRAJ_VEL_SFACTOR")==0) { 
		inctl.getline(opt,1024);
		ss.clear();
		ss.str(opt);
		ss>>null;
		if(strcmp(null,"CP2K")==0) ana_vac_vsfactor = 2188491.52E-2;
		else if (strcmp(null,"NAMD")==0) ana_vac_vsfactor = 20.45482706;
		else if (strcmp(null,"LAMMPS")==0) ana_vac_vsfactor = 1000;
		else ana_vac_vsfactor = atof(null);
	}
    else {
        cout<<" Error in ana_roient "<<opt<<endl;
        return unregkey;
    }
    inctl>>opt;
  }

  //error checking
  if( ! in_strt_flag && ana_vac_mol_flag) {
    cout<<" Error: Must specify structure file when using molecular 2PT option"<<endl;
    return filenotopen;
  } else if(! in_trj_n) {
    cout<<" Error: No trajectory file specified"<<endl;
    return filenotopen;
  } else if (in_trj_flag==ambertrj && (ana_vac_mol_flag && ! in_coord_amber_flag)) {
    cout<<" Error: Need to specify the coordinates for AMBER trajectory when using molecular option"<<endl;
    return filenotopen;
  } else if (in_trj_flag==charmmtrj && (ana_vac_mol_flag && ! in_coord_charmm_flag)) {
    cout<<" Error: Need to specify the coordinates for CHARMM trajectory when using molecular option"<<endl;
    return filenotopen;
  } else if (in_trj_flag==xyztrj && (ana_vac_mol_flag && ! in_coord_xyz_flag)) {
    cout<<" Error: Need to specify the coordinates for XYZ trajectory when using molecular option"<<endl;
    return filenotopen;
  } else if ((in_trj_flag != lmptrj && in_trj_flag != c2trj) && ! ana_vac_temp_flag) {
    cout<<" Error: Need to specify the system temperature"<<endl;
    return filenotopen;
  } else if (in_coord_n != in_trj_n) {
    cout<<" Error: Mismatch between number of trajectories with velocities "<<in_trj_n<<" and coordinates "<<in_coord_n<<endl;
    return countmismatch;
  }

  return none;
}

/* ---------------------------------------------------------------------- */
void CONTROL::out_ctl(ostream *outf)
{
  char null[1024];
  int i;

  if(in_strt_flag==strtlmp) {
    sprintf(null,"%-40s %s","IN_LMPDATA",in_strt);
    *outf<<null<<endl;
  } else if(in_strt_flag==strtbgf) {
    sprintf(null,"%-40s %s","IN_BGF",in_strt);
    *outf<<null<<endl;
  } else if(in_strt_flag==strtamber) {
    sprintf(null,"%-40s %s","IN_AMBER_PRMTOP",in_strt);
    *outf<<null<<endl;
  }

  if(in_trj_flag==c2trj) {
    sprintf(null,"%-40s","IN_C2TRJ");
  } else if(in_trj_flag==asctrj) {
    sprintf(null,"%-40s","IN_ASCTRJ");
  } else if(in_trj_flag==cpmdtrj) {
    sprintf(null,"%-40s","IN_CPMDTRJ");
  } else if(in_trj_flag==ambertrj) {
    sprintf(null,"%-40s","IN_USRTRJ");
  } else if(in_trj_flag==ambertrj) {
    sprintf(null,"%-40s","IN_AMBER_VEL");
  } else if(in_trj_flag==xyztrj) {
    sprintf(null,"%-40s","IN_XYZ_VEL");
  } else if(in_trj_flag==charmmtrj) {
    sprintf(null,"%-40s","IN_CHARMM_VEL");
  } else if(in_trj_flag==lmptrj) {
    for(i=0;i<in_trj_n;i++) {
      if(! in_lmp_thermo_flag) {
        strcpy(in_trj[i].thermo,in_trj[i].vel);
        strcat(in_trj[i].thermo,".eng");
      }
      if(! in_lmp_atomeng_flag) {
	strcpy(in_trj[i].atom_eng,in_trj[i].vel);
	strcat(in_trj[i].atom_eng,".atom.eng");
      }
      strcat(in_trj[i].vel,".lammps");
    }
    sprintf(null,"%-40s","IN_LAMMPS_TRJ");
  }
  for(i=0;i<in_trj_n;i++) sprintf(null,"%s %s",null,in_trj[i].vel);
  *outf<<null<<endl;

  if(in_trj_flag==lmptrj) {
    sprintf(null,"%-40s","IN_LAMMPS_THERMO");
    for(i=0;i<in_trj_n;i++) sprintf(null,"%s %s",null,in_trj[i].thermo);
    *outf<<null<<endl;
    sprintf(null,"%-40s","IN_LAMMPS_ATOMENG");
    for(i=0;i<in_trj_n;i++) sprintf(null,"%s %s",null,in_trj[i].atom_eng);
    *outf<<null<<endl;
  }else if (in_trj_flag != c2trj && ana_vac_mol_flag) {
    if(in_trj_flag==ambertrj) sprintf(null,"%-40s","IN_AMBER_COORD");
    else if(in_trj_flag==charmmtrj) sprintf(null,"%-40s","IN_CHARMM_COORD");
    else if(in_trj_flag==xyztrj) sprintf(null,"%-40s","IN_XYZ_COORD");
    for(i=0;i<in_trj_n;i++) sprintf(null,"%s %s",null,in_trj[i].coord);
    *outf<<null<<endl;
  }

  if(in_grp_flag) {
    sprintf(null,"%-40s %s","IN_GROUPFILE",in_grp);
    *outf<<null<<endl;
  }

  sprintf(null,"%-40s %s","OUT_PREFIX",ana_out);
  *outf<<null<<endl;

  sprintf(null,"%-40s %d","ANALYSIS_FRAME_INITIAL",ana_iframe);
  *outf<<null<<endl;
  sprintf(null,"%-40s %d","ANALYSIS_FRAME_FINAL",ana_fframe);
  *outf<<null<<endl;
  sprintf(null,"%-40s %d","ANALYSIS_FRAME_STEP",ana_sframe);
  *outf<<null<<endl;
  sprintf(null,"%-40s %f","ANALYSIS_VAC_CORLENGTH",ana_vac_corlen);
  *outf<<null<<endl;
  sprintf(null,"%-40s %f","ANALYSIS_VAC_MEMORYMB",ana_vac_mem);
  *outf<<null<<endl;
  sprintf(null,"%-40s %d","ANALYSIS_MOLECULE_OPTION",ana_vac_mol_flag);
  *outf<<null<<endl;
  if(ana_vac_method==0) sprintf(null,"%-40s %s","ANALYSIS_2PT_TYPE","1PT");
  else if(ana_vac_method==1) sprintf(null,"%-40s %s","ANALYSIS_2PT_TYPE","2PT");
  else if(ana_vac_method==2) sprintf(null,"%-40s %s","ANALYSIS_2PT_TYPE","2PT-MF");
  *outf<<null<<endl;
  ana_vac_2pt = ana_vac_method;
  if(ana_vac_mol_flag) ana_vac_2pt += 3;
  if(ana_vac_rotsym_flag) {
    sprintf(null,"%-40s %s","ANALYSIS_VAC_ROTN_SYMMETRY",ana_vac_rotsym);
    *outf<<null<<endl;
  }
  if(ana_vac_const_flag==1) {
    sprintf(null,"%-40s %s","ANALYSIS_VAC_FIXED_DF",ana_vac_const);
    *outf<<null<<endl;
  }
  if(ana_vac_linear_flag) {
    sprintf(null,"%-40s %s","ANALYSIS_VAC_LINEAR_MOL",ana_vac_linear);
    *outf<<null<<endl;
  }
  if(ana_vac_vol_flag) {
    sprintf(null,"%-40s %g A^3","MD_CELL_VOLUME",ana_vac_vol);
    *outf<<null<<endl;
  }
  if(ana_vac_void_vol > 0.0) {
    sprintf(null,"%-40s %g (A^3)","VOID_VOLUME",ana_vac_void_vol);
    *outf<<null<<endl;
  }
  if(ana_vac_tstep_flag) {
    sprintf(null,"%-40s %g (ps)","MD_TIMESTEP",ana_vac_tstep);
    *outf<<null<<endl;
  }
    if(ana_vac_eng_flag==1) {
	sprintf(null,"%-40s Reading from Group file","MD_AVG_STRAIN_ENERGY");
	*outf<<null<<endl;
    } else if(ana_vac_eng_flag==2) {
	sprintf(null,"%-40s %lf +/- %lf (kcal/mol)","MD_AVG_STRAIN_ENERGY",ana_vac_eng[0],ana_vac_eng[1]);
	*outf<<null<<endl;
    }
  if(ana_vac_press_flag) {
    sprintf(null,"%-40s %g (GPa)","MD_AVG_PRESSURE",ana_vac_press);
    *outf<<null<<endl;
  }
  if(ana_vac_temp_flag) {
    sprintf(null,"%-40s %g (K)","MD_AVG_TEMPERATURE",ana_vac_temp);
    *outf<<null<<endl;
  }
  if(ana_vac_dump_freq_flag) {
    sprintf(null,"%-40s %i (steps)","TRAJ_DUMP_FREQUENCY",ana_vac_dump_freq);
    *outf<<null<<endl;
  }
  if(ana_vac_check_grp_eng) {
    sprintf(null,"%-40s %d","CHECK_GROUP_ENERGY",ana_vac_check_grp_eng);
    *outf<<null<<endl;
  }
  if(ana_vac_per_mol_flag) {
    sprintf(null,"%-40s %d","THERMO_PER_MOLECULE",ana_vac_per_mol_flag);
    *outf<<null<<endl;
  }
  if(ana_vac_vsfactor != 0.0) {
	sprintf(null,"%-40s %f","TRJ_VELOCITY_SCALING_FACTOR",ana_vac_vsfactor);  
    *outf<<null<<endl;
  }
}

/* ---------------------------------------------------------------------- */

void CONTROL::init_trj(int trj_n)
{

  in_trj = new trjItem [trj_n];

/*  int i;
  for(i=0;i<trj_n;i++) {
    strcpy(in_trj[i].vel,"");
    strcpy(in_trj[i].coord,"");
    strcpy(in_trj[i].thermo,"");
    strcpy(in_trj[i].atom_eng,"");
    strcpy(in_trj[i].stress,"");
   } */ 
}
