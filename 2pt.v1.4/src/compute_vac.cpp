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

//does the actual work
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include "constant.h"
#include "compute_vac.h"
#include "memory.h"
#include "utility.h"
#include "statistics.h"
//#include <omp.h>

/* ---------------------------------------------------------------------- */
ComputeVAC::~ComputeVAC()
{
	int i,j;

	for(i=0;i<ngrp;i++) {
		for(j=0;j<vactype;j++) delete [] Pvac[i][j];
		delete [] Pvac[i];
	}
	delete [] Pvac;

	delete [] vacatom;
	delete [] tatom;

	delete [] vacread;

	memory->destroy_2d_double_array(vacE);
	memory->destroy_2d_double_array(vacT);
	memory->destroy_2d_double_array(vacDF);
	memory->destroy_2d_double_array(Cvc);
	if (cmol)
		memory->destroy_2d_double_array(pI);
		
	delete [] wep;
	delete [] wsp;
	delete [] wap;
	delete [] wcvp;
	delete [] wspd;
	delete [] wsr;
	delete [] war;
	delete [] wer;
	delete [] wcvr;
	delete [] vacV;
	memory->destroy_2d_double_array(f2pt);
	memory->destroy_2d_double_array(K2pt);
	memory->destroy_2d_double_array(hsdf);
	memory->destroy_3d_double_array(thermo);

	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);
	fftw_free(in);
	fftw_free(out);
	void fftw_cleanup_threads(void);
}

/* ---------------------------------------------------------------------- */
ComputeVAC::ComputeVAC()
{

	Pvac = NULL;
	vacatom = NULL;
	tatom = NULL;
	vacread = NULL;

	vacE = NULL;
	vacT = NULL;
	vacDF = NULL;
	pI = NULL;

	wep = NULL;
	wsp = NULL;
	wap = NULL;
	wcvp = NULL;
	wspd = NULL;
	wsr = NULL;
	war = NULL;
	wer = NULL;
	wcvr = NULL;

	f2pt = NULL;
	K2pt = NULL;
	hsdf = NULL;
	thermo = NULL;

	Cvc = NULL;

	vacipass=0;
	vactatmpassed=0;	//total atoms passed
	vacatmpass=0;		 //# atom in current pass
	vacmolpass=0;		 //# molecule in current pass
	vactmolpassed=0;	//total molecules passed
	vaclastmol=0;		 //last molecule id from previous pass
	vacatmperpass = 0;

	vacV=NULL;
	vacP=0;
	nframe = 0;
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::init(MODEL *model)
{
	int i;
	char null[1024];

	ngrp = model->ngrp;
	nmol = model->nmol;
	natom = model->natom;

	nframe = (model->ctl.ana_fframe - model->ctl.ana_iframe + 1)/model->ctl.ana_sframe;

	cmol=model->ctl.ana_vac_2pt;

	switch(cmol) {
		case 0: cout<<" Mode analysis with 2PT corrections will be performed"<<endl; break;
		case 1: cout<<" Mode analysis with molecular and 2PT corrections will be performed"<<endl; break;
	}
	if(!cmol) vactype=1; 
	else vactype=VELTYPE+1; 
	model->find_molingrp();
	fix_missing_grp_mol_entry(model);

	//set some variables
	vacnsteps=(int)((model->ctl.ana_fframe-model->ctl.ana_iframe+1.0)/(model->ctl.ana_sframe));
	vacmaxf=(int)((model->ctl.ana_vac_corlen)*(vacnsteps-1))+1;
	tot_N=vacnsteps+vacmaxf;
	nused = tot_N/2;

	//calc the memory needed to read in coordinates of an atom from the whole trj
	mem_usage(model);

	vacatom=new ATOM [natom];
	for(i=0;i<natom;i++) {
		vacatom[i].mass=model->atom[i].mass;
		vacatom[i].mol=model->atom[i].mol;
		vacatom[i].id=model->atom[i].id;
		vacatom[i].grp=model->atom[i].grp;
	}

	dump_freq(model);

	//print some info
	cout<<" VAC Parms"<<endl;
	sprintf(null,"  %-10s : %i\n  %-10s : %i\n  %-10s : %i\n  %-10s : %i\n  %-10s : %i\n",
			"start",model->ctl.ana_iframe,"end",model->ctl.ana_fframe,"step",model->ctl.ana_sframe,
			"nsteps",vacnsteps,"corr_steps",vacmaxf);
	cout<<null;
	sprintf(null,"  %-10s : %g (ps)\n  %-10s : %g (ps)\n  %-10s : %g (ps)\n  %-10s : %g (ps)\n",
			"dump_freq",vacdtime,"tstep",model->ctl.ana_vac_tstep,"traj_len",vacdtime*vacnsteps,"corr_len",vacdtime*vacmaxf);
	cout<<null;

	allocate();
	setup(model);
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::dump_freq(MODEL *model)
{
	int i,j;

	tatom=model->atom;
	model->atom=vacatom;
	for(i=0;i<model->trj.ntrjf;i++) {
		model->trj.strj[i].atom=model->atom;
		model->trj.strj[i].cell=&vaccell;
		model->trj.strj[i].prp =&vacprp;
	}

	if(!model->ctl.ana_vac_dump_freq_flag) {
		model->trj.rd_frame(model->ctl.ana_iframe+2);
		vacdtime=vacprp.time;
		model->trj.rd_frame(model->ctl.ana_iframe+model->ctl.ana_sframe+2);
		vacdtime=vacprp.time-vacdtime;
	} else {
		vacdtime = model->ctl.ana_vac_dump_freq * model->ctl.ana_vac_tstep;
	}

	model->atom=tatom;
	trj_atom_eng = 0;
	for(i=0;i<model->trj.ntrjf;i++) {
		model->trj.strj[i].atom=model->atom;
		for(j=0;j<natom;j++)
			if(model->trj.strj[i].header.eng_flag != 0.0) 
				trj_atom_eng = has_atom_eng = 1;
		model->trj.strj[i].cell=&(model->cell);
		model->trj.strj[i].prp =&(model->prp);
	}
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::mem_usage(MODEL *model)
{
	int i;

	//memory usage statistics
	//Pvac: ngrp*vactype*tot_N*sizeof(STAT)
	//vacT: ngrp*vactype*sizeof(double)
	//vacDF: ngrp*vactype*sizeof(double)
	//vacatom: natom*sizeof(ATOM)
	//pI: nmol*3*sizeof(double)
	//in: tot_N*sizeof(fftw_complex)
	//out: tot_N*sizeof(fftw_complex)
	//pwr: 2*ngrp*tot_N/2*sizeof(double);
	
	double memusage= ngrp*vactype*( tot_N*sizeof(STAT) +sizeof(double)*2 )
				 + natom*sizeof(ATOM)
				 + model->nmol*3*sizeof(double)
				 + 2*tot_N*sizeof(fftw_complex)
				 + ngrp*( tot_N*sizeof(double) );
	memusage/=MEGABYTE;
	//vacread: vacatmperpass*sizeof(int)
	//vacvv: vactype*vacnsteps*vacatmperpass*3*sizeof(float)

	//calc the memory needed to read in coordinates of an atom from the whole trj
	double memperatm=vactype*vacnsteps*sizeof(float)*3.0/MEGABYTE;
	vacatmperpass=(int)( (model->ctl.ana_vac_mem-memusage)/memperatm);

	if(vacatmperpass<1) {
		cout<<" Allocated memory ("<<model->ctl.ana_vac_mem<<" MB) insufficient for the calculation"<<endl;
		cout<<" Memory needed for reading one atom: "<<memperatm<<" MB"<<endl;
		cout<<" Memory needed for additional arrays: "<<memusage<<" MB"<<endl;
		cout<<" Minimum memory needed is "<<memusage+memperatm<<" MB"<<endl;
		cout<<" Please increase ANALYSIS_VAC_MEMORYMB"<<endl;
		cout<<vacatmperpass<<endl;
		exit(0);
	}
	cout<<" VAC extra memory needed: "<<memusage<<" MB"<<endl;
	cout<<"  Pvac[ngrp][vactype][tot_N]: "<<ngrp*vactype*tot_N*sizeof(STAT)/MEGABYTE<<endl;
	cout<<"  vacT[ngrp][vactype]       : "<<ngrp*vactype*sizeof(double)/MEGABYTE<<endl;
	cout<<"  vacDF[ngrp][vactype]      : "<<ngrp*vactype*sizeof(double)/MEGABYTE<<endl;
	cout<<"  vacatom[natom]            : "<<natom*sizeof(ATOM)/MEGABYTE<<endl;
	cout<<"  pI[nmol][3]               : "<<model->nmol*3*sizeof(double)/MEGABYTE<<endl;
	cout<<"  fftw in out               : "<<2*tot_N*sizeof(fftw_complex)/MEGABYTE<<endl;
	cout<<"  pwr[2][ngrp][tot_N/2]     : "<<ngrp*tot_N*sizeof(double)/MEGABYTE<<endl;
	cout<<"  ngroup "<<ngrp<<" vactype "<<vactype<<" tot_N "<<tot_N<<" STAT "<<sizeof(STAT)/MEGABYTE<<endl;
	cout<<" VAC memory needed for one atom: "<<memperatm<<" MB"<<endl;
	cout<<" VAC memory needed for all atoms: "<<natom*memperatm<<" MB"<<endl;
	vacnpass=(int)(natom/vacatmperpass+1.0);
	vacatmperpass=(vacatmperpass>natom?natom:vacatmperpass);
	vacnpass=(int)(natom/vacatmperpass)+(natom%vacatmperpass==0?0:1);
	printf(" ALLOC MEM %5.1f MB, Max Atom per Pass %d, nPass %d\n",model->ctl.ana_vac_mem,vacatmperpass,vacnpass);

	for(i=0;i<model->nmol;i++) {
		if(model->mol[i].natom>vacatmperpass) {
			cout<<" Error: Max allowed atoms per Pass "<<vacatmperpass<<" is less than the number of atoms "<<model->mol[i].natom<<" contained in molecule "<<i+1<<endl;
			cout<<" Please increase ANALYSIS_VAC_MEMORYMB"<<endl;
			exit(0);
		}
	}
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::allocate()
{
	int i,j,k;

	//fftw_init_threads();
	//fftw_plan_with_nthreads(omp_get_max_threads());

	vacread= new int [vacatmperpass]; //id of atoms read in each pass
	vacvv= memory->create_4d_float_array(vactype,vacnsteps,vacatmperpass,3,"ComputeVAC::vacvv");

	vacE= memory->create_2d_double_array(ngrp, vactype,"ComputeVAC::vacE");
	vacT= memory->create_2d_double_array(ngrp, vactype,"ComputeVAC::vacT");
	vacDF= memory->create_2d_double_array(ngrp, vactype,"ComputeVAC::vacDF");
	if (cmol) pI = memory->create_2d_double_array(nmol, 3,"ComputeVAC::vacpI");

	k=ngrp;
	wep =	new double [k];
	wsp =	new double [k];
	wap =	new double [k];
	wcvp = new double [k];
	wspd = new double [k];
	wsr = new double [k];
	war = new double [k];
	wer = new double [k];
	wcvr = new double [k];
	
	f2pt = memory->create_2d_double_array(2,k,"ComputeVAC::vacf2pt");
	K2pt = memory->create_2d_double_array(2,k,"ComputeVAC::vacK2pt");
	hsdf = memory->create_2d_double_array(2,k,"ComputeVAC::vachsdf");
		Cvc = memory->create_2d_double_array(ngrp,vactype,"ComputeVAC::Cvc");
	thermo = memory->create_3d_double_array(14,ngrp,vactype,"ComputeVAC::thermo");

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tot_N);
	out= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tot_N);

	p1 = fftw_plan_dft_1d(tot_N,in,out,FFTW_BACKWARD,FFTW_ESTIMATE); //backward(+1)
	p2 = fftw_plan_dft_1d(tot_N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);	//forward(-1)

	Pvac=new STAT **[ngrp];
	for(i=0;i<ngrp;i++) {
		Pvac[i]=new STAT *[vactype];
		for(j=0;j<vactype;j++) Pvac[i][j]=new STAT [tot_N];
	}

	for(i=0;i<14;i++)
		for(j=0;j<ngrp;j++)
			for(k=0;k<vactype;k++) {
				thermo[i][j][k] = 0.0;
				vacE[j][k] = 0.0;
			}

	 vacV = new double [ngrp];
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::setup(MODEL *model)
{
	int i,j,k,tmpi;

	if(cmol==0) {
		if(vacatmperpass>natom) j=natom;
		else j=vacatmperpass;
		for(i=0;i<j;i++) vacread[vacatmpass++]=i;
	} else {
		tmpi=0;
		for(i=vaclastmol;i<model->nmol;i++) {
			tmpi += model->mol[i].natom;
			if(tmpi<=vacatmperpass) {
				for(j=0;j<model->mol[i].natom;j++) vacread[vacatmpass++]=model->mol[i].atm[j];
				vacmolpass++;
			} else break; //break i-loop
		}
		vactmolpassed+=vacmolpass;
		vaclastmol=i;
	}
	vactatmpassed +=vacatmpass;

	//account for translation and rotational degrees of freedom removed
	if(model->ctl.in_trj_flag==usrtrj) model->periodic=1;
	if (!model->periodic){
		if(natom==1) trdf=3;
		else if(natom==2) trdf=5;
		else trdf=6;
	} else trdf = 3;
	trjT=0;

	//now specify if molecule in each group is linear
	if(model->ctl.ana_vac_linear_flag && (model->ctl.ana_vac_linear[0]=='G'||model->ctl.ana_vac_linear[0]=='g')) ;
	else for(i=0;i<ngrp;i++) model->grp[i].linear = atoi(model->ctl.ana_vac_linear);

	//determine the degree of freedom in each group
	getDOF(model);

	//now specify the rotational symmetry for each group
	if(model->ctl.ana_vac_rotsym_flag && (model->ctl.ana_vac_rotsym[0]=='G'||model->ctl.ana_vac_rotsym[0]=='g')) ;
	else for(i=0;i<ngrp;i++) model->grp[i].rotsym = atoi(model->ctl.ana_vac_rotsym);

	if(cmol) 
		for(i=0;i<model->nmol;i++) 
			for(k=0;k<3;k++) 
				pI[i][k]+=model->mol[i].pI[k];

	cout<<flush;

	//initialize the vacvv array
	if(cmol==0) 
		for(i=0;i<nframe;i++)
			for(j=0;j<vacatmpass;j++)
				for(k=0;k<3;k++) vacvv[0][i][j][k]= 0.0;
	else 
		for(i=0;i<nframe;i++)
			for(j=0;j<vacatmpass;j++)
				for(k=0;k<3;k++) {
					vacvv[vtrans][i][j][k]=0.0;
					vacvv[vrotat][i][j][k]=0.0;
					vacvv[vimvib][i][j][k]=0.0;
					vacvv[vangul][i][j][k]=0.0;
					vacvv[vtotal][i][j][k]=0.0;
				}

	//energy
	if(model->ctl.ana_vac_eng_flag==2) {
		model->grp[ngrp-1].eng.sum = model->grp[ngrp-1].eng.avg = model->ctl.ana_vac_eng[0];
		model->grp[ngrp-1].eng.std = model->ctl.ana_vac_eng[1];
		model->grp[ngrp-1].eng.sigma2x = model->grp[ngrp-1].eng.std*model->grp[ngrp-1].eng.std;
		model->grp[ngrp-1].eng.nXY =	model->grp[ngrp-1].eng.npt = 1;
	}

}

/* ---------------------------------------------------------------------- */
void ComputeVAC::getDOF(MODEL *model)
{
	int i,j,tmpi,tmpj,has_mols;

	for(i=0;i<ngrp;i++) {
		//cout<<"grp "<<i<<" nmol "<<model->grp[i].nmol<<endl;;
		for(j=0;j<vactype;j++) vacT[i][j]=0;
		if(cmol==0) vacDF[i][vactype-1]=3.0*model->grp[i].natom;
		else {
			vacDF[i][vtrans]=0;
			vacDF[i][vrotat]=0;
			vacDF[i][vangul]=0;
			vacDF[i][vimvib]=0;
			vacDF[i][vtotal]=0;
			has_mols = 1;
			if(! model->grp[i].nmol) { has_mols = 0; model->grp[i].nmol = 1;} //hack so that can assign dof even if using part of molecule
			for(j=0;j<model->grp[i].nmol;j++) {
				if(has_mols)
					tmpi=model->mol[ model->grp[i].mol[j] ].natom; /*number of atoms in this molecule*/
				else
					tmpi=model->grp[i].natom;
				tmpj=0; 
				if(model->grp[i].linear>0) model->mol[ model->grp[i].mol[j] ].linear=tmpj=1; 
				switch(tmpi) {/*check number of atoms in each molecule*/
					case 1: //monoatomic
						vacDF[i][vtrans]+=3.0;
						vacDF[i][vrotat]+=0.0;
						vacDF[i][vangul]+=0.0;
						vacDF[i][vimvib]+=0.0;
						vacDF[i][vtotal]+=3.0;
						break;
					case 2: //diatomic
						vacDF[i][vtrans]+=3.0;
						vacDF[i][vrotat]+=2.0;
						vacDF[i][vangul]+=2.0;
						vacDF[i][vimvib]+=1.0;
						vacDF[i][vtotal]+=6.0;
						break;
					default: //polyatomic
						vacDF[i][vtrans]+=3.0;
						vacDF[i][vrotat]+=(3.0-tmpj);
						vacDF[i][vangul]+=(3.0-tmpj);
						vacDF[i][vimvib]+=(3.0*tmpi-6.0+tmpj);
						vacDF[i][vtotal]+=(3.0*tmpi);
						break;
				}
				//cout<<"group "<<i+1<<" molecule "<<j<<" vacDF_trans "<<vacDF[i][vtrans]<<" vacDF_tot "<<vacDF[i][vtotal]<<endl;
			}
		}
	}

	//now correct the degrees of freedom for stretching constraints
	if(model->ctl.ana_vac_const[0]=='G'||model->ctl.ana_vac_const[0]=='g') ;
	else if(model->ctl.ana_vac_const_flag==1) model->grp[ngrp-1].constraint= atoi(model->ctl.ana_vac_const);
	//else model->grp[ngrp-1].constraint= 0;

	for(i=0;i<ngrp;i++) {
		if(! cmol) vacDF[i][vactype-1] -= model->grp[i].constraint;
		else {
			vacDF[i][vimvib] -= model->grp[i].constraint;
			vacDF[i][vtotal] -= model->grp[i].constraint;
		}
		//cout<<"group "<<i+1<<" final_vacDF "<<vacDF[i][vtotal]<<endl;
	}

}

/* ---------------------------------------------------------------------- */
void ComputeVAC::rd_frame(MODEL *model, int iframe) //mass weighted vac (vac_driver approach)
{ 
	int i,k,id,grpid;
	double atomEng,tot;
	double grpEng[ngrp];

	if(cmol==0) {
		for(i=0;i<vacatmpass;i++) {
			id=vacread[i];
			for(k=0;k<3;k++) { 
				vacvv[0][iframe][i][k]=(float)(model->atom[id].vel[k]);
			}
		}
	} else {
		model->cal_vcomp();
		for(i=0;i<vacatmpass;i++) {
			id=vacread[i];
			for(k=0;k<3;k++) {
				vacvv[vtrans][iframe][i][k]=(float)(model->atom[id].vt[k]); 
				vacvv[vrotat][iframe][i][k]=(float)(model->atom[id].vr[k]); 
				vacvv[vimvib][iframe][i][k]=(float)(model->atom[id].vv[k]); 
				vacvv[vtotal][iframe][i][k]=(float)(model->atom[id].vel[k]); 
			}
		} 
		id=vaclastmol-vacmolpass; //first molecular id
		for(i=id;i<vaclastmol;i++) //i is the molecular id
			for(k=0;k<3;k++) 
				vacvv[vangul][iframe][i-id][k]=(float)model->mol[i].anguv[k];
		if (vacipass==0)
			for(i=0;i<model->nmol;i++) 
				for(k=0;k<3;k++) pI[i][k]+=model->mol[i].pI[k];
	}
	
	if (vacipass==0) {
		trjT+=model->prp.T;
		vacE[ngrp-1][vactype-1] +=model->prp.Et;
		if(! model->ctl.ana_vac_eng_flag) model->grp[ngrp-1].eng.add_data(model->prp.Et);
		vacP+=model->prp.P;
		vacV[ngrp-1]+=model->prp.V;
		for(i=0;i<ngrp;i++) grpEng[i] = 0.0;
	}

	if (trj_atom_eng) {
		tot = 0.0;
		for(i=0;i<vacatmpass;i++) {
			id=vacread[i];
			grpid=model->atom[id].grp;
			atomEng = model->atom[id].eng_tmp;
			model->atom[id].eng.add_data(atomEng);
			grpEng[grpid] += atomEng;
			tot += atomEng;
		}
		for (i=0;i<ngrp-1;i++)
			model->grp[i].eng.add_data(grpEng[i]);
		model->grp[ngrp-1].eng.add_data(tot);
	}
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::compute(ostream *outf, MODEL *model) 
{ 
	int i,j,k,l,c,tmpi;
	char null[1024];

	tatom=model->atom;
	model->atom=vacatom;

	if(cmol) 
		for(i=0;i<model->nmol;i++) 
			for(j=0;j<model->mol[i].natom;j++) 
				model->mol[i].atom[j] = & vacatom[ model->mol[i].atm[j] ];

	if(! trj_atom_eng) {
		for(i=0;i<model->trj.ntrjf;i++) {
			if(model->trj.strj[i].has_atom_eng) {
				model->trj.strj[i].atom=model->atom;
				model->trj.strj[i].ngrp=ngrp;
				model->trj.strj[i].grp=model->grp;
				model->trj.strj[i].natom=natom;
				model->trj.strj[i].cell=&(model->cell);
				model->trj.strj[i].prp =&(model->prp);
				model->trj.strj[i].rd_atom_eng();
				has_atom_eng = 1;
			}
		}
	}
	if(cmol)
		for(i=0;i<model->nmol;i++) 
			for(k=0;k<3;k++) pI[i][k]/=nframe; //average principle moment of inertia

	vacipass++;
	calc_vacv(model);

	for (c=1;c<vacnpass;c++) {
		vacatmpass =0;
		vacmolpass =0;
		//determine which atoms to read
		if(cmol==0) {
			j= natom - vactatmpassed;
			if(j>vacatmperpass) j=vacatmperpass;
			for(i=0;i<j;i++) vacread[vacatmpass++]=vactatmpassed+i;
		} else {
			tmpi=0;
			for(i=vaclastmol;i<model->nmol;i++) {
				tmpi += model->mol[i].natom;
				if(tmpi<=vacatmperpass) {
					for(j=0;j<model->mol[i].natom;j++) vacread[vacatmpass++]=model->mol[i].atm[j];
					vacmolpass++;
				} else break; //break i-loop
			}
			vactmolpassed+=vacmolpass;
			vaclastmol=i;
		}

		//read the trj file again
		for(i=0;i<vactype;i++) 
			for(j=0;j<vacnsteps;j++) 
				for(k=0;k<vacatmperpass;k++) 
					for(l=0;l<3;l++) 
						vacvv[i][j][k][l]=0;
		l=0;
		for(j=model->ctl.ana_iframe;j<=model->ctl.ana_fframe;j+=model->ctl.ana_sframe) {
			cout<<" Reading frame "<<j<<"/"<<model->ctl.ana_fframe<<"\xd";
			fflush(stdout);
			model->atom=tatom;
			if(cmol) 
				for(i=0;i<model->nmol;i++) 
					for(k=0;k<model->mol[i].natom;k++) 
						model->mol[i].atom[k] = & tatom[ model->mol[i].atm[k] ];
			for(i=0;i<model->trj.ntrjf;i++) {
				model->trj.strj[i].atom=model->atom;
				model->trj.strj[i].cell=&vaccell;
				model->trj.strj[i].prp =&vacprp;
			}
			model->trj.rd_frame(j);
			rd_frame(model,l);
			l++;
		}
		cout<<" Reading frame "<<model->ctl.ana_fframe<<"/"<<model->ctl.ana_fframe<<endl;
		vacipass++;
		vactatmpassed +=vacatmpass;
		calc_vacv(model);
	}
	memory->destroy_4d_float_array(vacvv);

	cout<<" Calculating VAC Thermodynamics"<<endl;
	calc_thermo_vals(outf,model);
	compute_vel_temp();
	//calculate constant volume heat capacity
	for(i=0;i<ngrp;i++) {
		Cvc[i][vactype-1] = model->grp[i].eng.sigma2x*1E6*caltoj*caltoj/(R*vacT[i][vactype-1]*vacT[i][vactype-1]);
		for (j=0;j<vactype-1;j++) Cvc[i][j] = Cvc[i][vactype-1] * vacDF[ngrp-1][j]/vacDF[ngrp-1][vactype-1];
	}
	calc_pwr_spectrum();

	if(vacV[ngrp-1]<=0) {
		sprintf(null,"**WARNING: original volume %lf ",vacV[ngrp-1]);
		cout<<null<<endl;
		vacV[ngrp-1]= 9.9668*natom;
		sprintf(null,"Estimate volume to be %lf from water 9.9668 A3/atom **",vacV[ngrp-1]);
		cout<<null<<endl;
	}

	//apply 2PT correction
	do_2pt(outf,model);

	//vibrational analysis
	do_vib_analysis();

}

/* ---------------------------------------------------------------------- */
void ComputeVAC::calc_thermo_vals(ostream *outf, MODEL *model) 
{

	int i,j,k,ngmol;
	double ediff,scale_f,grp_tot,eng_tot;
	char null[1024];

	j = 0;
	//volume
	if(!model->ctl.ana_vac_vol_flag) 
		vacV[ngrp-1]/=nframe;
	else 
		vacV[ngrp-1] = model->ctl.ana_vac_vol;
	vacV[ngrp-1]-=model->ctl.ana_vac_void_vol; //adjust the volume based on fictitious void volume
	for(i=0;i<ngrp-1;i++) {
		ngmol= model->grp[i].nmol; 
		vacV[i] = ngmol*vacV[ngrp-1]/nmol; 
		if(model->grp[i].vol>0.0) { 
			vacV[i] = model->grp[i].vol; 
			if (DEBUG) cout<<"grp "<<i<<" volume "<<vacV[i]<<endl;
		}
	}

	//pressure
	if(model->ctl.ana_vac_press_flag) vacP = model->ctl.ana_vac_press;
	else vacP/=nframe; //average pressure (GPa)

	//temperature
	if(model->ctl.ana_vac_temp_flag) trjT = model->ctl.ana_vac_temp;
	else trjT/=nframe; //average temperature (K)

	//trdf will be evenly distributed to each degrees of freedom
	for(i=0;i<ngrp;i++) 
		for(j=0;j<vactype;j++) 
			vacDF[i][j]-= trdf*vacDF[i][j]/vacDF[ngrp-1][vactype-1];

	//energy
	k = vactype - 1;
	eng_tot = vacE[ngrp-1][vactype-1]/nframe;
	if(ngrp==1||(!has_atom_eng)) {
		if(! model->ctl.ana_vac_eng_flag) model->grp[ngrp-1].eng.cal_avg();
		vacE[ngrp-1][k] = model->grp[ngrp-1].eng.avg;
		for(i=0;i<ngrp-1;i++) { //distribute total energy to each grp based on dof (ideal gas/perfect solid approximation)
			model->grp[i].eng.init();
			model->grp[i].eng.dataflag = 1;
			model->grp[i].eng.nXY = model->grp[i].eng.npt = 1;
			vacE[i][vactype-1] = model->grp[i].eng.sum = vacE[ngrp-1][k] * vacDF[i][k]/vacDF[ngrp-1][k];
			model->grp[i].eng.sigma2x = model->grp[ngrp-1].eng.sigma2x * vacDF[i][k]/vacDF[ngrp-1][k];
		}
	} else {
		for(i=0;i<ngrp;i++) {
			if (! model->ctl.ana_vac_eng_flag) model->grp[i].eng.cal_avg();
			vacE[i][vactype-1] = model->grp[i].eng.avg;
		}
	}
	//get energy of each group and distribute by dof
	grp_tot = 0.0;
	for(i=0;i<ngrp;i++) {
		if(i<(ngrp-1)) grp_tot += vacE[i][vactype-1];
		for(k=0;k<vactype-1;k++) 
			vacE[i][k] = vacE[i][vactype-1] * vacDF[i][k]/vacDF[i][vactype-1];
	}
	ediff = (eng_tot-grp_tot)/eng_tot;
	scale_f = 1.0;
	if(fabs(ediff)>0.001&&ngrp>1) {
		scale_f = vacE[ngrp-1][vactype-1]/grp_tot;
		cout<<"WARNING: Total energy from trajectory: "<<eng_tot<<" and from groups ";
		cout<<grp_tot<<" differ by "<<(100*ediff)<<"%"<<endl<<" Will rescale group energies by "<<scale_f<<endl;
		for(i=0;i<ngrp-1;i++)
			for(j=0;j<vactype;j++)
				vacE[i][j] *= scale_f;
	}	

	sprintf(null,"%-20s :","Avg Thermo Vals:");
	*outf<<null<<endl; cout<<null<<endl;
	sprintf(null,"	%-12s : %.3f %s\n	%-12s : %.3f +/- %.3f %s\n	%-12s : %.3f %s\n	%-12s : %.3f %s",
		"Temperature",trjT,"(K)",
		"Energy",vacE[ngrp-1][vactype-1]*4.184,model->grp[ngrp-1].eng.std*scale_f*4.184,"(kJ/mol)",
		"Pressure",vacP,"(GPa)","Volume",vacV[ngrp-1],"(A^3)");
	*outf<<null<<endl;cout<<null<<endl;
	if(ngrp > 1) {
		//Energy
		sprintf(null,"%-12s :","Group Eg");
		*outf<<null<<endl; cout<<null<<endl;
		for(i=0;i<ngrp-1;i++) {
			sprintf(null,"	Group %d: %.3f +/- %.3f kJ/mol",(i+1),vacE[i][vactype-1]*caltoj,model->grp[i].eng.std*scale_f*caltoj); *outf<<null<<endl; cout<<null<<endl;
		}
		//Volume
		sprintf(null,"%-12s :","Group Vol");
		*outf<<null<<endl; cout<<null<<endl;
		for(i=0;i<ngrp-1;i++) {
			sprintf(null,"	Group %d: %.3f A^3",(i+1),vacV[i]); *outf<<null<<endl; cout<<null<<endl;
		}
	}
}
 
/* ---------------------------------------------------------------------- */
void ComputeVAC::calc_vacv(MODEL *model) 
{
	int i,id,tp;

	char null[1024];

	sprintf(null," Calculating VAC... Pass %d of %d: atoms %d/%d ",vacipass,vacnpass,vactatmpassed,natom);
	cout<<null;
	if(cmol) { 
		sprintf(null,"mols %d/%d ",vactmolpassed,model->nmol);
		cout<<null;
	}
	fflush(stdout);
	//determine VAC
	for(tp=0;tp<vactype;tp++) {
		if(tp!=vangul) //translation, rotation, vibration, total
			//#pragma omp parallel for 
			for(i=0;i<vacatmpass;i++)
				calc_acf(i,tp,model);
		else { //angular velocity
			id=vaclastmol-vacmolpass; //first molecular id
			//#pragma omp parallel for
			for(i=id;i<vaclastmol;i++) //i is the molecular id
				calc_acf_ang(i,id,model);
		}
	}
	cout<<"Done"<<endl;
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::calc_acf(int i, int tp, MODEL *model) //compute the atomic autocorreclation func
{
	int j,k,g1,g2;

	int id = vacread[i];
	double mass = model->atom[id].mass;
	STAT *atmvac= new STAT [tot_N];

	for(k=0;k<3;k++) {
		for(j=0;j<vacnsteps;j++) {
			in[j][0]= vacvv[tp][j][i][k];	//velocity in k direction of atom i in the j_th step
			in[j][1]= 0;
		}
		for(;j<tot_N;j++) in[j][0]=in[j][1]=0;
		fftw_execute(p1); // do backward(+1) fft

		for(j=0;j<tot_N;j++) {
			in[j][0]= (out[j][0]*out[j][0] + out[j][1]*out[j][1])/tot_N; //work(j)=work(j)*dconjg(work(j))/tot_N
			in[j][1]= 0;
		}
		fftw_execute(p2); // do forward(-1) fft

		for(j=0;j<tot_N;j++) 
			atmvac[j].add_data(mass*out[j][0]/vacnsteps);
	}

	for(g1=0;g1<ngrp;g1++) 
		for(g2=0;g2<model->grp[g1].natom;g2++) 
			if(id==model->grp[g1].atm[g2]) 
				for(j=0;j<tot_N;j++)
					Pvac[g1][tp][j].add_data(atmvac[j].sum);

	delete [] atmvac;
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::calc_acf_ang(int i, int id, MODEL *model) //compute the atomic autocorreclation func
{
	int j,k,g1,g2,tp;

	tp = vangul;
	STAT *atmvac= new STAT [tot_N];

	for(k=0;k<3;k++) {
		for(j=0;j<vacnsteps;j++) {
			in[j][0]= vacvv[tp][j][i-id][k];	//velocity in k direction of molecule i in the j_th step
			in[j][1]= 0;
		}
		for(;j<tot_N;j++) in[j][0]=in[j][1]=0;
		fftw_execute(p1); // do backward(+1) fft

		for(j=0;j<tot_N;j++) {
			in[j][0]= (out[j][0]*out[j][0] + out[j][1]*out[j][1])/tot_N; //work(j)=work(j)*dconjg(work(j))/tot_N
			in[j][1]= 0;
		}
		fftw_execute(p2); // do forward(-1) fft

		for(j=0;j<tot_N;j++) 
			atmvac[j].add_data(1.0*out[j][0]/vacnsteps); //use unity for mass, angular vel has been weighted by sqrt(I).
	}

	for(g1=0;g1<ngrp;g1++) 
		for(g2=0;g2<model->grp[g1].nmol;g2++) 
			if(i==model->grp[g1].mol[g2]) 
				for(j=0;j<tot_N;j++)	
					Pvac[g1][tp][j].add_data(atmvac[j].sum);

	delete [] atmvac;
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::compute_vel_temp() //mass weighted vac (vac_driver approach)
{ 
		int i,j,k;

		cout<<"Calculating Temperature from velocities...";
		fflush(stdout);
		for(i=0;i<vacmaxf;i++) 
			 for(j=0;j<ngrp;j++) 
					for(k=0;k<vactype;k++) 
						Pvac[j][k][i].cal_avg();
						//Pvac[j][k][i].sum is the sum of vac of all atoms 
						//Pvac[j][k][i].npt is the number of atoms

		//determine the real temperature of each group
		//real temp from df*kT/2=sum(mv^2/2), sum(mv^2)=Pvac[j][0]
		for(j=0;j<ngrp;j++) 
				for(i=0;i<vactype;i++) {
					if((int)(vacDF[j][i])==0)	vacT[j][i] = 0;
					else vacT[j][i] = Pvac[j][i][0].sum/(0.1*R*vacDF[j][i]);
				}
		
		i=ngrp-1; j=vactype-1;
		cout<<"Done"<<endl;
		if(fabs(vacT[i][j]-trjT)/trjT > 0.0001) {
			cout<<" Trajectory temperature "<<trjT<<" and velocity temperature "<<vacT[i][j]<<" differ by "<<(vacT[i][j]-trjT)/trjT*100.0<<"%"<<endl;
			cout<<" Check the degrees of freedom calculation in the vac code"<<endl;
		}
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::calc_pwr_spectrum()
{ 
	//determine the power spectrum
	int i,j,k;
 
	cout<<"Calculating Power Spectrum using FFTs...";
	fflush(stdout);

	pwrfreq=1.0E10/(vacdtime*tot_N*vlight);	//in cm-1
	for(j=0;j<ngrp;j++) 
		for(k=0;k<vactype;k++) {
			if(vacT[j][k]<=0) for(i=0;i<tot_N;i++) in[i][0]=in[i][1]=0;
			else 
				for(i=0;i<tot_N;i++)	{
					in[i][0]=Pvac[j][k][i].sum*20.0/(R*vacT[j][k]); //work(i)=dcmplx(vac(i)*20.0/(gascon*ave_temp),0.d0)
					in[i][1]=0; // in[] is dimensionless
				 } 
			fftw_execute(p1);

			for(i=0;i<nused;i++) 
				Pvac[j][k][i].max=out[i][0]*vacdtime*vlight*1.0E-10; //pwr(i)=dble(work(i))*time_step*cspeed

		Pvac[j][k][0].min = Pvac[j][k][0].max*pwrfreq*0.5;	//use min to store PWR integration
		for(i=1;i<nused-1;i++) Pvac[j][k][i].min = Pvac[j][k][i-1].min + Pvac[j][k][i].max*pwrfreq;
		Pvac[j][k][i].min = Pvac[j][k][i-1].min + Pvac[j][k][i].max*pwrfreq*0.5;
	 }
		
	cout<<"Done"<<endl;
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::do_2pt(ostream *outf,MODEL *model)
{ 

	int i,j;
	char null[1024];

	int ntmol,ngmol;
	double gmass,rT[3],rs,alpha;
	double ry[2],hs_sigma[2];
	double tmpg,tmps,gvol;

	double ***pwr;
	double **y,**ttdf;
	double wsehs[ngrp];

	pwr = memory->create_3d_double_array(2,ngrp,nused,"vacpwr");
	y = memory->create_2d_double_array(2,ngrp,"vacy");
	ttdf = memory->create_2d_double_array(2,ngrp,"vacttdf");

	cout<<"Applying 2PT correction...";
	fflush(stdout);

	ntmol= model->nmol; //number of molecules in the whole system
	for(j=0;j<ngrp;j++) {
		ngmol= model->grp[j].nmol; //number of molecules in the group of interest
		if(! cmol) ngmol = 1;
		gmass= model->grp[j].mass; //total mass of atoms (molecules) in the group of interest
		rs=model->grp[j].rotsym*1.0; //rotational symmetry number
		gvol = vacV[j];
		rT[0]=rT[1]=rT[2]=0.0;
		if(cmol) {
			get_rot_temp(model,j,rT);
			for(i=0;i<nused;i++) pwr[1][j][i]=Pvac[j][vangul][i].max; //rotation in cm
		}
		//cout<<"grp "<<j+1<<" ngmol "<<ngmol<<" gmass "<<gmass<<" rs "<<rs<<" gvol "<<gvol<<endl;
		for(i=0;i<nused;i++) pwr[0][j][i]=Pvac[j][vtrans][i].max; //translation in cm
		//now obtain translational diffusivity K2pt - eqn 33
		K2pt[0][j]=(pwr[0][j][0]/vlight*1E-2)/ngmol*sqrt(PI*Na*kb*vacT[j][vtrans]/(gmass*1e-3/ngmol))*2.0/9.0*pow(ngmol/gvol,1.0/3.0)*1E10*pow(6.0/PI,2.0/3.0); 
		//now obtain fluidicity self consitently by solving cubic equation: eqn 34 
		f2pt[0][j]=search2PT(K2pt[0][j]);
		ry[0]=pow(f2pt[0][j]/K2pt[0][j],1.5); //packing fraction y, eqn 32
		hs_sigma[0]=pow(ry[0]*6.0/PI/ngmol*gvol,1.0/3.0); //hard sphere diameter sigma: y = (pi/6)*rho*hs_sigma^3
		hsdf[0][j]=HSDF(pwr[0][j],nused,pwrfreq,ngmol,f2pt[0][j]); /*determine hard sphere degrees of freedom hsdf*/
		ttdf[0][j]=TTDF(pwr[0][j],nused,pwrfreq); /*determine total degrees of freedom ttdf*/
		y[0][j]=ry[0]*(hsdf[0][j]/ttdf[0][j]);
		alpha=12*hsdf[0][j]/3.0/Pvac[j][vtrans][i].max; //Eskog friction constant, from eqn 23
		if(rT[0]>0.0) {
			K2pt[1][j]=(pwr[1][j][0]/vlight*1E-2)/ngmol*sqrt(PI*Na*kb*vacT[j][vangul]/(gmass*1e-3/ngmol))*2.0/9.0*pow(ngmol/gvol,1.0/3.0)*1E10*pow(6.0/PI,2.0/3.0); //rotational diffusivity
			f2pt[1][j]=search2PT(K2pt[1][j]);
			ry[1]=pow(f2pt[1][j]/K2pt[1][j],1.5);
			hs_sigma[1]=pow(ry[1]*6.0/PI/ngmol*gvol,1.0/3.0);
			hsdf[1][j]=HSDF(pwr[1][j],nused,pwrfreq,ngmol,f2pt[1][j]); /*determine hard sphere degrees of freedom hsdf*/
			ttdf[1][j]=TTDF(pwr[1][j],nused,pwrfreq); /*determine total degrees of freedom ttdf*/
			y[1][j]=ry[1]*(hsdf[1][j]/ttdf[1][j]);
		} else K2pt[1][j]=f2pt[1][j]=ry[1]=hs_sigma[1]=hsdf[1][j]=ttdf[1][j]=y[1][j]=0;

		if(y[0][j]>0.74) f2pt[0][j]=f2pt[1][j]=0; /*packing fraction too large, ignore partition*/
		HSweighting(&wep[j],&wap[j],&wcvp[j],&wsehs[j],&wsp[j],&wspd[j],y[0][j],gmass/ngmol,hsdf[0][j]/3.0,vacT[j][vtrans],vacT[j][vangul],gvol,&wsr[j],&war[j],rT,rs);

		if(DEBUG) {
			sprintf(null,"\n rotational temperatures %e %e %e K\n",rT[0],rT[1],rT[2]);
			*outf<<null<<endl;
			sprintf(null,"trans s0=%lf T=%lf mass=%lf V=%lf ngmol=%d ntmol=%d K=%lf",pwr[0][j][0],vacT[j][vtrans],gmass/ngmol,gvol,ngmol,ntmol,K2pt[0][j]);
			*outf<<null<<endl;
			sprintf(null,"rotat s0=%lf T=%lf mass=%lf V=%lf ngmol=%d ntmol=%d K=%lf",pwr[1][j][0],vacT[j][vangul],gmass/ngmol,gvol,ngmol,ntmol,K2pt[1][j]);
			*outf<<null<<endl;
			sprintf(null,"trans fludi packfrac hsdf tdf	%lf %lf %lf %lf",f2pt[0][j],hsdf[0][j]/ttdf[0][j],hsdf[0][j],ttdf[0][j]);
			*outf<<null<<endl;
			sprintf(null,"rotat fludi packfrac hsdf tdf	%lf %lf %lf %lf",f2pt[1][j],hsdf[1][j]/ttdf[1][j],hsdf[1][j],ttdf[1][j]);
			*outf<<null<<endl;
			sprintf(null,"constant K  %lf %lf",K2pt[0][j],K2pt[1][j]);
			*outf<<null<<endl;
			sprintf(null,"fludicity %lf %lf",f2pt[0][j],f2pt[1][j]);
			*outf<<null<<endl;
			sprintf(null,"refernce Hard Sphere packing  %lf %lf",ry[0],ry[1]);
			*outf<<null<<endl;
			sprintf(null,"Hard Sphere diameter  %lf %lf",hs_sigma[0],hs_sigma[1]);
			*outf<<null<<endl;
			sprintf(null,"Hard Sphere packing fraction %lf %lf",y[0][j],y[1][j]);
			*outf<<null<<endl;
			sprintf(null,"Effective N for hard sphere  %lf %lf",ngmol*(hsdf[0][j]/(3.0*ntmol-3.0)),ngmol*(hsdf[1][j]/(3.0*ntmol-3.0)));
			*outf<<null<<endl;
			sprintf(null,"Hard Sphere reduced density %lf %lf",y[0][j]*6.0/PI,y[1][j]*6.0/PI);
			*outf<<null<<endl;
			sprintf(null,"True Hard Sphere dof %lf %lf",hsdf[0][j],hsdf[1][j]);
			*outf<<null<<endl;
			sprintf(null,"temperature rotT3    %lf %lf %lf %lf",vacT[j][vangul],rT[0],rT[1],rT[2]);
			*outf<<null<<endl;
			sprintf(null,"Weighting Stran Srot %lf %lf\n",wsp[j],wsr[j]);
			*outf<<null<<endl;
		}
		for(i=0;i<nused;i++) {
			twoPT(&tmpg,&tmps,pwr[0][j][0],pwr[0][j][i],pwrfreq*i,ngmol,f2pt[0][j]); //get frequency dependent hs DoS - eqn 24
			Pvac[j][vtrans][i].a=tmpg; //use a to store hard sphere gas contribution
			Pvac[j][vtrans][i].b=tmps; //use b to store solid contribution
			if(rT[0]>0) twoPT(&tmpg,&tmps,pwr[1][j][0],pwr[1][j][i],pwrfreq*i,ngmol,f2pt[1][j]); 
			else tmpg=tmps=0;
			if (cmol) { 
				Pvac[j][vangul][i].a=tmpg; //use a to store hard sphere gas contribution
				Pvac[j][vangul][i].b=tmps; //use b to store solid contribution
			}
		}
		Pvac[j][vtrans][0].SEa = Pvac[j][vtrans][0].a*pwrfreq*0.5;	//use SEa to store PWR integration of hs 
		Pvac[j][vtrans][0].SEb = Pvac[j][vtrans][0].b*pwrfreq*0.5;	//use SEb to store PWR integration of solid 
		if(cmol) {
			Pvac[j][vangul][0].SEa = Pvac[j][vangul][0].a*pwrfreq*0.5;	//use SEa to store PWR integration of hs 
			Pvac[j][vangul][0].SEb = Pvac[j][vangul][0].b*pwrfreq*0.5;	//use SEb to store PWR integration of solid 
		}
		for(i=1;i<nused-1;i++) { 
			Pvac[j][vtrans][i].SEa = Pvac[j][vtrans][i-1].SEa + Pvac[j][vtrans][i].a*pwrfreq;
			Pvac[j][vtrans][i].SEb = Pvac[j][vtrans][i-1].SEb + Pvac[j][vtrans][i].b*pwrfreq;
			if(cmol) {
				Pvac[j][vangul][i].SEa = Pvac[j][vangul][i-1].SEa + Pvac[j][vangul][i].a*pwrfreq;
				Pvac[j][vangul][i].SEb = Pvac[j][vangul][i-1].SEb + Pvac[j][vangul][i].b*pwrfreq;
			}
		}
		Pvac[j][vtrans][i].SEa = Pvac[j][vtrans][i-1].SEa + Pvac[j][vtrans][i].a*pwrfreq*0.5;
		Pvac[j][vtrans][i].SEb = Pvac[j][vtrans][i-1].SEb + Pvac[j][vtrans][i].b*pwrfreq*0.5;
		if(cmol) {
			Pvac[j][vangul][i].SEa = Pvac[j][vangul][i-1].SEa + Pvac[j][vangul][i].a*pwrfreq*0.5;
			Pvac[j][vangul][i].SEb = Pvac[j][vangul][i-1].SEb + Pvac[j][vangul][i].b*pwrfreq*0.5;
		}
	}
	cout<<"Done"<<endl;

	memory->destroy_3d_double_array(pwr);
	memory->destroy_2d_double_array(y);
	memory->destroy_2d_double_array(ttdf);
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::do_vib_analysis()
{ 
	int i,j,k,l;
	double tmpg,tmps;
	double scaled_temp; // kb*vacT/(h*vlight*100); //scaled_temp= k*T/h/cspeed, unit 1/cm 
	double u; //scaled frequency hcv/kT= hc pwrfreq*i /kT
	double weq,wec,wsq,waq,wcvq,wcvc;
	double tmpTT = 0.0;

	weq=wsq=waq=wcvq=0;
	wcvc=wec=0;
	//property for zero frequency
	for(k=0;k<vactype;k++) { 
		i=0;
		for(j=0;j<ngrp;j++) { 
			tmpTT=vacT[j][k];
			scaled_temp = kb*tmpTT/(h*vlight*100);
			if(vacT[j][k]>0) u=pwrfreq*0.5/scaled_temp;
			else u=0;
			if(k==vtrans || k == vangul) {
				tmpg=Pvac[j][k][i].a; //gas contribution
				tmps=Pvac[j][k][i].b; //solid contribution
			} else { 
				tmpg=0; //gas-like fraction, from 2pt
				tmps=Pvac[j][k][i].max; //solid-like fraction
				if(tmps!=tmps) tmps=Pvac[j][k][i].max=u=0; //??
			}
			thermo[0][j][k] =0.5*(tmps+tmpg);	 //degrees of freedom
			thermo[1][j][k] =0.5*(tmps*0.0);	//ZPE
			thermo[2][j][k] =0.5*(tmps*1.0+tmpg*wep[j]);	 //Eq
			thermo[3][j][k] =0.5*(tmps*1.0+tmpg*wep[j]);	 //Ec
			thermo[4][j][k] =0.5*(tmps*sqweighting(u)+tmpg*wsp[j]);	 //Sq
			thermo[6][j][k] =0.5*(tmps*(1-sqweighting(u))+tmpg*wap[j]); //Aq
			thermo[8][j][k] =0.5*(tmps*1.0+tmpg*wcvp[j]);	//Cvq
			thermo[9][j][k] =0.5*(tmps*1.0+tmpg*wcvp[j]);	//Cvc
			thermo[13][j][k]=0.5*(tmps*sqweighting(u)+tmpg*wspd[j]);	//Sd
			if(k==vangul) {
				thermo[2][j][k] += 0.5*tmpg*(wer[j]-wep[j]);	//Eq
				thermo[3][j][k] += 0.5*tmpg*(wer[j]-wep[j]);	//Ec
				thermo[4][j][k] += 0.5*tmpg*(wsr[j]-wsp[j]);	//Sq
				thermo[6][j][k] += 0.5*tmpg*(war[j]-wap[j]);	//Aq
				thermo[8][j][k] += 0.5*tmpg*(wcvr[j]-wcvp[j]);	//Cvq
				thermo[9][j][k] += 0.5*tmpg*(wcvr[j]-wcvp[j]);	//Cvc
			}
		}
			 
		//all other frequencies
		for(i=1; i<nused; i++) {
			for(j=0;j<ngrp;j++) {
				tmpTT=vacT[j][k];
				scaled_temp = kb*tmpTT/(h*vlight*100);
				if(scaled_temp>0.0) {
					u=pwrfreq*i/scaled_temp;
					weq=u/2.0+u/(exp(u)-1.0);
					wec=1.0;
					wsq=u/(exp(u)-1)-log(1-exp(-u));
					waq=log((1-exp(-u))/exp(-u/2));
					wcvq=(u*u*exp(u))/(1-exp(u))/(1-exp(u));
					wcvc=1.0;
				} else u=weq=wec=wsq=waq=wcvq=wcvc=0;

				if(k==vtrans) {
					tmpg=Pvac[j][k][i].a; //gas contribution
					tmps=Pvac[j][k][i].b; //solid contribution
				} else if(k==vangul) { 
					tmpg=Pvac[j][k][i].a; //gas contribution
					tmps=Pvac[j][k][i].b; //solid contribution
				} else { 
					tmpg=0; //gas-like fraction, from 2pt
					tmps=Pvac[j][k][i].max; //solid-like fraction
					if(tmps!=tmps) tmps=u=weq=wsq=waq=wcvc=wcvq=Pvac[j][k][i].max=0;
				}
				if(i==nused-1) { tmpg*=0.5; tmps*=0.5; }
				thermo[0][j][k] +=(tmps+tmpg);	//degrees of freedom
				thermo[1][j][k] +=(tmps*u*0.5);	//ZPE
				thermo[2][j][k] +=(tmps*weq+tmpg*wep[j]);	//Eq
				thermo[3][j][k] +=(tmps*wec+tmpg*wep[j]);	//Ec
				thermo[4][j][k] +=(tmps*wsq+tmpg*wsp[j]);	//Sq
				thermo[6][j][k] +=(tmps*waq+tmpg*wap[j]);	//Aq
				thermo[8][j][k] +=(tmps*wcvq+tmpg*wcvp[j]);	//Cvq
				thermo[9][j][k] +=(tmps*wcvc+tmpg*wcvp[j]);	//Cvc
				thermo[13][j][k]+=(tmps*wsq+tmpg*wspd[j]);	//Sdq
				if(k==vangul) {
					thermo[2][j][k] += tmpg*(wer[j]-wep[j]);	//Eq
					thermo[3][j][k] += tmpg*(wer[j]-wep[j]);	//Ec
					thermo[4][j][k] += tmpg*(wsr[j]-wsp[j]);	//Sq
					thermo[6][j][k] += tmpg*(war[j]-wap[j]);	//Aq
					thermo[8][j][k] += tmpg*(wcvr[j]-wcvp[j]);	//Cvq
					thermo[9][j][k] += tmpg*(wcvr[j]-wcvp[j]);	//Cvc
				}
			}
		}
		
		for(j=0;j<ngrp;j++) {
			thermo[0][j][k]*= pwrfreq;	//total degrees of freedom
			thermo[1][j][k]*= pwrfreq*tmpTT*R*1.0e-3;	//zero point energy kJ/mol/SimBox
			thermo[2][j][k]*= pwrfreq*tmpTT*R*1.0e-3;	//quantum energy  kJ/mol/SimBox
			thermo[3][j][k]*= pwrfreq*tmpTT*R*1.0e-3;	//classical energy  kJ/mol/SimBox
			thermo[4][j][k]*= pwrfreq*R;	//quantum entropy  J/(K mol)/SimBox
			thermo[6][j][k]*= pwrfreq*tmpTT*R*1.0e-3;	//quantum Helmholtz free energy kJ/mol/SimBox
			thermo[8][j][k]*= pwrfreq*R;	//quantum constant volume heat capacity J/(K mol)/SimBox
			thermo[9][j][k]*= pwrfreq*R;	//classical constant volume heat capacity J/(K mol)/SimBox
			thermo[10][j][k]= vacE[j][k]*caltoj;	//MD total (strain) energy kJ/mol/SimBox
			thermo[11][j][k]= thermo[10][j][k]-thermo[3][j][k];	//reference energy Eo kJ/mol/SimBox
			thermo[8][j][k]= Cvc[j][k]-(thermo[9][j][k]-thermo[8][j][k]);	//Cvq = Cvc - Cv'
			thermo[9][j][k]= Cvc[j][k];
			thermo[2][j][k]+= thermo[10][j][k]-thermo[3][j][k];	//Eq= Eq'+Eo
			thermo[6][j][k]+= thermo[10][j][k]-thermo[3][j][k];	//Aq= Aq'+Eo
			thermo[12][j][k]= Pvac[j][k][0].max;	//DoS at zero frequency cm
			thermo[13][j][k]*= pwrfreq*R;	//quantum entropy distinguishable particles J/(K mol)/SimBox
		}
	}
		
	if(cmol)
		//apply 2pt correction to system properties
 
		for(j=0;j<ngrp;j++) {
			i=vrotat; //store the 2pt results in vrotat - alternative way of getting rotations
			for(l=0;l<14;l++) thermo[l][j][i]=0;
			//ZPE,Eq,Ec,Sq,Sc,Aq,Ac,Cvq,Cvc,Emd,Vo,So
			for(k=0;k<i;k++) 
				for(l=0;l<14;l++)
					thermo[l][j][i]+=thermo[l][j][k];
			vacT[j][i]=vacT[j][vtotal];
		}
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::report(ostream *outf, MODEL *model) //mass weighted vac (vac_driver approach)
{ 

	int i,j,k;
	double x;
	int show2pt = model->ctl.ana_vac_show_2pt;
	char null[1024],filename[1024];

	strcpy(filename,model->ctl.ana_out);
	strcat(filename,".vac");
	ofstream outvac(filename,ios::out);

	strcpy(filename,model->ctl.ana_out);
	strcat(filename,".pwr");
	ofstream outpwr(filename,ios::out);

	strcpy(filename,model->ctl.ana_out);
	strcat(filename,".3n");
	ofstream out3n(filename,ios::out);

	strcpy(filename,model->ctl.ana_out);
	strcat(filename,".thermo");
	ofstream outthermo(filename,ios::out);

	//print some general info
	copyright(outf);
	copyright(&outthermo);

	compute(outf,model);

	sprintf(null," Velocity autocorrelation function (g/mol*A^2/ps^2) analysis");
	*outf<<null<<endl; outvac<<null<<endl;
	sprintf(null," %13s","time(ps)"); 
	*outf<<null; outvac<<null;
	for(j=0;j<ngrp;j++) {
		if(! cmol) { 
			sprintf(null,"	%6s[G%03d]","VAC",j+1); 
			*outf<<null; outvac<<null;
		} else {
			sprintf(null,"	%6s[G%03d]	%6s[G%03d]	%6s[G%03d]	%6s[G%03d]","VACcmt",j+1,"VACcmo",j+1,"VACimv",j+1,"VACcmr",j+1); 
			*outf<<null; outvac<<null;
			sprintf(null,"	%6s[G%03d]","VACtot",j+1); 
			*outf<<null; outvac<<null;
		}
	}
	*outf<<endl; outvac<<endl;

	for(i=0;i<vacmaxf;i++) {
		sprintf(null," %13.3f",vacdtime*i); 
		*outf<<null; outvac<<null;
		for(j=0;j<ngrp;j++) {
			for(k=0;k<vactype;k++) {
				sprintf(null," %13.3f",Pvac[j][k][i].sum); 
				*outf<<null; outvac<<null;
				//Pvac[j][k][i].sum is the sum of vac of all atoms 
				//Pvac[j][k][i].npt is the number of atoms
			}
		}
		*outf<<endl; outvac<<endl;
	}

	sprintf(null,"\n Power Spectrum (Density of State) (cm)");
	*outf<<null<<endl; outpwr<<null<<endl;
	sprintf(null," %13s","freq(cm-1)"); 
	*outf<<null; outpwr<<null;
	for(j=0;j<ngrp;j++) {	
		sprintf(null,"  %6s[G%03d]  %6s[G%03d]","PWR_hs",j+1,"PWR_st",j+1); 
		*outf<<null; outpwr<<null; 
		if(cmol) { 
			sprintf(null,"  %6s[G%03d]  %6s[G%03d]","PWR_fr",j+1,"PWR_sr",j+1); *outf<<null; outpwr<<null; 
			sprintf(null,"  %6s[G%03d]  %6s[G%03d]  %6s[G%03d]  %6s[G%03d]","PWRcmt",j+1,"PWRcmo",j+1,"PWRimv",j+1,"PWRcmr",j+1); 
			*outf<<null; outpwr<<null; 
		}
		sprintf(null,"  %6s[G%03d]","PWRtot",j+1); *outf<<null; outpwr<<null;
	}
	*outf<<endl; outpwr<<endl;
	for(i=0;i<nused;i++) { 
		sprintf(null," %13.4f",pwrfreq*i); 
		*outf<<null; outpwr<<null;
		for(j=0;j<ngrp;j++) {
			k=0;
			sprintf(null," %13.4f %13.4f",Pvac[j][k][i].a,Pvac[j][k][i].b); 
			*outf<<null; outpwr<<null; 
			k=vangul;
			if(cmol) { 
				sprintf(null," %13.4f %13.4f",Pvac[j][k][i].a,Pvac[j][k][i].b); 
				*outf<<null; outpwr<<null; 
			}
			for(k=0;k<vactype;k++) { 
				sprintf(null," %13.4f",Pvac[j][k][i].max); 
				*outf<<null; outpwr<<null;
			}
		}
		*outf<<endl; outpwr<<endl;
	}

	cout<<"Saving results...";
	fflush(stdout);
	sprintf(null,"\n Integration of Power Spectrum (cm)");
	*outf<<null<<endl; out3n<<null<<endl;

	sprintf(null," %13s","freq(cm-1)"); 
	*outf<<null; out3n<<null;
	for(j=0;j<ngrp;j++) {
		sprintf(null,"  %6s[G%03d]  %6s[G%03d]","INT_hs",j+1,"INT_st",j+1); 
		*outf<<null; out3n<<null;
		if(cmol) {
			sprintf(null,"  %6s[G%03d]  %6s[G%03d]","INT_fr",j+1,"INT_sr",j+1); 
			*outf<<null; out3n<<null;
			sprintf(null,"  %6s[G%03d]  %6s[G%03d]  %6s[G%03d]  %6s[G%03d]","INTcmt",j+1,"INTcmo",j+1,"INTimv",j+1,"INTcmr",j+1); 
			*outf<<null; out3n<<null; 
		}
		sprintf(null,"  %6s[G%03d]","INTtot",j+1); 
		*outf<<null; out3n<<null;
	}
	*outf<<endl; out3n<<endl;
	for(i=0;i<nused;i++) {
		sprintf(null," %13.4f",pwrfreq*i); 
		*outf<<null; out3n<<null;
		for(j=0;j<ngrp;j++) {
			k=0;
			sprintf(null," %13.4f %13.4f",Pvac[j][k][i].SEa,Pvac[j][k][i].SEb); 
			*outf<<null; out3n<<null;
			k=vangul;
			if(cmol) { 
				sprintf(null," %13.4f %13.4f",Pvac[j][k][i].SEa,Pvac[j][k][i].SEb); 
				*outf<<null; out3n<<null; 
			}
			for(k=0;k<vactype;k++) {
				sprintf(null," %13.4f",Pvac[j][k][i].min); 
				*outf<<null; out3n<<null;
			}
		}
		*outf<<endl; out3n<<endl;
	}

	sprintf(null,"\n Calculation of Thermodynamic Properties"); *outf<<null<<endl; outthermo<<null<<endl;
	if(cmol) vactype--;

	sprintf(null," %20s","property"); *outf<<null; outthermo<<null;
	for(j=0;j<ngrp;j++) {	
		if(show2pt) { 
			if(cmol) sprintf(null," %4s[G%03d] %4s[G%03d] %4s[G%03d] %4s[G%03d]","Tgas",j+1,"Tsol",j+1,"Rgas",j+1,"Rsol",j+1); 
			else sprintf(null," %4s[G%03d] %4s[G%03d]","Tgas",j+1,"Tsol",j+1); 
			*outf<<null; outthermo<<null; 
		}
		if(!cmol) { 
			sprintf(null," %4s[G%03d]","Tot",j+1); 
			*outf<<null; outthermo<<null; 
		} else {
			 sprintf(null," %4s[G%03d] %4s[G%03d] %4s[G%03d] %4s[G%03d]","Trns",j+1,"Rot",j+1,"Ivib",j+1,"Tot",j+1); 
			 *outf<<null; outthermo<<null;
		}
	}
	*outf<<endl; outthermo<<endl;

	if(cmol) {
		sprintf(null," %20s","nmolecules"); *outf<<null; outthermo<<null;
		for(j=0;j<ngrp;j++) {
			i=model->grp[j].nmol;
			k=0;
			if(show2pt) {
				if(cmol) { sprintf(null," %10d %10d %10d %10d",i,i,i,i); *outf<<null; outthermo<<null; }
				else { sprintf(null," %10d %10d",i,i); *outf<<null; outthermo<<null; }
			}
			for(;k<vactype;k++) {
				sprintf(null," %10d",i); 
				*outf<<null; outthermo<<null;
			}
		}
		*outf<<endl; outthermo<<endl;
	}

	sprintf(null," %20s","natom"); *outf<<null; outthermo<<null;
	for(j=0;j<ngrp;j++) {
		i=model->grp[j].natom;
		k=0;
		if(show2pt) {
			if(cmol) { sprintf(null," %10d %10d %10d %10d",i,i,i,i); *outf<<null; outthermo<<null; }
			else { sprintf(null," %10d %10d",i,i); *outf<<null; outthermo<<null; }
		}
		for(;k<vactype;k++) {
			sprintf(null," %10d",i); 
			*outf<<null; outthermo<<null;
		}
	}
	*outf<<endl; outthermo<<endl;
	print_thermo((char *)"dof",0,show2pt,outf,&outthermo);

	sprintf(null," %20s","temperature_____(K)"); *outf<<null; outthermo<<null;
	for(j=0;j<ngrp;j++) {
		k=0;
		if(show2pt) {
			if(cmol) sprintf(null," %10.2f %10.2f %10.2f %10.2f",vacT[j][k],vacT[j][k],vacT[j][k],vacT[j][k]); 
			else sprintf(null," %10.2f %10.2f",vacT[j][k],vacT[j][k]); 
			*outf<<null; outthermo<<null;
		}
		for(;k<vactype;k++) {
			sprintf(null," %10.2f",vacT[j][k]); 
			*outf<<null; outthermo<<null;
		}
	}
	*outf<<endl; outthermo<<endl;

	sprintf(null," %20s","pressure______(GPa)"); *outf<<null; outthermo<<null;
	for(j=0;j<ngrp;j++)	{
		k=0;
		if(show2pt) {
			if(cmol) sprintf(null," %10.2f %10.2f %10.2f %10.2f",vacP,vacP,vacP,vacP); 
			else sprintf(null," %10.2f %10.2f",vacP,vacP); 
			*outf<<null; outthermo<<null;
		}
		for(;k<vactype;k++) {
			sprintf(null," %10.2f",vacP); 
			*outf<<null; outthermo<<null;
		}
	}
	*outf<<endl; outthermo<<endl;

	sprintf(null," %20s","volume___(A^3)"); *outf<<null; outthermo<<null;
	for(j=0;j<ngrp;j++) {
		k=0;
		if(show2pt) {
			if(cmol) sprintf(null," %10.2f %10.2f %10.2f %10.2f",vacV[j],vacV[j],vacV[j],vacV[j]); 
			else sprintf(null," %10.2f %10.2f",vacV[j],vacV[j]); 
			*outf<<null; outthermo<<null;
		}
		for(;k<vactype;k++) {
			sprintf(null," %10.2f",vacV[j]);
			*outf<<null; outthermo<<null;
		}
	}
	*outf<<endl; outthermo<<endl;

	print_thermo((char *)"ZPE_(kJ/mol/SimBox)", 1,show2pt,outf,&outthermo);
	print_thermo((char *)"Emd_(kJ/mol/SimBox)",10,show2pt,outf,&outthermo);
	print_thermo((char *)"Eq__(kJ/mol/SimBox)", 2,show2pt,outf,&outthermo);
	print_thermo((char *)"Sq_(J/mol_K/SimBox)", 4,show2pt,outf,&outthermo);
	print_thermo((char *)"Aq__(kJ/mol/SimBox)", 6,show2pt,outf,&outthermo);
	print_thermo((char *)"Cvq(J/mol_K/SimBox)", 8,show2pt,outf,&outthermo);
	print_thermo((char *)"S(0)(cm/mol/SimBox)",12,show2pt,outf,&outthermo);

	sprintf(null," %20s","Diffus(cm2/s_,_1/s)"); *outf<<null; outthermo<<null;
	for(j=0;j<ngrp;j++) {
		if(show2pt) {
			k = 0;
			x=(Pvac[j][k][0].a*R*vacT[j][k])/(12.0*vlight*model->grp[j].mass)*1.0E5;
			sprintf(null," %10.4e %10.4e",x,0.0); *outf<<null; outthermo<<null;
			if (cmol) {
				k = 1;
				x=(Pvac[j][k][0].a*R*vacT[j][k])/(12.0*vlight*model->grp[j].mass)*1.0E5;
				sprintf(null," %10.4e %10.4e",x,0.0); *outf<<null; outthermo<<null;
			}
		}
		for(k=0;k<vactype;k++) {
			x=(thermo[12][j][k]*R*vacT[j][k])/(12.0*vlight*model->grp[j].mass)*1.0E5;
			sprintf(null," %10.4e",x); *outf<<null; outthermo<<null;
		}
	}
	*outf<<endl; outthermo<<endl;

	sprintf(null," %20s","fluidicity_________"); *outf<<null; outthermo<<null;
	for(j=0;j<ngrp;j++) {
		if(show2pt) {
			k=0;
			sprintf(null," %10.4e %10.4e",f2pt[k][j],0.0); *outf<<null; outthermo<<null;
			if (cmol) {
				k=1;
				sprintf(null," %10.4e %10.4e",f2pt[k][j],0.0); *outf<<null; outthermo<<null;
			}
		}
		for(k=0;k<vactype;k++) {
			x=0.0;
			if(k==vtrans) x=f2pt[0][j];
			else if(k==vangul) x=f2pt[1][j];
			sprintf(null," %10.4e",x); *outf<<null; outthermo<<null;
		}
	}
	*outf<<endl; outthermo<<endl;
 
	outvac.close();
	outpwr.close();
	out3n.close();
	outthermo.close();

	cout<<"Done"<<endl;
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::print_thermo(char * label, int i, int show2pt, ostream *outf, ostream *outs)
{
	int j,k;
	double x;
	char null[1024];

	x = 0.0;
	sprintf(null," %20s",label); *outf<<null; *outs<<null;
	for(j=0;j<ngrp;j++) {
		k=0;
		if (std::isnan(thermo[i][j][k]) || std::isinf(thermo[i][j][k])) thermo[i][j][k] = 0.0;
		if(show2pt) { 
			x = get_2pt_val(i,j,k);
			sprintf(null," %10.2f %10.2f",x,thermo[i][j][k]-x); *outf<<null;*outs<<null;
			if(cmol) {
				k=1;
				x = get_2pt_val(i,j,k);
				sprintf(null," %10.2f %10.2f",x,thermo[i][j][k]-x); *outf<<null;*outs<<null;
			}
		}
		k=0;
		for(;k<vactype;k++) {
			if (std::isnan(thermo[i][j][k]) || std::isinf(thermo[i][j][k])) thermo[i][j][k] = 0.0;
			sprintf(null," %10.2f",thermo[i][j][k]); *outf<<null; *outs<<null;
		}
	}
	*outf<<endl; *outs<<endl;
}

/* ---------------------------------------------------------------------- */
double ComputeVAC::get_2pt_val(int i, int j, int k)
{
	double x = 0.0;

	if(i ==0) x = hsdf[k][j];
	else if (i ==1) x = 0;
	else if (i ==2) x=hsdf[k][j]/thermo[0][j][k]*thermo[10][j][k];
	else if (i ==4) x=hsdf[k][j]*wsp[j];
	else if (i ==6) x=hsdf[k][j]/thermo[0][j][k]*thermo[10][j][k]-vacT[j][k]*hsdf[0][j]*wsp[j]*1E-3;
	else if (i ==8) x=hsdf[k][j]*wcvp[j];
	else if (i==10) x=hsdf[k][j];
	else if (i==12) x=Pvac[j][k][0].a;

	return x;
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::get_rot_temp(MODEL *model, int j, double *rT)
{
	int i,id;
	double tmp;

	for(i=0;i<model->grp[j].nmol;i++) {
		id=model->grp[j].mol[i];
		tmp=(h*h/(8.0*PI*PI*pI[id][0]*kb)); if(tmp<0) tmp=0; rT[0]+=tmp;
		tmp=(h*h/(8.0*PI*PI*pI[id][1]*kb)); if(tmp<0) tmp=0; rT[1]+=tmp;
		tmp=(h*h/(8.0*PI*PI*pI[id][2]*kb)); if(tmp<0) tmp=0; rT[2]+=tmp;
		if(tmp<0) tmp=0;
	}
	rT[0]/=model->grp[j].nmol;
	rT[1]/=model->grp[j].nmol;
	rT[2]/=model->grp[j].nmol;
	if((int)fabs(model->grp[j].linear)) rT[2]=-999;
}

/* ---------------------------------------------------------------------- */
void ComputeVAC::fix_missing_grp_mol_entry(MODEL *model)
{
	int i,j,id;

	for (i=0;i<model->ngrp;i++) {
		if(model->grp[i].nmol==0) {
			MOLECULE *tmol = new MOLECULE [nmol];
			for(j=0;j<nmol;j++) tmol[j] = model->mol[j];
			model->mol = NULL;
			model->mol = new MOLECULE [nmol+1];
			for(j=0;j<nmol;j++) model->mol[j] = tmol[j];
			//save inflated mol array
			model->mol[nmol].natom=model->grp[i].natom;
			model->mol[nmol].mass=model->grp[i].mass;
			model->mol[nmol].atm = new int [model->grp[i].natom];
			model->mol[nmol].atom = new ATOM* [model->grp[i].natom];
			for(j=0;j<model->grp[i].natom;j++) {
				id = model->grp[i].atm[j];
				model->mol[nmol].atom[j] = & model->atom[id];
			}
			//model->grp[i].mol = new int [1];
			model->grp[i].mol[0] = nmol;
			nmol++;
			model->nmol++;
			model->grp[i].nmol=1;
		}
	}
}
