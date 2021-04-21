/*********************************************************************************
*					 Two-Phase Thermodynamics (2PT) Program					 *
*						  Shiang-Tai Lin (stlin@ntw.edu.tw)					 *
* Department of Chemical Engineering, National Tiawan University, Taipei, Taiwan *
*					 Prabal K. Maiti (maiti@physics.iisc.ernet.in)			  *
*  Department of Physics, Indian Institute of Science, Bangalore, India, 560012  *
*						  Tod A Pascal (tpascal@wag.caltech.edu)				*
*									 and										*
*						William A Goddard III (wag@wag.caltech.edu)			 *
*		 Materials and Process Simulation Center, Caltech, Pasadena, CA USA	 *
*								Copyright (c) 2010							  *
*								All rights Reserved							 *
*							 *************************						  *
*									  Cite:									 *
* S.T. Lin, M. Blanco and W.A. Goddard, J. Chem. Phys., 2003, 119, 11792-11805   *
* S.T. Lin, P.K. Maiti and W.A. Goddard, J. Phys. Chem. B., 2010, 114, 8191-8198 *
* T.A. Pascal, S.T. Lin and W.A. Goddard, PCCP, 2011, 13(1), 169-181			 *
***********************************************************************************/

//top level control routine for analysis

#include <iostream>
#include "constant.h"
#include "utility.h"
#include "analysis.h"
#include "compute_vac.h"

/* ---------------------------------------------------------------------- */
ANALYSIS::ANALYSIS () {  
	nframe=0;
	stat_nXY=0;
}

/* ---------------------------------------------------------------------- */
ANALYSIS::~ANALYSIS () { 

	compute_list.clear();
}

/* ---------------------------------------------------------------------- */
void ANALYSIS::setup ()
{
	Compute * compute_pointer;

	if(model->ctl.ana_vac_flag) {
		compute_pointer = new ComputeVAC;
		compute_list.push_back( compute_pointer );
	}

	if(compute_list.size() == 0) {
		cerr<<" Nothing to analyze!"<<endl;
		exit(1);
	}
}

/* ---------------------------------------------------------------------- */
void ANALYSIS::doit(MODEL *in_model)
{
	int i;
	char filename[1024];

	model=in_model;
	if(model->ctl.ana_flag==0) return;
	if(model->ctl.in_trj_flag==filenotopen) return;

	i = 0;
	setup();
	ana_init();

	cout<<"Reading trajectory file"<<endl;
	int tot_frames = int((model->ctl.ana_fframe-model->ctl.ana_iframe+1.0)/model->ctl.ana_sframe);
	nframe = 0;
	for(i=model->ctl.ana_iframe;i<=model->ctl.ana_fframe;i+=model->ctl.ana_sframe) {
		cout<<" Reading frame "<<nframe<<"/"<<tot_frames<<"\xd";
		fflush(stdout);
		model->trj.rd_frame(i);
		ana_ana();
		nframe++;
	}
	cout<<" Reading frame "<<model->ctl.ana_fframe<<"/"<<model->ctl.ana_fframe<<endl;

	strcpy(filename,model->ctl.ana_out);
	strcat(filename,".out.log");
	ofstream out(filename);
	if ( out ) clog.rdbuf(out.rdbuf());
	else cerr << "Error! Cannot open log file. Printing all data to screen"<<endl;

	ana_report(&clog);

}

/* ---------------------------------------------------------------------- */
void ANALYSIS::ana_init()
{
	
	unsigned int j;

	nframe=0;
	stat_nXY=(model->ctl.ana_fframe-model->ctl.ana_iframe)/model->ctl.ana_sframe+1;
	model->rd_grp(model->ctl.in_grp);

	for(j=0;j<compute_list.size();j++) compute_list[j]->init(model);
}

/* ---------------------------------------------------------------------- */
void ANALYSIS::ana_ana()
{
	unsigned int i;
	for(i=0;i<compute_list.size();i++) compute_list[i]->rd_frame(model,nframe);

}

/* ---------------------------------------------------------------------- */
void ANALYSIS::ana_report(ostream *outf)
{
	unsigned int i;
	for(i=0;i<compute_list.size();i++) compute_list[i]->report(outf,model);

}

