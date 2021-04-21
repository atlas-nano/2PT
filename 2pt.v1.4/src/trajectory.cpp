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

//read atom trajectories

#include <iostream>
#include <vector>
#include <sstream>
#include "trajectory.h"
//#include <sys/mman.h>
//#include <sys/stat.h>

/* ---------------------------------------------------------------------- */
TRAJECTORY::TRAJECTORY() 
{
    strj=NULL;
    iframe=NULL;
    fframe=NULL;
}

/* ---------------------------------------------------------------------- */
TRAJECTORY::~TRAJECTORY() 
{
    if(strj!=NULL) delete [] strj;
    if(iframe!=NULL) delete [] iframe;
    if(fframe!=NULL) delete [] fframe;
}

/* ---------------------------------------------------------------------- */
int TRAJECTORY::init_trj(int ntrj,trjItem *trjfname,int iformat,int oformat,CELL *in_cell,ATOM *in_atom,PROPERTY *in_prp, int natom)
{
    int i;
    int filestatus=0;

    ntrjf=ntrj;
    itrjf=iformat;
    otrjf=oformat;

    if(ntrjf==0) return filenotopen;
    if(strj!=NULL) delete [] strj;
    if(iframe!=NULL) delete [] iframe;
    if(fframe!=NULL) delete [] fframe;

    strj=new TRJ [ntrjf];
    iframe=new int [ntrjf];
    fframe=new int [ntrjf];
    totframe=0;

    for(i=0;i<ntrjf;i++) {
	strj[i].itrjformat=itrjf;
        strj[i].otrjformat=otrjf;
        strj[i].cell=in_cell;
        strj[i].prp=in_prp;
        strj[i].atom=in_atom;
        if(iformat==ambertrj) strj[i].header.totmov = strj[i].header.nmovatm = natom; //amber traj fix
        filestatus+=strj[i].init_trj(trjfname[i]);
        totframe+=strj[i].totframe;
        fframe[i]=totframe;
    }
    if(filestatus % filenotopen) return filenotopen;
    if(iformat==usrtrj) {
    } else if(iformat==xyztrj) {
    } else if(iformat==lmptrj) {
	for(i=1;i<ntrjf;i++) {
	    fframe[i]-= i;
       }
    } else {
       //correction for lammps trj (remove the 1st frame for 2nd,3rd,..trj)
       for(i=1;i<ntrjf;i++) {
          //fframe[i]-= i;
	}
    }

    totframe -= (ntrjf-1);
    return itrjf;
}

/* ---------------------------------------------------------------------- */
void TRAJECTORY::rd_frame(int ith)
{
    int i;

    cframe=ith;
    //find current trj file
    for(i=0;i<ntrjf;i++) {
	if( fframe[i]>=ith) break;
    }
    ctrjf=i;
    if(ctrjf==0)  strj[ctrjf].rd_frame(ith);
    else if(ctrjf==lmptrj) strj[ctrjf].rd_frame(ith);
    else if(ctrjf==xyztrj) strj[ctrjf].rd_frame(ith);
    else if(ctrjf==cpmdtrj) strj[ctrjf].rd_frame(ith);
    else strj[ctrjf].rd_frame(ith-fframe[ctrjf-1]+1); //adding one for lammps trj
}

/* ---------------------------------------------------------------------- */
void TRAJECTORY::rd_atom_eng()
{
    int i;

    for(i=0;i<ntrjf;i++)
	if(strj[i].has_atom_eng) 
	    strj[i].rd_atom_eng();
}
