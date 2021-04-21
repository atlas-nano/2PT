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

//top level routines for program
#include <iostream>
#include "model.h"
#include "analysis.h"
#include "timing.h"
#include "utility.h"

int main(int argc, char *argv[])
{

  int err;

  if(argc<2) {
     cout<<" Usage: "<<argv[0]<<" Control.in"<<endl;
     return 1;
  }

  copyright(&cout);
  MODEL model;
  ANALYSIS analysis;
  TIME time;

  if(model.ctl.rd_ctl(argv[1])) {
    cout<<" Error reading file"<<endl;
    return 1;
  }
  cout<<"Done"<<endl;

  model.ctl.out_ctl(&cout); //print options to screen
  err = model.rd_strt(); //read molecular structure
  if(err) return err;
  model.init_trj(); //get trj file info
  analysis.doit(&model); //do thermodynamic analysis
  //time.gettimetotal(&cout);

  return 0;
}

