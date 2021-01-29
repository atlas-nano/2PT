#ifndef TRJREADER_H
#define TRJREADER_H
#include <fstream>
#include <iostream>
#include "constant.h"
#include "property.h"
#include "structure.h"
#include "trj_header.h"

class TRJ //double precision for version 2010
{
public:

  CELL *cell;
  PROPERTY *prp;
  ATOM *atom;
  GROUP *grp;

  TRJheader header;

  double avetem;	//
  double  dstep;	//time step
  double  firstt;	//initial set temperature K
  double  finalt;	//final set temperature
  double  ec;		//constraint
  double  eu;		//user energy
  double  tint;		//total interanl
  double  tnb;		//total nonbond
  double  ea;		//trj average properties
  double  eba;
  double  eta;
  double  epa;
  double  eia;
  double  enba;
  double  eela;
  double  ehba;
  double  eca;
  double  eua;
  double  tinta;
  double  tnba;	//trj average properties
  double  totea;
  double  tkea;
  double  dum[12];	//unknown meanings
  double  duma[12];	//nuknown averaged properties

  int  iconmp;
  int imstep;
  int iconfs;
  int icstep;
  int ngrp;
  int natom;
  logical lvelwr;	//velocity output
  logical lfrcwr;	//force output

  double  pressur; //, pressura;	//pressure in GPa
  double  vol; //, vola;		//volume A3
  double  pvtot, pvtota;
  double  pvkin, pvkina;
  double  pvpot, pvpota;
  double  radgyr, radgyra;

  double  signose;
  double  zfrict;
  double  zprfrict;
  double  snose,snoseh,ssdot;
  double  qcanon;
  double  sigdyn[2];
  double  gamtmp;

  double  tcel,ucel,tcela,ucela;
  double  s2ra[6];//,s2r[6];
  double  s2rdot[6];
    
 int  natmcel;
 double  strsa[6];//,strs[6]; xx,yy,zz,xy,xz,yz
 double  extstrs,extstrsa; //! added after 220 but version called 210

 double  eabtota;
 double  eabvala;
 double  eabelha;
 double  eabnba;
 double  eabmisa;
 double  dltfaba;
 double  expprta;
 //double *x,*y,*z; //coordinates in Angstroms
 //double *velx,*vely,*velz; //velocities in Angstroms/ps

 int itrjformat;     //in trj format 0:binary 1:ascii
 int otrjformat;     //out trj format 0:binary 1:ascii
 ifstream intrj;     //input trj (velocity) stream
 ifstream incoord;   //input coordinate stream for xyz/amber/charmm trjs
 ifstream ineng;     //input energy stream for lammps trjs
 ifstream instr;     //input stress stream
 ifstream inatomeng; //input lammps atom energies
 logical has_atom_eng; //flag for whether we have atom energies

 unsigned long datasize,fcount,totframe;
 unsigned long location0,location1,location2;//after header,after 1st frame,end
 unsigned long ascfid,datasize2fid; //frame id, datasize upto frame fid
 unsigned long location0e,datasize2fide; //after header, datasize upto frame lf
 unsigned long location0s,datasize2fids; //after header, datasize upto frame lf

 TRJ();
 ~TRJ();
 int  init_trj(trjItem );
 int  init_c2bintrj(trjItem );
 int  init_c2asctrj(trjItem );
 int  init_usrtrj(trjItem ); //user
 int  init_xyztrj(trjItem ); //user
 int  init_lmptrj(trjItem );  //lammps
 int  init_ambertrj(trjItem ); //amber
 int  init_cpmdtrj(trjItem ); //cpmd
 int  rd_header(ATOM *);    //read header
 void rd_frame(int ); //read frame

 int  init_contentDouble();
 void ReadBinEng();
 void ReadBinFrame();
 void ReadAscFrame();
 void ReadCPMDFrame();
 void ReadUSRFrame();
 void ReadXYZFrame(ifstream *, int);
 void CleanTRJcontentDouble();
 void ReadCHARMMFrame();
 void ReadAMBERFrame(ifstream *,double, int, int, int);
 void ReadLMPFrame(int buflen, int iframe);
 void ReadLMPThermo();
 void rd_atom_eng(); //read atom energies
 void ReadLMPAtomEng(int buflen);

};

#endif
