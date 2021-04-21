#ifndef TRJHEADER_H
#define TRJHEADER_H
#include <fstream>
#include <iostream>
#include "constant.h"
#include "property.h"
#include "structure.h"

const int MAXFILE=24;

class byteoffset {
public:
    long int xyz; //trajectory file
    long int coord; //coordinates file
    long int eng; //energy file
    logical vamberbox; //box information store at end of amber frame
    logical camberbox; //box information store at end of amber frame
};

class TRJheader
{
public:
    byteoffset *byte_offset; //holds lammps offset info
    double *** lmp_data_index; //hold pointer to field for each atom
    char hdr[5];	//Header
    int icntrl[20];	//Control (icntrl[0]=trj file version, default=2010)
    int ntrjti;		//number of comments
    char trjtic[10][80];	//comments
    int neexti;		//number of EEX comments
    char eextic[10][80];	//EEX comments
    logical period;		//Periodicity (1,2,3 =1D,2D,3D)
    logical molxtl;		//Molecular Crystal(unknown usage)
    logical lcanon;		//Canonical (NVT)
    logical defcel;		//Cell definition (make it true to be safe)
    logical prtthrm;	//Perturbation Theory
    logical lnose;		//NoseorHoover
    logical lnpecan;	//NPT Canonical
    logical ltmpdamp;	//Temperature damping
    int  nflusd;		//number of files unsed (make it 1)
    int  mvatmpfu[MAXFILE];	//Moveable atoms (must give mvatmpfu[0])
    int  natmpfu[MAXFILE];	//Total atoms (must give natmpfu[0])
    char  decusd[MAXFILE][8];//Descriptor
    int  totmov;		//Total moveable atoms (totmov=mvatmpfu[0] when nflusd=1)
    int  *mvatmofst; 	//movable atom ids (starts from 1)
    int leexti;		//length of EEX title
    char eextit[80];	//EEX title
    int lparti;		//length of parameter file
    char partit[80];	//parameter file
    int natom;		//total number of atoms		
    int nmovatm;
    int movatm1;
    int movatmn;
    int version;
    int lmp_cell_format; //lammps cell format

    int xyzfreq,nxyz;  //coordinate dump frequency, skips when reading traj
    int velfreq,nvel;  //velocity dump frequency
    int engfreq,neng;  //energy dump frequency
    int strfreq,nstr;  //stress dump frequency
    double timestep; //simulation time step in ps
    int totaccustep;   //total accumulated steps (needed for CPMD trj)
    int lmptmp; //temp storage for lammps
    int lmp_data_len; //length of dump items
    int xyz_flag,vel_flag,stress_flag,image_flag,force_flag,charge_flag,eng_flag;

    TRJheader();
    ~TRJheader();
    int init_header();
    void CleanTRJheader(int);
    int ReadBinHeader(ifstream *);
    int ReadAscHeader(ifstream *);
    int ReadCPMDHeader(ifstream *,ifstream *,ifstream *);
    int ReadUSRHeader(ifstream *);
    int ReadXYZHeader(ifstream *,int);
    int ReadLMPHeader(ifstream *,ifstream *,ATOM *);
    int ReadAMBERHeader(ifstream *, int);
    int ReadCHARMMHeader(ifstream *,ifstream *);
private:
    void SetDefLMPFormat(ATOM *,int); //set the default lammps dump format
};

#endif
