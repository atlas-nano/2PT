#ifndef TRAJECTORY_H
#define TRAJECTORY_H
#include <fstream>
#include <iostream>
#include "constant.h"
#include "property.h"
#include "trj_reader.h"

class TRAJECTORY //multiple trjs
{
public:
    TRAJECTORY();
    ~TRAJECTORY(); 
    TRJ *strj;

    int itrjf; //in trj format
    int otrjf; //out trj format
    int ntrjf; //number of trj files
    int ctrjf;  //current trj file
    int cframe; //current frame id
    int *iframe;      
    int *fframe;      
    int totframe;

    int  init_trj(int, trjItem *,int,int,CELL*,ATOM*,PROPERTY*,int);
    void rd_frame(int );
    void rd_atom_eng(); //read atom energies
};

#endif
