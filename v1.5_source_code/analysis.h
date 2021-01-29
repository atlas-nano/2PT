#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include "constant.h"
#include "statistics.h"
#include "model.h"
#include <vector>
#include "compute.h"

using namespace std;

class ANALYSIS {
  public:
   ANALYSIS();
   ~ANALYSIS();
   void doit(MODEL *);

  private:
   MODEL *model;

   int nframe;
   int stat_nXY;   //number of data points

   int ncomputes;
   vector<Compute *> compute_list;

   ofstream trjf;      //output trj 

   void setup();
   void ana_init();
   void ana_ana();
   void ana_report(ostream *); //output analysis results

};

#endif
