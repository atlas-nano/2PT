#ifndef STATISTICS_H
#define STATISTICS_H
#include <stdlib.h>
#include "constant.h"

using namespace std;

class STAT
{
public:
  STAT() {
     max=bmax=-1e99; min=bmin=-max;
     sum=sum2=avg=std=bsize=sigma2=0;
     npt=nbin=nbinpt=nbinexceed=nXY=0;
     bin=NULL;
     binflag=0;
     dataflag=0;
     XYdata=NULL;
     lsfitflag=0;
  }
  ~STAT() {
     if(bin!=NULL) delete [] bin;
     if(XYdata!=NULL) delete [] XYdata;
  }
  double max;
  double min;
  double sum;
  double sum2;
  double avg;
  double std;
  double bmax;
  double bmin;
  double bsize;
  double sigma2;
  double *bin;
  double (*XYdata)[2]; //data for least square fit calc
  int    nbin;
  int    nbinpt;
  int    nbinexceed;
  int    npt;
  int    binflag;   //flag for histogram analysis
  int    nXY;
  int    dataflag;  //flag for storage of original data
  int    lsfitflag; //flag for least square fitting
  double Ax,Ay,Sxx,Sxy,Syy; //avg X,Y,
  double sigma2x,sigma2y,covxy,R2,s2,s; //variance,correlation coefficient, variance
  double a,b;   //Y= a + b X
  double SEa,SEb;   //standard errors in aXY and bXY
  int init();
  int reset();
  int cal_avg();
  int add_data(double);
  int add_data(double,double);
  int normalize(double);
  int add_Ydata(double);
  int init_bin();
  int reset_bin();
  int add_bin(double);
  int normalize_bin();
  int init_XYdata();
  int reset_XYdata();
  int add_XYdata(double,double);
  int lsfit_XYdata();
};

#endif
