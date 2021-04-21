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

#include <iostream>
#include <math.h>
#include "structure.h"
#include "utility.h"
#include "memory.h"

/* MOLECULE Functions */
MOLECULE::MOLECULE () 
{
  int i;

  atm=NULL;
  mass=chg=mu[0]=mu[1]=mu[2]=0;
  head=tail=-1;
  linear=0;
  atom=NULL;
  inertia=new double *[3];
  for(i=0;i<3;i++) {
    inertia[i]=new double [3];
  }
}

MOLECULE::~MOLECULE () 
{
  if(atm!=NULL) delete [] atm;
  if(atom!=NULL) {
    delete [] atom;
    atom=NULL;
  }
  if(inertia!=NULL) {
    delete [] inertia[0];
    delete [] inertia[1];
    delete [] inertia[2];
    delete [] inertia;
  }
}

int MOLECULE::cal_mass()
{
    int i;
    mass=0;
    for(i=0;i<natom;i++) mass += atom[i]->mass;
    return 0;
}

int MOLECULE::cal_chgmu()
{
    int i,k;
    chg=mu[0]=mu[1]=mu[2]=0;
    for(i=0;i<natom;i++) {
        chg += atom[i]->chg;
        for(k=0;k<3;k++) mu[k] += (atom[i]->chg)*(atom[i]->pv[k])*eA2debye;
    }
    return 0;
}

int MOLECULE::init_atom()
{
    if(natom==0) return 1;
    if(atm!=NULL) delete [] atm;
    atm= new int [natom];
    atom = new ATOM *[natom];
    return 0;
}

int MOLECULE::find_cm()
{
    int i,k;
    mass=pv[0]=pv[1]=pv[2]=0;
    for(i=0;i<natom;i++) {
        mass += atom[i]->mass;
        for(k=0;k<3;k++) pv[k] += (atom[i]->mass)*(atom[i]->pv[k]);
    }
    for(k=0;k<3;k++) pv[k]/=mass;
    return 0;
}

int MOLECULE::cal_vc()
{
// translation = center of mass motion
// M cmvel_k = sum_i sum_k ( mi vi_k )
//
// angular velocity
// omega = sum_i ( ri x vi )
// vrot = omega x ri
//
// vibration
// vib = vi - cmvel - vrot

    int i,k;
    double m;

    if(natom==1) {
      i=0;
      for(k=0;k<3;k++) {
         atom[i]->vt[k]=atom[i]->vel[k];
         atom[i]->vr[k]=atom[i]->vv[k]=anguv[k]=0;
         pI[k]=-1;
      }
      return 0;
    }

    for(i=0;i<natom;i++) for(k=0;k<3;k++)  atom[i]->vt[k]=atom[i]->vr[k]=atom[i]->vv[k]=0; 
     
    for(i=0;i<3;i++) {
       omega[i]=angmom[i]=0;
       for(k=0;k<3;k++) inertia[i][k]=0;
       pv[i]=vel[i]=0;
    }

    if(DEBUG==2) {
      printf("atomic coordinates (A)\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->pv[0],atom[i]->pv[1],atom[i]->pv[2]);
      printf("atomic velocities (A/ps)\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vel[0],atom[i]->vel[1],atom[i]->vel[2]);
    }
    
    //determine cm pv and vel
    for(i=0;i<natom;i++) {
        for(k=0;k<3;k++) {
            m=atom[i]->mass;
            pv[k]  += m*(atom[i]->pv[k]);
            vel[k] += m*(atom[i]->vel[k]);
        }
    }
    for(k=0;k<3;k++) { pv[k]/=mass; vel[k]/=mass;}

    double relpv[natom][3],relvel[natom][3];
    //set the center of mass to the origin
    for(i=0;i<natom;i++) {
        for(k=0;k<3;k++){
            relpv[i][k] = atom[i]->pv[k]-pv[k];
            relvel[i][k]= atom[i]->vel[k]-vel[k];
            atom[i]->vt[k]=vel[k];  //translational velocity
        }
    }

    if(DEBUG==2) {
      printf("translation velocities\n");
      printf("cm %8.4f %8.4f %8.4f\n",vel[0],vel[1],vel[2]);
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vt[0],atom[i]->vt[1],atom[i]->vt[2]);
    }

    //calculate angular momentum 
    for(i=0;i<natom;i++) { // L = sum_i mi (ri x vi) 
        m=atom[i]->mass;
        angmom[0] += m*(relpv[i][1]*relvel[i][2]-relpv[i][2]*relvel[i][1]);
        angmom[1] += m*(relpv[i][2]*relvel[i][0]-relpv[i][0]*relvel[i][2]);
        angmom[2] += m*(relpv[i][0]*relvel[i][1]-relpv[i][1]*relvel[i][0]);
    }

    //calculate inertia tensor of molecule
    for(i=0;i<natom;i++) {
        m=atom[i]->mass;
        inertia[0][0] += m*(pow(relpv[i][1],2.0) + pow(relpv[i][2],2.0));
        inertia[1][1] += m*(pow(relpv[i][0],2.0) + pow(relpv[i][2],2.0));
        inertia[2][2] += m*(pow(relpv[i][0],2.0) + pow(relpv[i][1],2.0));
        inertia[0][1] -= m* relpv[i][0] * relpv[i][1];
        inertia[0][2] -= m* relpv[i][0] * relpv[i][2];
        inertia[1][2] -= m* relpv[i][1] * relpv[i][2];
    }
    inertia[1][0] = inertia[0][1];
    inertia[2][0] = inertia[0][2];
    inertia[2][1] = inertia[1][2];

    if(DEBUG==2) {
      printf("\nInertia Tensor\n");
      for(i=0;i<3;i++) {
          for(k=0;k<3;k++) printf("%16.8lf",inertia[i][k]);
              printf("\n");
      }
      for(i=0;i<natom;i++) {
         printf("atom %d mass %f pv %f %f %f\n",i+1,atom[i]->mass,atom[i]->pv[0],atom[i]->pv[1],atom[i]->pv[2]);
      }
    }

    //calculate principle moments of inertia
    double evl[3],**evt,pomega[3];
    evt=new double *[3]; for(i=0;i<3;i++) evt[i]=new double [3];
    for(i=0;i<3;i++) evl[i]=evt[i][0]=evt[i][1]=evt[i][2]=pomega[i]=0;
    jacobi(inertia,3,evl,evt); //inertia tensor,dimension,eigenvalues,eigenvectors
    eigsrt(evl,evt,3); //sort eigenvalues
    if( linear ) evl[2]=0; //stlin 2010/08/04, added to correct for flexible CO2
    for(i=0;i<3;i++) pI[i]=evl[i]/(Na*1e23); //principle moment of inertia
    //for(i=0;i<3;i++) rotT[i]=h*h*Na*1e23/(8.0*PI*PI*evl[i]*kb); rotational temperatures

    if(DEBUG==2) {
      printf("\nPrinciple moments\n");
      for(i=0;i<3;i++) printf("%16.8lf vector %16.8lf %16.8lf %16.8lf\n",evl[i],evt[0][i],evt[1][i],evt[2][i]);
      printf("\n");
      printf("\nRotational temperatures\n");
      for(i=0;i<3;i++) printf("I %16.8lf (g/mol*A2) theta %16.8lf K\n",evl[i],h*h/(8.0*PI*PI*pI[i]*kb));
      printf("\n");
    }

    //calculate angular velocities along principle axises
    for(i=0;i<3;i++) pomega[i]=0;
    for(i=0;i<3;i++) for(k=0;k<3;k++) { if(evl[i]>0) pomega[i]+=angmom[k]*evt[k][i]/evl[i]; }

    //calculate inertia weighted angular velocities
    for(i=0;i<3;i++) {
        anguv[i]=omega[i]=0;
        for(k=0;k<3;k++) {
           if(evl[k]>0) {
              //angular velocity
              omega[i] += pomega[k]*evt[i][k];
              //angular vel weighted by principle moments of inertia 
              anguv[i] += pomega[k]*evt[i][k]*sqrt(evl[k]);
           }

        }
    }

    if(DEBUG==2) {
      printf("angular velocities along principle axises\n");
      printf("w %8.4f %8.4f %8.4f\n",pomega[0],pomega[1],pomega[2]);
      for(i=0;i<3;i++) printf("%8.4f v %8.4f %8.4f %8.4f\n",pomega[i],pomega[i]*evt[0][i],pomega[i]*evt[1][i],pomega[i]*evt[2][i]);
    }

    //calculate velocity due to rotation
    for(i=0;i<natom;i++) { //vr = w x r
        atom[i]->vr[0]=(omega[1]*relpv[i][2]-omega[2]*relpv[i][1]);
        atom[i]->vr[1]=(omega[2]*relpv[i][0]-omega[0]*relpv[i][2]);
        atom[i]->vr[2]=(omega[0]*relpv[i][1]-omega[1]*relpv[i][0]);
    }

    if(DEBUG==2) {
      printf("angular velocities\n");
      printf("w %8.4f %8.4f %8.4f\n",omega[0],omega[1],omega[2]);
      printf("weighted angular velocities\n");
      printf("w %8.4f %8.4f %8.4f\n",anguv[0],anguv[1],anguv[2]);
      printf("rotational velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vr[0],atom[i]->vr[1],atom[i]->vr[2]);
    }

    //calculate vibrational velocities
    for(i=0;i<natom;i++) {
       for(k=0; k<3; k++) atom[i]->vv[k]= relvel[i][k]-(atom[i]->vr[k]);
    }

    if(DEBUG==2) {
      printf("vibration velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vv[0],atom[i]->vv[1],atom[i]->vv[2]);
    }

    if(DEBUG==2) {
      printf("translational velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vt[0],atom[i]->vt[1],atom[i]->vt[2]);
      printf("rotational velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vr[0],atom[i]->vr[1],atom[i]->vr[2]);
      printf("vibrational velocities\n");
      for(i=0;i<natom;i++) printf("%d %8.4f %8.4f %8.4f\n",i+1,atom[i]->vv[0],atom[i]->vv[1],atom[i]->vv[2]);
    }
    for(i=0;i<3;i++) delete [] evt[i];
    delete [] evt;

    return 0;
}

/* GROUP functions */
GROUP::GROUP() 
{
  natom=0;
  nmol=0;
  constraint=0;
  rotsym=1;
  linear=0;
  chg=mass=0;
  vol = 0.0;
  atm=NULL;
  mol=NULL;
  color=-1;
  atomsize=bondsize=labelsize=1;
  mu[0]=mu[1]=mu[2]=0;
}

GROUP::~GROUP() 
{
  natom=0;
  if(atm!=NULL) delete [] atm;
  if(mol!=NULL) delete [] mol;
}

int GROUP::init_atom()
{
   if(natom==0) return 1;
   if(atm!=NULL) delete [] atm;
   atm= new int [natom];
   return 0;
}

/* CELL Functions */
CELL::CELL () 
{
  H[0][0]=H[0][1]=H[0][2]=0;
  H[1][0]=H[1][1]=H[1][2]=0;
  H[2][0]=H[2][1]=H[2][2]=0;
  la=lb=lc=alpha=beta=gamma=volume=0;
  a[0]=a[1]=a[2]=0;
  b[0]=b[1]=b[2]=0;
  c[0]=c[1]=c[2]=0;
  o[0]=o[1]=o[2]=0;
  Hinv[0][0]=Hinv[0][1]=Hinv[0][2]=0;
  Hinv[1][0]=Hinv[1][1]=Hinv[1][2]=0;
  Hinv[2][0]=Hinv[2][1]=Hinv[2][2]=0;
  cb[0]=cb[1]=cb[2]=0;
}

CELL::~CELL () 
{
  ;
}

int CELL::s2r2H(double *s2r)
{
  H[0][0]=s2r[0];
  H[1][1]=s2r[1];
  H[2][2]=s2r[2];
  H[1][0]=s2r[5];
  H[2][0]=s2r[4];
  H[2][1]=s2r[3];
  H[0][1]=0;
  H[0][2]=0;
  H[1][2]=0;

  return 0;
}

int CELL::cal_Hinv()
{
  //inverse of H matrix
  Hinv[0][0]=1/H[0][0];
  Hinv[1][1]=1/H[1][1];
  Hinv[2][2]=1/H[2][2];
  Hinv[1][0]=-H[1][0]/H[0][0]/H[1][1];
  Hinv[2][0]=(H[1][0]*H[2][1]-H[1][1]*H[2][0])/H[0][0]/H[1][1]/H[2][2];
  Hinv[2][1]=-H[2][1]/H[1][1]/H[2][2];
  Hinv[0][1]=0;
  Hinv[1][2]=0;
  Hinv[0][2]=0;

  return 0;
}

int CELL::H2others()
{
  //  the H matrix
  //  h00  h01 h02     a1 b1 c1   h00   0   0
  //  h10  h11 h12     a2 b2 c2   h10 h11   0
  //  h20  h21 h22     a3 b3 c3   h20 h21 h22

   double fac=180.0/acos(-1.0);
   a[0]=H[0][0]; a[1]=H[1][0]; a[2]=H[2][0];
   b[0]=0;       b[1]=H[1][1]; b[2]=H[2][1];
   c[0]=0;       c[1]=0;       c[2]=H[2][2];
   la = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
   lb = sqrt(b[1]*b[1] + b[2]*b[2]);
   lc = c[2];
   alpha = fac * acos( b[2] / lb );
   beta =  fac * acos( a[2] / la );
   gamma = fac * acos( (a[1]*b[1] + a[2]*b[2])/la/lb );

   //center of box
   cb[0]=0.5*(a[0]+b[0]+c[0]);
   cb[1]=0.5*(a[1]+b[1]+c[1]);
   cb[2]=0.5*(a[2]+b[2]+c[2]);

  //inverse of H matrix
  cal_Hinv();

   cal_volume();
   return 0;
}

int CELL::labc2H()
{
   //given la,lb,lc,alpha,beta,gamma, get H2
   double aa,bb,cc;
   double fac=acos(-1.0)/180.0;
   aa=alpha*fac; bb=beta*fac; cc=gamma*fac;
   //place c along Z
   c[0]=0; c[1]=0; c[2]=lc;
   //place b in the yz plane
   b[0]=0; b[1]=lb*sin(aa); b[2]=lb*cos(aa);
   //find vector for a
   a[2]=la*cos(bb);
   a[1]=(la*lb*cos(cc)-a[2]*b[2])/b[1];  //la*(cos(cc)-cos(aa)*cos(bb))/sin(aa);
   a[0]=sqrt(la*la-a[1]*a[1]-a[2]*a[2]);

   H[0][0]=a[0]; H[1][0]=a[1]; H[2][0]=a[2];
   H[0][1]=b[0]; H[1][1]=b[1]; H[2][1]=b[2];
   H[0][2]=c[0]; H[1][2]=c[1]; H[2][2]=c[2];

   return 0;
}

int CELL::cal_volume()
{
   volume=0;
   volume=fabs(H[0][0]*H[1][1]*H[2][2]);
   return 0;
}

/* BOND functions */
BOND::BOND () 
{
  ffid=id=-1;
  ffp=NULL;
  lbd=-1;
  len=len0=-1;
}

BOND::~BOND () 
{
  ffp=NULL;     
}

int BOND::cal_len()
{
    len=0;
    double v[3]; //v = R(C1-C2)
    int k;
    double ip=0;
    for(k=0;k<3;k++) {
       v[k]=(atom[0]->pv[k])-(atom[1]->pv[k]);
       ip+=(v[k]*v[k]);
    }
    len=sqrt(ip);
    return 0;
}
