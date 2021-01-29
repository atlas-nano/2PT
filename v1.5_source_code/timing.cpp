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
#include <fstream>
#include "timing.h"

void TIME::settimezero()
{
    myclock();
    timelast=timezero=timecurrent;
    clocklast=clockzero=clockcurrent;
}   

void TIME::myclock()
{
    clockcurrent=clock();
    timecurrent=time(NULL);
}

void TIME::gettimetotal(ostream *outf)
{
    myclock();
    outtime(timezero,timecurrent,clockzero,clockcurrent,(char *)"Total",outf);
}

void TIME::outtime(time_t itime,time_t ftime,clock_t iclock,clock_t fclock,char *str,ostream *outf)
{
    time_t difft;
    clock_t diffc;

    difft=ftime-itime;
    diffc=fclock-iclock;

    day=  difft/cday; 
    hr = (difft%cday)/chr;
    min=((difft%cday)%chr)/cmin;
    sec=(((difft%cday)%chr)%cmin);
    dsec=sec+(diffc%csec)/(double)csec;
    *outf<<str<<" time used "<<day<<" day "<<hr<<" hrs "<<min<<" min "<<dsec<<" sec"<<endl;
}
