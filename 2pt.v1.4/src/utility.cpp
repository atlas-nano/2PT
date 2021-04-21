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

/* General routines, some from Numerical Recipes (slightly modified). */
#include "utility.h"
#include <math.h>

/* ---------------------------------------------------------------------- */
void ROTATE(double **a,int i,int j,int k,int l,double *tau,double *s)
{
	double g,h;

	g=a[i][j];
	h=a[k][l];
	a[i][j]=g-(*s)*(h+g*(*tau));
	a[k][l]=h+(*s)*(g-h*(*tau));
}

/* ---------------------------------------------------------------------- */
void jacobi(double **a, int n, double *d, double **v)
{
	int j,iq,ip,i;
	int nrot=0;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b=new double [n];
	z=new double [n];
	for (ip=0;ip<n;ip++) {
		v[ip][ip]=1.0;
		for (iq=ip+1;iq<n;iq++) v[ip][iq]=v[iq][ip]=0.0;
	}
	for (ip=0;ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	for (i=0;i<50;i++) {
		sm=0.0;
		for (ip=0;ip<n-1;ip++) 
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		if (sm == 0.0) {
			delete [] z;
			delete [] b;
			break;
		}
		if (i < 3) tresh=0.2*sm/(n*n);
		else tresh=0.0;
		for (ip=0;ip<n-1;ip++) {
			for (iq=ip+1;iq<n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 3 && (fabs(d[ip])+g) == fabs(d[ip])
					&& (fabs(d[iq])+g) == fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((fabs(h)+g) == fabs(h)) t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0;j<ip;j++)
						ROTATE(a,j,ip,j,iq,&tau,&s);
					for (j=ip+1;j<iq;j++)
						ROTATE(a,ip,j,j,iq,&tau,&s);
					for (j=iq+1;j<n;j++)
						ROTATE(a,ip,j,iq,j,&tau,&s);
					for (j=0;j<n;j++) 
						ROTATE(v,j,ip,j,iq,&tau,&s);
					++(nrot);
				}
			}
		}
		for (ip=0;ip<n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	if(sm != 0.0) {
		printf("Too many iterations in routine jacobi %d ",nrot);
		if(DEBUG) {
			printf("original matrix\n");
			for(i=0;i<n;i++) {
				for(j=0;j<n;j++) printf("%8.4e ",a[i][j]);
					printf("\n");
			}
		}
		delete [] z;
		delete [] b;
	}
}

/* ---------------------------------------------------------------------- */
void eigsrt(double *d, double **v, int n)
{
	int k,j,i;

	for (i=0;i<n-1;i++) {
		for (j=i+1;j<n;j++) {
			if (d[j] > d[i]) {
				SWAP(&d[j],&d[i]);
				for(k=0;k<n;k++) SWAP(&v[k][i],&v[k][j]);
			}
		}
	}
}

/* ---------------------------------------------------------------------- 
 The following subroutines calculates the integration of weighting
 function for entropy from 0+ to 1/2 fmin. The algorithm is based on 
 mathematica result:
 w[u_] := (u/ (Exp[u] - 1) - Log[1 - Exp[-u]] )
 Integrate[  w[u], {u, 0, x}] =
 (Pi^2)/2 -x^2 -x Log[1-Exp[-x]] + 2 PolyLog[2,Exp[-x]])              */
double scweighting(double upper)
{
	if(upper==0 || isnan(upper)) return 0;
	return (2*upper-upper*log(upper))/upper;
}

/* ---------------------------------------------------------------------- */
double sqweighting(double upper)
{
	if(upper==0 || isnan(upper)) return 0;
	double pi=3.1415926535897932385;
	double wsq= pi*pi/3.0 -upper*upper +upper*log(-1+exp(upper))
				- 2*polylog(2.0,exp(-upper));

	return wsq/upper;
}

/* ---------------------------------------------------------------------- */
double polylog(double n, double z)
{
	int k;
	double sum,sum_old;

	k=1;
	sum=0.0;
	sum_old=0.0;

	do 
	{
		sum_old=sum;
		sum+= pow(z,k)/pow(k,n);
		k++;
	} while (fabs(sum-sum_old)!=0);

	return sum;
}

/* ---------------------------------------------------------------------- */
double HSDF(double *pwr,int lenth,double fmin,int nmol,double fract_f)
{
	double hsdf,tmpg,tmps;
	int j;

	hsdf=tmpg=tmps=0;
	for(j=0;j<lenth;j++) {
		twoPT(&tmpg,&tmps,pwr[0],pwr[j],fmin*j,nmol,fract_f);
		if(j==0 || j==lenth-1) { hsdf += tmpg*fmin*0.5;}
		else { hsdf += tmpg*fmin; }
	}
	return hsdf;
}

/* ---------------------------------------------------------------------- */
double TTDF(double *pwr,int lenth,double fmin)
{
	double ttdf;
	int j;

	ttdf=0;
	for(j=0;j<lenth;j++) {
		if(j==0 || j==lenth-1) { ttdf += pwr[j]*fmin*0.5;}
		else { ttdf += pwr[j]*fmin; }
	}
	return ttdf;
}

/* ---------------------------------------------------------------------- */
void HSweighting(double *wep,double *wap,double *wcvp,double *wsehs,double *wsp,double *wspd,double y,double mass,double nmol,double transT,double rotT,double volume,double *wsr,double *war,double *rT,double rs)
{
	*wep=*wap=*wcvp=*wsehs=*wsp=*wspd=*wsr=*war=0;
	//Ideal gas
	if(y<=0.74) { //packing fraction too large, ignore hard sphere correction
		*wsp=  5.0/2.0 + log(pow(2.0*PI*mass*1e-3/Na*kb*transT/h/h,1.5)*volume*1e-30/(nmol)); //ideal gas - indistinguishable particles 
		*wspd= 3.0/2.0 + log(pow(2.0*PI*mass*1e-3/Na*kb*transT/h/h,1.5)*volume*1e-30); //ideal gas - distinguishable particles 

		//Carnahan-Starling
		*wsehs=log((1.0+y+y*y-y*y*y)/pow(1.0-y,3.0))+y*(3.0*y-4.0)/pow(1.0-y,2.0); //eqn 36 in original paper

		//Hard sphere gas
		*wsp  = (*wsp  + *wsehs)/3.0; //eqn 35b
		*wspd = (*wspd + *wsehs)/3.0;
		*wep  = 0.5;
		*wap  = *wep-*wsp; /*wa=we-ws, no need of factor T in front of ws*/
		*wcvp = 0.5;

		//rigid rotor
		if(rT[0]>0.0001) { //for none-atomic molecules
			if(rT[2]<0) { //linear
				*wsr  = (1.0 + log(rotT/sqrt(rT[0]*rT[1])/rs))/2.0;
			} else { //non-linear
				*wsr  = (3.0/2.0 + log(sqrt(PI*rotT*rotT*rotT/(rT[0]*rT[1]*rT[2]))/rs))/3.0;
			}

			*war  = *wep-*wsr;
		} else *wsr = *war = 0;
	}
}

/* ---------------------------------------------------------------------- */
void twoPT(double *tmpg,double *tmps,double s0,double sv,double v,int nmol,double fract_f)
{
	if(fract_f==0) { *tmps=sv; *tmpg=0; }
	else if(v==0) { *tmpg=s0; *tmps=0; }
	else {
		*tmpg=s0/(1.0+pow(PI*s0*v/(6.0*nmol*fract_f),2.0)); //eqn 24
		//*tmpg=3*fract_f*nmol*(s0/(1.0+pow(PI*s0*v/2,2.0)));
		if(*tmpg>sv) *tmpg=sv;
		*tmps=sv-*tmpg;
	}
}

/* ---------------------------------------------------------------------- 
   determine fraction factor f from constant K                            */
double search2PT(double K)
{
	double P,fold,fnew,dPdf,tol;
	int count;

	fold=0.0;
	fnew=0.7293*pow(K,0.5727); /*initial guess*/
	if(fnew>0.5) fnew=0.5;
	tol=1e-10;
	count=0;

	while( fabs(fnew-fold)>tol && count<999) {
		dPdf=0.0;
		fold=fnew;
		P= 2.0*pow(K,-4.5)*pow(fnew,7.5)-6.0*pow(K,-3.0)*pow(fnew,5.0)-pow(K,-1.5)*pow(fnew,3.5)
			+ 6.0*pow(K,-1.5)*pow(fnew,2.5)+2.0*fnew-2;
		dPdf = 15.0*pow(K,-4.5)*pow(fnew,6.5)-30.0*pow(K,-3.0)*pow(fnew,4.0)-3.5*pow(K,-1.5)*pow(fnew,2.5)
			+ 15.0*pow(K,-1.5)*pow(fnew,1.5)+2.0;
		fnew=fold-(P)/dPdf;
		count++;
	}
	return fnew;
}

/* ---------------------------------------------------------------------- */
ele * loadelements() 
{
	ele *element = new ele [110];
	element[1].symbol = "h";	element[1].mass =  1.00790;	element[1].name = "Hydrogen";
	element[2].symbol = "he";	element[2].mass =  4.00260;	element[2].name = "Helium";
	element[3].symbol = "li";	element[3].mass =  6.94100;	element[3].name = "Lithium";
	element[4].symbol = "be";	element[4].mass =  9.01220;	element[4].name = "Beryllium";
	element[5].symbol = "b";	element[5].mass = 10.81100;	element[5].name = "Boron";
	element[6].symbol = "c";	element[6].mass = 12.01070;	element[6].name = "Carbon";
	element[7].symbol = "n";	element[7].mass = 14.00670;	element[7].name = "Nitrogen";
	element[8].symbol = "o";	element[8].mass = 15.99940;	element[8].name = "Oxygen";
	element[9].symbol = "f";	element[9].mass = 18.99840;	element[9].name = "Fluorine";
	element[10].symbol = "ne";	element[10].mass = 20.17970;	element[10].name = "Neon";
	element[11].symbol = "na";	element[11].mass = 22.98970;	element[11].name = "Sodium";
	element[12].symbol = "mg";	element[12].mass = 24.30500;	element[12].name = "Magnesium";
	element[13].symbol = "al";	element[13].mass = 26.98150;	element[13].name = "Aluminum";
	element[14].symbol = "si";	element[14].mass = 28.08550;	element[14].name = "Silicon";
	element[15].symbol = "p";	element[15].mass = 30.97380;	element[15].name = "Phosphorus";
	element[16].symbol = "s";	element[16].mass = 32.06500;	element[16].name = "Sulfur";
	element[17].symbol = "cl";	element[17].mass = 35.45300;	element[17].name = "Chlorine";
	element[19].symbol = "k";	element[19].mass = 39.09830;	element[19].name = "Potassium";
	element[18].symbol = "ar";	element[18].mass = 39.94800;	element[18].name = "Argon";
	element[20].symbol = "ca";	element[20].mass = 40.07800;	element[20].name = "Calcium";
	element[21].symbol = "sc";	element[21].mass = 44.95590;	element[21].name = "Scandium";
	element[22].symbol = "ti";	element[22].mass = 47.86700;	element[22].name = "Titanium";
	element[23].symbol = "v";	element[23].mass = 50.94150;	element[23].name = "Vanadium";
	element[24].symbol = "cr";	element[24].mass = 51.99610;	element[24].name = "Chromium";
	element[25].symbol = "mn";	element[25].mass = 54.93800;	element[25].name = "Manganese";
	element[26].symbol = "fe";	element[26].mass = 55.84500;	element[26].name = "Iron";
	element[28].symbol = "ni";	element[28].mass = 58.69340;	element[28].name = "Nickel";
	element[27].symbol = "co";	element[27].mass = 58.93320;	element[27].name = "Cobalt";
	element[29].symbol = "cu";	element[29].mass = 63.54600;	element[29].name = "Copper";
	element[30].symbol = "zn";	element[30].mass = 65.39000;	element[30].name = "Zinc";
	element[31].symbol = "ga";	element[31].mass = 69.72300;	element[31].name = "Gallium";
	element[32].symbol = "ge";	element[32].mass = 72.64000;	element[32].name = "Germanium";
	element[33].symbol = "as";	element[33].mass = 74.92160;	element[33].name = "Arsenic";
	element[34].symbol = "se";	element[34].mass = 78.96000;	element[34].name = "Selenium";
	element[35].symbol = "br";	element[35].mass = 79.90400;	element[35].name = "Bromine";
	element[36].symbol = "kr";	element[36].mass = 83.80000;	element[36].name = "Krypton";
	element[37].symbol = "rb";	element[37].mass = 85.46780;	element[37].name = "Rubidium";
	element[38].symbol = "sr";	element[38].mass = 87.62000;	element[38].name = "Strontium";
	element[39].symbol = "y";	element[39].mass = 88.90590;	element[39].name = "Yttrium";
	element[40].symbol = "zr";	element[40].mass = 91.22400;	element[40].name = "Zirconium";
	element[41].symbol = "nb";	element[41].mass = 92.90640;	element[41].name = "Niobium";
	element[42].symbol = "mo";	element[42].mass = 95.94000;	element[42].name = "Molybdenum";
	element[43].symbol = "tc";	element[43].mass = 98.00000;	element[43].name = "Technetium";
	element[44].symbol = "ru";	element[44].mass = 101.07000;	element[44].name = "Ruthenium";
	element[45].symbol = "rh";	element[45].mass = 102.90550;	element[45].name = "Rhodium";
	element[46].symbol = "pd";	element[46].mass = 106.42000;	element[46].name = "Palladium";
	element[47].symbol = "ag";	element[47].mass = 107.86820;	element[47].name = "Silver";
	element[48].symbol = "cd";	element[48].mass = 112.41100;	element[48].name = "Cadmium";
	element[49].symbol = "in";	element[49].mass = 114.81800;	element[49].name = "Indium";
	element[50].symbol = "sn";	element[50].mass = 118.71000;	element[50].name = "Tin";
	element[51].symbol = "sb";	element[51].mass = 121.76000;	element[51].name = "Antimony";
	element[53].symbol = "i";	element[53].mass = 126.90450;	element[53].name = "Iodine";
	element[52].symbol = "te";	element[52].mass = 127.60000;	element[52].name = "Tellurium";
	element[54].symbol = "xe";	element[54].mass = 131.29300;	element[54].name = "Xenon";
	element[55].symbol = "cs";	element[55].mass = 132.90550;	element[55].name = "Cesium";
	element[56].symbol = "ba";	element[56].mass = 137.32700;	element[56].name = "Barium";
	element[57].symbol = "la";	element[57].mass = 138.90550;	element[57].name = "Lanthanum";
	element[58].symbol = "ce";	element[58].mass = 140.11600;	element[58].name = "Cerium";
	element[59].symbol = "pr";	element[59].mass = 140.90770;	element[59].name = "Praseodymium";
	element[60].symbol = "nd";	element[60].mass = 144.24000;	element[60].name = "Neodymium";
	element[61].symbol = "pm";	element[61].mass = 145.00000;	element[61].name = "Promethium";
	element[62].symbol = "sm";	element[62].mass = 150.36000;	element[62].name = "Samarium";
	element[63].symbol = "eu";	element[63].mass = 151.96400;	element[63].name = "Europium";
	element[64].symbol = "gd";	element[64].mass = 157.25000;	element[64].name = "Gadolinium";
	element[65].symbol = "tb";	element[65].mass = 158.92530;	element[65].name = "Terbium";
	element[66].symbol = "dy";	element[66].mass = 162.50000;	element[66].name = "Dysprosium";
	element[67].symbol = "ho";	element[67].mass = 164.93030;	element[67].name = "Holmium";
	element[68].symbol = "er";	element[68].mass = 167.25900;	element[68].name = "Erbium";
	element[69].symbol = "tm";	element[69].mass = 168.93420;	element[69].name = "Thulium";
	element[70].symbol = "yb";	element[70].mass = 173.04000;	element[70].name = "Ytterbium";
	element[71].symbol = "lu";	element[71].mass = 174.96700;	element[71].name = "Lutetium";
	element[72].symbol = "hf";	element[72].mass = 178.49000;	element[72].name = "Hafnium";
	element[73].symbol = "ta";	element[73].mass = 180.94790;	element[73].name = "Tantalum";
	element[74].symbol = "w";	element[74].mass = 183.84000;	element[74].name = "Tungsten";
	element[75].symbol = "re";	element[75].mass = 186.20700;	element[75].name = "Rhenium";
	element[76].symbol = "os";	element[76].mass = 190.23000;	element[76].name = "Osmium";
	element[77].symbol = "ir";	element[77].mass = 192.21700;	element[77].name = "Iridium";
	element[78].symbol = "pt";	element[78].mass = 195.07800;	element[78].name = "Platinum";
	element[79].symbol = "au";	element[79].mass = 196.96650;	element[79].name = "Gold";
	element[80].symbol = "hg";	element[80].mass = 200.59000;	element[80].name = "Mercury";
	element[81].symbol = "tl";	element[81].mass = 204.38330;	element[81].name = "Thallium";
	element[82].symbol = "pb";	element[82].mass = 207.20000;	element[82].name = "Lead";
	element[83].symbol = "bi";	element[83].mass = 208.98040;	element[83].name = "Bismuth";
	element[84].symbol = "po";	element[84].mass = 209.00000;	element[84].name = "Polonium";
	element[85].symbol = "at";	element[85].mass = 210.00000;	element[85].name = "Astatine";
	element[86].symbol = "rn";	element[86].mass = 222.00000;	element[86].name = "Radon";
	element[87].symbol = "fr";	element[87].mass = 223.00000;	element[87].name = "Francium";
	element[88].symbol = "ra";	element[88].mass = 226.00000;	element[88].name = "Radium";
	element[89].symbol = "ac";	element[89].mass = 227.00000;	element[89].name = "Actinium";
	element[91].symbol = "pa";	element[91].mass = 231.03590;	element[91].name = "Protactinium";
	element[90].symbol = "th";	element[90].mass = 232.03810;	element[90].name = "Thorium";
	element[93].symbol = "np";	element[93].mass = 237.00000;	element[93].name = "Neptunium";
	element[92].symbol = "u";	element[92].mass = 238.02890;	element[92].name = "Uranium";
	element[95].symbol = "am";	element[95].mass = 243.00000;	element[95].name = "Americium";
	element[94].symbol = "pu";	element[94].mass = 244.00000;	element[94].name = "Plutonium";
	element[96].symbol = "cm";	element[96].mass = 247.00000;	element[96].name = "Curium";
	element[97].symbol = "bk";	element[97].mass = 247.00000;	element[97].name = "Berkelium";
	element[98].symbol = "cf";	element[98].mass = 251.00000;	element[98].name = "Californium";
	element[99].symbol = "es";	element[99].mass = 252.00000;	element[99].name = "Einsteinium";
	element[100].symbol = "fm";	element[100].mass = 257.00000;	element[100].name = "Fermium";
	element[101].symbol = "md";	element[101].mass = 258.00000;	element[101].name = "Mendelevium";
	element[102].symbol = "no";	element[102].mass = 259.00000;	element[102].name = "Nobelium";
	element[103].symbol = "lr";	element[103].mass = 262.00000;	element[103].name = "Lawrencium";
	element[104].symbol = "rf";	element[104].mass = 261.00000;	element[104].name = "Rutherfordium";
	element[105].symbol = "db";	element[105].mass = 262.00000;	element[105].name = "Dubnium";
	element[106].symbol = "sg";	element[106].mass = 266.00000;	element[106].name = "Seaborgium";
	element[107].symbol = "bh";	element[107].mass = 264.00000;	element[107].name = "Bohrium";
	element[108].symbol = "hs";	element[108].mass = 277.00000;	element[108].name = "Hassium";
	element[109].symbol = "mt";	element[109].mass = 268.00000;	element[109].name = "Meitnerium";

	return element;
}

/* ---------------------------------------------------------------------- */
char* toLowerCase(char* str)
{	
	int differ = 'A'-'a';
	char ch;
	int ii = strlen(str);
	for (int i=0; i <ii;i++)
	{
		strncpy(&ch,str+i,1);
		if (ch>='A' && ch<='Z')
		{
			ch = ch-differ;
			memcpy(str+i,&ch,1);
		}
	}
	return str;
}

/* ---------------------------------------------------------------------- */
void copyright(ostream *outf)
{
  char null[1024];
	strcpy(null,"**********************************************************************************");
	*outf<<null<<endl;
	strcpy(null,"*                     Two-Phase Thermodynamics (2PT) Program                     *");
	*outf<<null<<endl;
	strcpy(null,"*                          Shiang-Tai Lin (stlin@ntw.edu.tw)                     *");
	*outf<<null<<endl;
	strcpy(null,"* Department of Chemical Engineering, National Tiawan University, Taipei, Taiwan *");
	*outf<<null<<endl;
	strcpy(null,"*                   Prabal K Maiti (maiti@physics.iisc.ernet.in)                 *");
	*outf<<null<<endl;
	strcpy(null,"*   Department of Physics, Indian Institute of Science, Bangalore, India 560012  *");
	*outf<<null<<endl;
	strcpy(null,"*                          Tod A Pascal (tpascal@wag.caltech.edu)                *");
	*outf<<null<<endl;
	strcpy(null,"*                                     and                                        *");
	*outf<<null<<endl;
	strcpy(null,"*                        William A Goddard III (wag@wag.caltech.edu)             *");
	*outf<<null<<endl;
	strcpy(null,"*         Materials and Process Simulation Center, Caltech, Pasadena, CA USA     *");
	*outf<<null<<endl;
	strcpy(null,"*                                Copyright (c) 2010                              *");
	*outf<<null<<endl;
	strcpy(null,"*                                All rights Reserved                             *");
	*outf<<null<<endl;
	strcpy(null,"*                             *************************                          *");
	*outf<<null<<endl;
	strcpy(null,"*                                      Cite:                                     *");
	*outf<<null<<endl;
	strcpy(null,"* S.T. Lin, M. Blanco and W.A. Goddard, J. Chem. Phys., 2003, 119, 11792-11805   *");
	*outf<<null<<endl;
	strcpy(null,"* S.T. Lin, P.K. Maiti and W.A. Goddard, J. Phys. Chem. B., 2010, 114, 8191-8198 *");
	*outf<<null<<endl;
	strcpy(null,"* T.A. Pascal, S.T. Lin and W.A. Goddard, PCCP, 2011, 13 (1), 169 - 181          *");
	*outf<<null<<endl;
	strcpy(null,"**********************************************************************************");
	*outf<<null<<endl;
}
