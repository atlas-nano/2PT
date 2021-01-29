/*********************************************************************************
*			 Two-Phase Thermodynamics (2PT) Program					 *
*				  Shiang-Tai Lin (stlin@ntw.edu.tw)					 *
* Department of Chemical Engineering, National Tiawan University, Taipei, Taiwan *
*		 Prabal K. Maiti (maiti@physics.iisc.ernet.in)			  *
*  Department of Physics, Indian Institute of Science, Bangalore, India, 560012  *
*		  Tod A Pascal (tpascal@wag.caltech.edu)				*
*					 and										*
*		William A Goddard III (wag@wag.caltech.edu)			 *
* Materials and Process Simulation Center, Caltech, Pasadena, CA USA	 *
*			Copyright (c) 2010							  *
*			All rights Reserved							 *
*		 *************************						  *
*		  Cite:									 *
* S.T. Lin, M. Blanco and W.A. Goddard, J. Chem. Phys., 2003, 119, 11792-11805   *
* S.T. Lin, P.K. Maiti and W.A. Goddard, J. Phys. Chem. B., 2010, 114, 8191-8198 *
* T.A. Pascal, S.T. Lin and W.A. Goddard, PCCP, 2011, 13(1), 169-181			 *
***********************************************************************************/

//read atom trajectories

#include <iostream>
#include <vector>
#include <sstream>
#include "trj_reader.h"
//#include <sys/mman.h>
//#include <sys/stat.h>

/* ---------------------------------------------------------------------- */
TRJ::TRJ() 
{

	itrjformat=0;
	otrjformat=0;
	ascfid=datasize2fid=0;
	location0e=datasize2fide=0;
	location0s=datasize2fids=0;
	has_atom_eng=0;
}

/* ---------------------------------------------------------------------- */
TRJ::~TRJ() 
{
	if(intrj.is_open()) intrj.close();
	if(ineng.is_open()) ineng.close();
	if(instr.is_open()) instr.is_open();
	if(inatomeng.is_open()) inatomeng.close();
}

/* ---------------------------------------------------------------------- */
int TRJ::init_contentDouble()
{
	
	if(header.nmovatm==0) return 1;

	return 0;
}

/* ---------------------------------------------------------------------- */
void TRJ::CleanTRJcontentDouble()
{
	int i,j;

	prp->time=0.0;
	prp->step=0;
	prp->T=avetem=dstep=firstt=finalt=0.0;
	prp->Ep=prp->Ebond=prp->Eangle=prp->Etorsion=prp->Einversion
		=prp->Evdw=prp->Eel=prp->Ehb=prp->Et=prp->Ek;
	//e=eb=et=ep=ei=enb=eel=ehb
	ec=eu=tint=tnb=ea=eba
		=eta=epa=eia=enba=eela=ehba
		=eca=eua=tinta=tnba   //=tote=totke
		=totea=tkea=0.0;
	for(i=0;i<12;i++) dum[i]=duma[i]=0.0;
	iconmp=imstep=iconfs=icstep=lvelwr=lfrcwr=0;

	prp->P=prp->V=pvtota=pvkina=pvpota=radgyra=0.0;
	pressur =vol =pvtot =pvkin =pvpot =radgyr =0.0;

	signose=zfrict=snose=snoseh=ssdot=sigdyn[0]=sigdyn[1]=0.0;
	zprfrict=qcanon=gamtmp=0.0;
	tcel=ucel=tcela=ucela=0.0;
	//for(i=0;i<6;i++) s2r[i]=s2ra[i]=s2rdot[i]=0.0;
	for(i=0;i<6;i++) s2ra[i]=s2rdot[i]=0.0;
	//for(i=0;i<3;i++) for(j=0;j<3;j++) cell->H[i][j]=0;

	natmcel=0;
	for(i=0;i<6;i++) strsa[i]=prp->strs[i]=0.0;
	extstrs=extstrsa=eabtota=eabvala=eabelha=eabnba=eabmisa=dltfaba=expprta=0.0;
	for(i=0;i<header.nmovatm;i++) //x[i]=y[i]=z[i]=velx[i]=vely[i]=velz[i]=0.0;
		for(j=0;j<3;j++)
		atom[i].pv[j]=atom[i].vel[j]=0.0;
}

/* ---------------------------------------------------------------------- */
void TRJ::ReadBinEng()
{
	int i;
	double dnull;
	 double crtn; //two carriage returns (each carriage return in fortran is 4 bytes)

	intrj.read((char *)&dnull,sizeof(dnull));
	prp->time=dnull; //current time (ps)
	intrj.read((char *)&(prp->step),sizeof(prp->step)); //number of steps
	intrj.read((char *)&(prp->T),sizeof(prp->T)); //instantaneous temperature
	intrj.read((char *)&avetem,sizeof(avetem)); //running averaged temperature
	intrj.read((char *)&dstep,sizeof(dstep)); //time step (ps)
	intrj.read((char *)&firstt,sizeof(firstt)); //initial temperature setting
	intrj.read((char *)&finalt,sizeof(finalt)); //final temperature setting
	intrj.read((char *)&(prp->Ep),sizeof(prp->Ep)); //potential energy
	intrj.read((char *)&(prp->Ebond),sizeof(prp->Ebond)); //bond
	intrj.read((char *)&(prp->Eangle),sizeof(prp->Eangle)); //angle
	intrj.read((char *)&(prp->Etorsion),sizeof(prp->Etorsion)); //torsion
	intrj.read((char *)&(prp->Einversion),sizeof(prp->Einversion)); //inversion
	intrj.read((char *)&(prp->Evdw),sizeof(prp->Evdw)); //vdw
	intrj.read((char *)&(prp->Eel),sizeof(prp->Eel)); //coulobm
	intrj.read((char *)&(prp->Ehb),sizeof(prp->Ehb)); //h-bond
	intrj.read((char *)&ec,sizeof(ec)); //restraint
	for(i=0;i<12;i++) intrj.read((char *)&dum[i],sizeof(dum[i])); //12 cross terms
	intrj.read((char *)&eu,sizeof(eu)); //user
	intrj.read((char *)&tint,sizeof(tint)); //valence 
	intrj.read((char *)&tnb,sizeof(tnb)); //nonbond
	intrj.read((char *)&ea,sizeof(ea));
	intrj.read((char *)&eba,sizeof(eba));
	intrj.read((char *)&eta,sizeof(eta));
	intrj.read((char *)&epa,sizeof(epa));
	intrj.read((char *)&eia,sizeof(eia));
	intrj.read((char *)&enba,sizeof(enba));
	intrj.read((char *)&eela,sizeof(eela));
	intrj.read((char *)&ehba,sizeof(ehba));
	intrj.read((char *)&eca,sizeof(eca));
	intrj.read((char *)&eua,sizeof(eua));
	intrj.read((char *)&tinta,sizeof(tinta));
	for(i=0;i<12;i++) intrj.read((char *)&duma[i],sizeof(duma[i]));
	intrj.read((char *)&tnba,sizeof(tnba));
	intrj.read((char *)&(prp->Et),sizeof(prp->Et));
	intrj.read((char *)&(prp->Ek),sizeof(prp->Ek));
	intrj.read((char *)&totea,sizeof(totea));
	intrj.read((char *)&tkea,sizeof(tkea));
	intrj.read((char *)&iconmp,sizeof(iconmp)); //flag
	intrj.read((char *)&imstep,sizeof(imstep)); //flag
	intrj.read((char *)&lvelwr,sizeof(lvelwr)); //write velocity
	intrj.read((char *)&lfrcwr,sizeof(lfrcwr)); //write force
	intrj.read((char *)&iconfs,sizeof(iconfs)); //flag
	intrj.read((char *)&icstep,sizeof(icstep)); //flag
	intrj.read((char *)&crtn,sizeof(crtn));

	if (header.version>=155) {
		intrj.read((char *)&pressur,sizeof(pressur));
		intrj.read((char *)&vol,sizeof(vol));
		intrj.read((char *)&pvtot,sizeof(pvtot));
		intrj.read((char *)&pvkin,sizeof(pvkin));
		intrj.read((char *)&pvpot,sizeof(pvpot));
		intrj.read((char *)&radgyr,sizeof(radgyr));
		intrj.read((char *)&(prp->P),sizeof(prp->P));
		intrj.read((char *)&(prp->V),sizeof(prp->V));
		intrj.read((char *)&pvtota,sizeof(pvtota));
		intrj.read((char *)&pvkina,sizeof(pvkina));
		intrj.read((char *)&pvpota,sizeof(pvpota));
		intrj.read((char *)&radgyra,sizeof(radgyra));
		intrj.read((char *)&crtn,sizeof(crtn));
	}

	if (header.lcanon) {
		if (header.version<300) {
			intrj.read((char *)&signose,sizeof(signose));
			intrj.read((char *)&zfrict,sizeof(zfrict));
			intrj.read((char *)&zprfrict,sizeof(zprfrict));
			intrj.read((char *)&crtn,sizeof(crtn));
		} else {
			if(header.lnose) {
				intrj.read((char *)&snose,sizeof(snose));
				intrj.read((char *)&snoseh,sizeof(snoseh));
				intrj.read((char *)&ssdot,sizeof(ssdot));
				intrj.read((char *)&qcanon,sizeof(qcanon));
				intrj.read((char *)&crtn,sizeof(crtn));
			} else {
				intrj.read((char *)&signose,sizeof(signose));
				intrj.read((char *)&zfrict,sizeof(zfrict));
				intrj.read((char *)&zprfrict,sizeof(zprfrict));
				intrj.read((char *)&qcanon,sizeof(qcanon));
				intrj.read((char *)&crtn,sizeof(crtn));
			}
		}
	}
						
	if(header.version>=220) {
		if(header.period) {
			intrj.read((char *)&tcel,sizeof(tcel));
			intrj.read((char *)&tcela,sizeof(tcela));
			//for(i=0;i<6;i++) intrj.read((char *)&s2r[i],sizeof(s2r[i]));
			intrj.read((char *)&(cell->H[0][0]),sizeof(cell->H[0][0]));
			intrj.read((char *)&(cell->H[1][1]),sizeof(cell->H[1][1]));
			intrj.read((char *)&(cell->H[2][2]),sizeof(cell->H[2][2]));
			intrj.read((char *)&(cell->H[2][1]),sizeof(cell->H[2][1]));
			intrj.read((char *)&(cell->H[2][0]),sizeof(cell->H[2][0]));
			intrj.read((char *)&(cell->H[1][0]),sizeof(cell->H[1][0]));
			for(i=0;i<6;i++) intrj.read((char *)&s2rdot[i],sizeof(s2rdot[i]));
			intrj.read((char *)&ucel,sizeof(ucel));
			intrj.read((char *)&ucela,sizeof(ucela));
			for(i=0;i<6;i++) intrj.read((char *)&s2ra[i],sizeof(s2ra[i]));
			intrj.read((char *)&crtn,sizeof(crtn));
		}
	} else if (header.defcel) {
		intrj.read((char *)&tcela,sizeof(tcela));
		//for(i=0;i<6;i++) intrj.read((char *)&s2r[i],sizeof(s2r[i]));
		intrj.read((char *)&(cell->H[0][0]),sizeof(cell->H[0][0]));
		intrj.read((char *)&(cell->H[1][1]),sizeof(cell->H[1][1]));
		intrj.read((char *)&(cell->H[2][2]),sizeof(cell->H[2][2]));
		intrj.read((char *)&(cell->H[2][1]),sizeof(cell->H[2][1]));
		intrj.read((char *)&(cell->H[2][0]),sizeof(cell->H[2][0]));
		intrj.read((char *)&(cell->H[1][0]),sizeof(cell->H[1][0]));
		intrj.read((char *)&crtn,sizeof(crtn));
	}

	if (header.period) {
		if ( (header.version==210) || (header.version>=300) ) {
			intrj.read((char *)&natmcel,sizeof(natmcel));
			for(i=0;i<6;i++) intrj.read((char *)&(prp->strs[i]),sizeof(prp->strs[i]));
			intrj.read((char *)&extstrs,sizeof(extstrs));
			for(i=0;i<6;i++) intrj.read((char *)&strsa[i],sizeof(strsa[i]));
			intrj.read((char *)&extstrsa,sizeof(extstrsa));
			intrj.read((char *)&crtn,sizeof(crtn));
		} else {
			intrj.read((char *)&natmcel,sizeof(natmcel));
			for(i=0;i<6;i++) intrj.read((char *)&(prp->strs[i]),sizeof(prp->strs[i]));
			intrj.read((char *)&crtn,sizeof(crtn));
		}
	}

	if ( header.version>=300 ) {
		if (header.period && header.lnpecan ) {
			intrj.read((char *)&sigdyn[0],sizeof(sigdyn[0]));
			intrj.read((char *)&sigdyn[1],sizeof(sigdyn[1]));
			intrj.read((char *)&qcanon,sizeof(qcanon));
			intrj.read((char *)&crtn,sizeof(crtn));
		}
		if (header.ltmpdamp) {
			intrj.read((char *)&gamtmp,sizeof(gamtmp));
			intrj.read((char *)&crtn,sizeof(crtn));
		}
	}

	if (header.prtthrm) {
		intrj.read((char *)&eabtota,sizeof(eabtota));
		intrj.read((char *)&eabvala,sizeof(eabvala));
		intrj.read((char *)&eabelha,sizeof(eabelha));
		intrj.read((char *)&eabnba,sizeof(eabnba));
		intrj.read((char *)&eabmisa,sizeof(eabmisa));
		intrj.read((char *)&dltfaba,sizeof(dltfaba));
		intrj.read((char *)&expprta,sizeof(expprta));
		intrj.read((char *)&crtn,sizeof(crtn));
	}

}

/* ---------------------------------------------------------------------- */
void TRJ::ReadBinFrame()
{
	int i;
	double crtn; //two carriage returns (each carriage return in fortran is 4 bytes)

	CleanTRJcontentDouble();
	ReadBinEng();

	for(i=0;i<header.nmovatm;i++) intrj.read((char *)&(atom[i].pv[0]),sizeof((atom[i].pv[0])));
	intrj.read((char *)&crtn,sizeof(crtn));
	for(i=0;i<header.nmovatm;i++) intrj.read((char *)&(atom[i].pv[1]),sizeof((atom[i].pv[1])));
	intrj.read((char *)&crtn,sizeof(crtn));
	for(i=0;i<header.nmovatm;i++) intrj.read((char *)&(atom[i].pv[2]),sizeof((atom[i].pv[2])));
	intrj.read((char *)&crtn,sizeof(crtn));

	//	 ----- velocities if needed -----
	if ( lvelwr )  {
		for(i=0;i<header.nmovatm;i++) intrj.read((char *)&(atom[i].vel[0]),sizeof(atom[i].vel[0]));
		intrj.read((char *)&crtn,sizeof(crtn));
		for(i=0;i<header.nmovatm;i++) intrj.read((char *)&(atom[i].vel[1]),sizeof((atom[i].vel[1])));
		intrj.read((char *)&crtn,sizeof(crtn));
		for(i=0;i<header.nmovatm;i++) intrj.read((char *)&(atom[i].vel[2]),sizeof((atom[i].vel[2])));
		intrj.read((char *)&crtn,sizeof(crtn));
	}
}

/* ---------------------------------------------------------------------- */
void TRJ::ReadAscFrame()
{
	int i;

	CleanTRJcontentDouble();

	intrj>>prp->time>>prp->step>>prp->T>>avetem>>dstep>>firstt
		>>finalt>>prp->Ep>>prp->Ebond>>prp->Eangle>>prp->Etorsion>>prp->Einversion
		>>prp->Evdw>>prp->Eel>>prp->Ehb>>ec>>eu>>tint
		>>tnb>>ea>>eba>>eta>>epa>>eia
		>>enba>>eela>>ehba>>eca>>eua>>tinta
		>>tnba>>prp->Et>>prp->Ek>>totea>>tkea;
	for(i=0;i<12;i++) intrj>>dum[i];
	for(i=0;i<12;i++) intrj>>duma[i];
	intrj>>iconmp>>imstep>>lvelwr>>lfrcwr>>iconfs>>icstep;

	if (header.version>=155) {
		intrj>>pressur>>vol>>pvtot
			>>pvkin>>pvpot>>radgyr;
		intrj>>prp->P>>prp->V>>pvtota
			>>pvkina>>pvpota>>radgyra;
	}

	if (header.lcanon) {
		if (header.version<300) {
			intrj>>signose>>zfrict>>zprfrict;
		} else {
			if(header.lnose)
				intrj>>snose>>snoseh>>ssdot>>qcanon;
			else
				intrj>>signose>>zfrict>>zprfrict>>qcanon;
		}
	}

	if(header.version>=220) {
		if(header.period) {
			intrj>>tcel>>tcela;
			//for(i=0;i<6;i++) intrj>>s2r[i];
			intrj>>cell->H[0][0]>>cell->H[1][1]>>cell->H[2][2]>>cell->H[2][1]>>cell->H[2][0]>>cell->H[1][0];
			for(i=0;i<6;i++) intrj>>s2rdot[i];
			intrj>>ucel>>ucela;
			for(i=0;i<6;i++) intrj>>s2ra[i];
		}
	} else if (header.defcel) {
		intrj>>tcela;
		//for(i=0;i<6;i++) intrj>>s2r[i];
		intrj>>cell->H[0][0]>>cell->H[1][1]>>cell->H[2][2]>>cell->H[2][1]>>cell->H[2][0]>>cell->H[1][0];
	}

	if (header.period) {
		if ( (header.version==210) || (header.version>=300) ) {
			intrj>>natmcel;
			for(i=0;i<6;i++) intrj>>prp->strs[i];
			intrj>>extstrs;
			for(i=0;i<6;i++) intrj>>strsa[i];
			intrj>>extstrsa;
		} else {
			intrj>>natmcel;
			for(i=0;i<6;i++) intrj>>prp->strs[i];
		}
	}

	if ( header.version>=300 ) {
		if (header.period && header.lnpecan ) intrj>>sigdyn[0]>>sigdyn[1]>>qcanon;
		if (header.ltmpdamp)  intrj>>gamtmp;
	}

	if (header.prtthrm) {
		intrj>>eabtota>>eabvala>>eabelha>>eabnba
			>>eabmisa>>dltfaba>>expprta;
	}

	for(i=0;i<header.nmovatm;i++) intrj>>atom[i].pv[0];
	for(i=0;i<header.nmovatm;i++) intrj>>atom[i].pv[1];
	for(i=0;i<header.nmovatm;i++) intrj>>atom[i].pv[2];

	//	 ----- velocities if needed -----

	if ( lvelwr ) {
		for(i=0;i<header.nmovatm;i++) intrj>>atom[i].vel[0];
		for(i=0;i<header.nmovatm;i++) intrj>>atom[i].vel[1];
		for(i=0;i<header.nmovatm;i++) intrj>>atom[i].vel[2];
	}
}

/* ---------------------------------------------------------------------- */
void TRJ::ReadUSRFrame()
{
	int i,j;

	CleanTRJcontentDouble();
 
	intrj.read((char *)&(prp->step),sizeof(prp->step)); //number of steps
	intrj.read((char *)&prp->time,sizeof(prp->time)); //time step (ps)
	intrj.read((char *)&(prp->T),sizeof(prp->T)); //instantaneous temperature

	for(i=0;i<header.nmovatm;i++) {
		for(j=0;j<3;j++){
		intrj.read((char *)&(atom[i].pv[j]),sizeof(atom[i].pv[j]));
		intrj.read((char *)&(atom[i].vel[j]),sizeof(atom[i].vel[j]));
		atom[i].vel[j]*=sqrt(418.4032);
		}
	}
}

/* ---------------------------------------------------------------------- */
void TRJ::ReadXYZFrame(ifstream *trj, int iflag)
{
	int i,j;
	char null[1024];

	trj->getline(null,1024);
	sscanf(null,"%*s %*s %d%*s %*s %*s %lf%*s",&prp->step,&prp->Et);
	prp->time=prp->step*header.timestep;
	trj->getline(null,1024);
	for(i=0;i<header.nmovatm;i++) {
		if(iflag)
			*trj>>null>>atom[i].pv[0]>>atom[i].pv[1]>>atom[i].pv[2];
		else 
			*trj>>null>>atom[i].vel[0]>>atom[i].vel[1]>>atom[i].vel[2];
		
	}
	prp->Et *= 627.509; //Hartrees to kcal/mol
}

/* ---------------------------------------------------------------------- */
void TRJ::ReadAMBERFrame(ifstream *trj, double buflen, int iframe, int iflag, int boxflag)
{
	int i,j,chunk;
	double box[6];
	char junk[8],null[1024];
	float crtn;

	box[0] = box[1] = box[2] = 0.0;
	box[3] = box[4] = box[5] = 90.0;

	istringstream snap,ss;

	chunk = 8 * sizeof(char);

	std::vector<char> buffer(buflen); //allocate a temporary buffer to slurp the entire timestep data
	trj->read(&buffer[0],buflen); //read into the buffer
	snap.rdbuf()->pubsetbuf(&buffer[0],buflen); //set the ss buffer to point to the start of our temporary buffer

	//get timestep
	dstep=header.timestep;
	prp->step = iframe;
	prp->time = iframe * header.timestep;

	i = 0;
	j = 0;
	ss.width(8);
	while (snap.good() && ! snap.eof()) {
		snap.getline(null,1024);
		ss.clear();
		ss.str(null);
		while (ss.good() && ! ss.eof()) {
			ss.read(junk,chunk);
			if(iflag)
				atom[i].pv[j % 3] = atof(junk);
			else
				atom[i].vel[j % 3] = atof(junk)*sqrt(418.4032);
			j++;
			if(j % 3 == 0) i++;
			if(j % 10 == 0) ss.read((char *)&crtn,sizeof(crtn));
			if(i==header.nmovatm) break;
		}
		if(i==header.nmovatm) break;
	}
	if(boxflag) {
		i = 0;
		while (ss.good() && ! ss.eof() && ss.read(junk,chunk)) {
			box[i] = atof(junk);
			j++; i++;
			if(i>5) break;
			if(j % 10 == 0) ss.read((char *)&crtn,sizeof(crtn));
		}
		if(snap.good() && ! snap.eof()) { //read next line
			snap.getline(null,1024);
			ss.clear();
			ss.str(null);
			while (ss.good() && ! ss.eof() && ss.read(junk,chunk)) {
				box[i] = atof(junk);
				i++;
				if(i>5) break;
			}
		}
		cell->la = box[0];
		cell->lb = box[1];
		cell->lc = box[2];
		cell->alpha = box[3];
		cell->beta = box[4];
		cell->gamma = box[5];
		cell->labc2H();
		cell->cal_volume();
		prp->V = cell->volume;
	}

}

/* ---------------------------------------------------------------------- */
void TRJ::ReadCPMDFrame()
{
	int i,j,skip;
	double tmp;  //dummy variable
	char null[1024];

	CleanTRJcontentDouble();
	//read xyz and vel
	dstep=header.timestep;
	for(skip=0;skip<header.nxyz;skip++) {
		for(i=0;i<header.nmovatm;i++) {
			intrj>>null;
			if(strcmp(null,"<<<<<<")==0) {
				intrj.getline(null,1024); intrj>>null;
			}
			prp->step=(int)atof(null);
			intrj>>atom[i].pv[0]>>atom[i].pv[1]>>atom[i].pv[2]>>atom[i].vel[0]>>atom[i].vel[1]>>atom[i].vel[2];
			//cout<<"atom "<<i;
			//for(j=0;j<3;j++) cout<<" "<<atom[i].pv[j];
			prp->step+=header.totaccustep;
			for(j=0;j<3;j++) {
				atom[i].pv[j]*=Bohr2A; //in A
				//atom[i].vel[j]*=2187691.254E-1; //in A/ps
				atom[i].vel[j]*=Bohr2A/(au2fs/1000); //in A/ps
			}
			prp->time=(prp->step-1.0)*(header.timestep);
		}
	}
}

/* ---------------------------------------------------------------------- */
void TRJ::ReadLMPFrame(int buflen, int iframe)
{
	int i,j,k;
	double xlo,xhi,ylo,yhi,zlo,zhi;
	double xy,xz,yz;
	char line[1024],null[1024];
	float val;
	//float data[header.lmp_data_len+1];
	istringstream snap;

	iframe=iframe;
	CleanTRJcontentDouble();
	std::vector<char> buffer(buflen); //allocate a temporary buffer to slurp the entire timestep data
	intrj.read(&buffer[0],buflen); //read into the buffer
	snap.rdbuf()->pubsetbuf(&buffer[0],buflen); //set the ss buffer to point to the start of our temporary buffer

	//get timestep
	dstep=header.timestep;
	snap>>null>>null>>null;
	i=(int)(atof(null));
	if(i<0) { cerr<<endl<<"Error: Could not determine the timestep! Got "<<null<<endl; exit(1); }
	prp->step = i;
	prp->time=i*dstep;
	snap.getline(line,1024);

	//get volume and set cell info
	while(! strstr(line,"ITEM: BOX BOUNDS")) { snap.getline(line,1024); };
	if (! header.lmp_cell_format) {
		snap>>xlo>>xhi>>ylo>>yhi>>zlo>>zhi;
		xy=xz=yz=0;
	} else {
		snap>>xlo>>xhi>>xy>>ylo>>yhi>>xz>>zlo>>zhi>>yz;
	}
	cell->H[0][0]=(xhi-xlo);
	cell->H[0][1]=xy;
	cell->H[0][2]=xz;
	cell->H[1][0]=0;
	cell->H[1][1]=(yhi-ylo);
	cell->H[1][2]=yz;
	cell->H[2][0]=cell->H[2][1]=0;
	cell->H[2][2]=(zhi-zlo);
	cell->H2others();
	cell->volume = cell->la*cell->lb*cell->lc;
	cell->cal_volume();
	prp->V = cell->volume;
	//get atom data
	while(! strstr(line,"ITEM: ATOMS")) { snap.getline(line,1024); } //navigate to atom line
		for(k=0;k<header.nmovatm;k++) {
		snap>>null;
		i = (int)atoi(null)-1;
		for(j=1;j<header.lmp_data_len;j++) {
			if(header.lmp_data_index[i][j]) snap>>(*header.lmp_data_index[i][j]);
			else snap>>val;
		}
	}
	for(i=0;i<header.nmovatm;i++) {
		for(j=0;j<3;j++)
		atom[i].vel[j] *= 1000;
		//todo:unwrap and unscale coordinates...
	}
}

/* ---------------------------------------------------------------------- */
void TRJ::ReadLMPThermo() {

	double tmp,tmp1;  //dummy variable
	char data[1024],junk[1024];

	tmp = tmp1 = 0.0;
	while(ineng>>data) {
		if(strstr(data,"Step")) break;
		else if(!strcmp(data,"TotEng")) { ineng>>junk>>prp->Et; }
		else if(!strcmp(data,"KinEng")) { ineng>>junk>>prp->Ek; }
		else if(!strcmp(data,"Temp")) { ineng>>junk>>prp->T; }
		else if(!strcmp(data,"PotEng")) { ineng>>junk>>prp->Ep; }
		else if(!strcmp(data,"E_bond")) { ineng>>junk>>prp->Ebond; }
		else if(!strcmp(data,"E_angle")) { ineng>>junk>>prp->Eangle; }
		else if(!strcmp(data,"E_dihed")) { ineng>>junk>>prp->Etorsion; }
		else if(!strcmp(data,"E_impro")) { ineng>>junk>>prp->Einversion; }
		else if(!strcmp(data,"E_vdwl")) { ineng>>junk>>prp->Evdw; }
		else if(!strcmp(data,"E_coul")) { ineng>>junk>>tmp; }
		else if(!strcmp(data,"E_long")) { ineng>>junk>>tmp1; }
		else if(!strcmp(data,"E_hbond")) { ineng>>junk>>prp->Ehb; }
		else if(!strcmp(data,"Press")) { ineng>>junk>>prp->P; }
		else if(!strcmp(data,"Volume")) { ineng>>junk>>prp->V; }
		else if(!strcmp(data,"PXX")) { ineng>>junk>>prp->strs[0]; }
		else if(!strcmp(data,"PYY")) { ineng>>junk>>prp->strs[1]; }
		else if(!strcmp(data,"PZZ")) { ineng>>junk>>prp->strs[2]; }
		else if(!strcmp(data,"PXY")) { ineng>>junk>>prp->strs[3]; }
		else if(!strcmp(data,"PXZ")) { ineng>>junk>>prp->strs[4]; }
		else if(!strcmp(data,"PYZ")) { ineng>>junk>>prp->strs[5]; }
	}
	prp->Eel = tmp+tmp1;
	prp->P *= (101.325/1e6);
}

/* ---------------------------------------------------------------------- */
void TRJ::ReadLMPAtomEng(int buflen) 
{
	int i,j,count,grpid,index;
	char line[1024],null[1024];
	double atomeng,avg,tot,snap_tot;
	double *grpEng = new double [ngrp];

	for(i=0;i<ngrp;i++) grpEng[i] = 0.0;
	istringstream snap;

	std::vector<char> buffer(buflen); //allocate a temporary buffer to slurp the entire timestep data
	inatomeng.read(&buffer[0],buflen); //read into the buffer
	snap.rdbuf()->pubsetbuf(&buffer[0],buflen); //set the ss buffer to point to the start of our temporary buffer

	count=0;
	while(snap.good() && ! snap.eof()) {
		while(! strstr(line,"ITEM: ATOMS") && snap.good() && ! snap.eof())
			snap.getline(line,1024); //navigate to atom line
		if(! snap.good() || snap.eof()) break;
		count++;
		snap.getline(line,1024);
	}

	for(i=0;i<ngrp;i++) {
		grp[i].eng.dataflag=1;
		grp[i].eng.nXY=count;
		grp[i].eng.init();
		grpEng[i] = 0.0;
	}

	avg = tot = snap_tot = 0.0;
	snap.clear();
	snap.rdbuf()->pubsetbuf(&buffer[0],buflen); //set the ss buffer to point to the start of our temporary buffer
	snap.getline(line,1024);
	index = 0;
	while(snap.good() && ! snap.eof()) {
		while(! strstr(line,"ITEM: ATOMS") && snap.good() && ! snap.eof())
		snap.getline(line,1024); //navigate to atom line
		if(! snap.good() || snap.eof()) break;
		for(j=0;j<header.nmovatm;j++) {
			snap>>null;
			i = (int)atoi(null)-1;
			grpid=atom[i].grp;
			snap>>null;
			atomeng = (double)atof(null);
			atom[i].eng_tmp = atomeng;
			grpEng[grpid] += atomeng;
			tot += atomeng;
			snap_tot += atomeng;
		}
		for (i=0;i<ngrp-1;i++)
			grp[i].eng.add_data(grpEng[i]);
		grp[ngrp-1].eng.add_data(snap_tot);
		for (i=0;i<ngrp;i++) grpEng[i] = 0.0; 
		snap_tot = 0.0;
		snap.getline(line,1024);
		index++;
	}

	//get averages
	avg = tot/count;
	cout<<" total: "<<(avg*4.184)<<" (kJ/mol)"<<endl;
	delete [] grpEng;
}

/* ---------------------------------------------------------------------- */
int TRJ::init_trj(trjItem trjfname)
{
	int i;

	switch(itrjformat) {
		case c2trj:
			i=init_c2bintrj(trjfname);
			break;
		case asctrj:
			i=init_c2asctrj(trjfname);
			break;
		case usrtrj:
			i=init_usrtrj(trjfname);
			break;
		case cpmdtrj:
			i=init_cpmdtrj(trjfname);
			break;
		case xyztrj:
			i=init_xyztrj(trjfname);
			break;
		case lmptrj:
			i=init_lmptrj(trjfname);
			break;
		case ambertrj:
			i=init_ambertrj(trjfname);
			break;
		case none:
		default:
			i=0;
			break;
	}
	return i;
}

/* ---------------------------------------------------------------------- */
int TRJ::init_c2bintrj(trjItem trjf)
{
	cout<<"Analyzing Binary Cerius2 trajectory file "<<trjf.vel<<endl;
	intrj.open(trjf.vel,ios::in);
	if(!intrj.is_open()) {
		cout<<" Error: Trj file "<<trjf.vel<<" cannot be opened."<<endl;
		return filenotopen;
	}

	int i = rd_header(atom);  //read header
	if(i) return i;

	location0=intrj.tellg();
	init_contentDouble();
	CleanTRJcontentDouble();
	ReadBinFrame(); 
	location1=intrj.tellg();
	datasize=location1-location0;

	//find the total frame number
	intrj.seekg(0, ios::end);
	location2=intrj.tellg();
	totframe=(location2-location0+4)/datasize;

	return 0;
}

/* ---------------------------------------------------------------------- */
int TRJ::init_c2asctrj(trjItem trjf)
{
	cout<<"Analyzing ASCII Cerius2 trajectory file "<<trjf.vel<<endl;
	intrj.open(trjf.vel,ios::in);
	if(!intrj.is_open()) {
		cout<<" Error: Trj file "<<trjf.vel<<" cannot be opened."<<endl;
		return filenotopen;
	}

	int i = rd_header(atom);  //read header
	if(i) return i;

	location0=intrj.tellg();
	init_contentDouble();
	CleanTRJcontentDouble();
	ReadAscFrame();
	
	location1=intrj.tellg();
	datasize=location1-location0;
	ascfid=1; datasize2fid=datasize;

	//find the total frame number
	totframe=1;
	while(!intrj.eof()) {
		ReadAscFrame();
		totframe++;
	}
	totframe--;

	return 0;
}

/* ---------------------------------------------------------------------- */
int TRJ::init_usrtrj(trjItem trjf)
{

	cout<<"Analyzing user trajectory file "<<trjf.vel<<endl;
	intrj.open(trjf.vel,ios::in);
	if(!intrj.is_open()) {
		cout<<" Error: Trj file "<<trjf.vel<<" cannot be opened."<<endl;
		return filenotopen;
	}

	int i = rd_header(atom);  //read header
	if(i) return i;

	location0=intrj.tellg();
	init_contentDouble();
	CleanTRJcontentDouble();
	ReadUSRFrame();
	location1=intrj.tellg();
	datasize=location1-location0;

	//find the total frame number
	intrj.seekg(0, ios::end);
	location2=intrj.tellg();
	totframe=(location2-location0+4)/datasize;

	return 0;
}

/* ---------------------------------------------------------------------- */
int TRJ::init_cpmdtrj(trjItem trjf)
{
	cout<<"Analyzing CPMD trajectory file "<<trjf.vel<<endl;
	intrj.open(trjf.vel,ios::in);
	if(!intrj.is_open()) {
		cout<<" Error: Trj file "<<trjf.vel<<" cannot be opened."<<endl;
		return filenotopen;
	}

	int i  = rd_header(atom);  //read header
	//if(i) return i;

	location0=intrj.tellg();
	location0e=ineng.tellg();
	location0s=instr.tellg();
	init_contentDouble();
	CleanTRJcontentDouble();
	ReadCPMDFrame();
	location1=intrj.tellg();
	datasize=location1-location0;
	ascfid=1; datasize2fid=datasize;
	location1=ineng.tellg();
	datasize2fide=location1-location0e;
	location1=instr.tellg();
	datasize2fids=location1-location0s;
	//find the total frame number
	totframe=1;
	int tmp=0;
	while(intrj.eof()+instr.eof()+ineng.eof()==0) {
		ReadCPMDFrame();
		if(prp->step > tmp) tmp=prp->step;
		totframe++;
	}
	totframe--;
	header.totaccustep=tmp;

	return 0;
}

/* ---------------------------------------------------------------------- */
int TRJ::init_xyztrj(trjItem trjf)
{
	cout<<"Analyzing XYZ trajectory "<<endl<<"  Velfile: "<<trjf.vel<<"...";
	fflush(stdout);
	intrj.open(trjf.vel,ios::in);
	if(!intrj.is_open()) {
		cout<<" Error: Vel file "<<trjf.vel<<" cannot be opened."<<endl;
		return filenotopen;
	}
	if(trjf.molopt) {
		cout<<" Coordfile: "<<trjf.coord<<"...";
		incoord.open(trjf.coord,ios::in);
		if(!incoord.is_open()) {
			cout<<" Error: Coord file "<<trjf.coord<<" cannot be opened."<<endl;
			return filenotopen;
		}
	}

	int i = rd_header(atom);  //read header
	if(i) return i;

	location0=intrj.tellg();
	init_contentDouble();
	CleanTRJcontentDouble();
	intrj.clear();
	intrj.seekg(header.byte_offset[0].xyz, ios::beg);
	ReadXYZFrame(&intrj,0);
	totframe=header.lmptmp;

	return 0;
}


/* ---------------------------------------------------------------------- */
int TRJ::init_ambertrj(trjItem trjf)
{
	cout<<"Analyzing AMBER trajectory "<<endl<<"  Velfile: "<<trjf.vel<<"...";
	fflush(stdout);
	intrj.open(trjf.vel,ios::in);
	if(!intrj.is_open()) {
		cout<<" Error: Vel file "<<trjf.vel<<" cannot be opened."<<endl;
		return filenotopen;
	}
	if(trjf.molopt) {
		cout<<" Coordfile: "<<trjf.coord<<"...";
		incoord.open(trjf.coord,ios::in);
		if(!incoord.is_open()) {
			cout<<" Error: Coord file "<<trjf.coord<<" cannot be opened."<<endl;
			return filenotopen;
		}
	}

	int i = rd_header(atom);  //read header
	if(i) return i;
	location0=intrj.tellg();
	init_contentDouble();
	CleanTRJcontentDouble();
	intrj.clear();
	intrj.seekg(header.byte_offset[0].xyz, ios::beg);
	double buflen = header.byte_offset[1].xyz - header.byte_offset[0].xyz;
	ReadAMBERFrame(&intrj,buflen,0,0,header.byte_offset[0].vamberbox);
	totframe=header.lmptmp;

	return 0;
}

/* ---------------------------------------------------------------------- */
int TRJ::init_lmptrj(trjItem trjf)
{
	cout<<"Analyzing trajectory files : "<<endl;

	cout<<" TRAJECTORY : "<<trjf.vel<<endl;
	intrj.open(trjf.vel,ios::in);
	if(!intrj.is_open()) {
		cout<<" Error: Trj file "<<trjf.vel<<" cannot be opened."<<endl;
		return filenotopen;
	}

	ineng.open(trjf.thermo,ios::in);
	if(ineng.is_open())
		cout<<" ENERGIES   : "<<trjf.thermo<<endl;

	int i = rd_header(atom);  //read header
	if(i) return i;
	CleanTRJcontentDouble();

	inatomeng.open(trjf.atom_eng,ios::in);
	if(inatomeng.is_open()) {
		cout<<" ATOM ENERGY : "<<trjf.atom_eng<<endl;
		has_atom_eng = 1;
	}

	totframe=header.lmptmp-1;

	return 0;
}

/* ---------------------------------------------------------------------- */
int TRJ::rd_header(ATOM * atom)
{
	int i = 0;

	switch(itrjformat) {
		case c2trj:
			i = header.ReadBinHeader(&intrj);
			break;
		case asctrj:
			i = header.ReadAscHeader(&intrj);
			break;
		case lmptrj:
			i = header.ReadLMPHeader(&intrj,&ineng,atom);
			break;
		case cpmdtrj:
			header.ReadCPMDHeader(&intrj,&ineng,&instr);
			i = 1;
			break;
		case usrtrj:
			i = header.ReadUSRHeader(&intrj);
			break;
		case xyztrj:
			header.CleanTRJheader(0);
			i = header.ReadXYZHeader(&intrj,0);
			if(incoord.is_open() && ! i) i = header.ReadXYZHeader(&incoord,1);
			cout<<"Done"<<endl;
			break;
		case ambertrj:
			header.CleanTRJheader(0);
			i = header.ReadAMBERHeader(&intrj,0);
			if(incoord.is_open() && ! i) i = header.ReadAMBERHeader(&incoord,1);
			cout<<"Done"<<endl;
			break;
		case none:
		default:
			break;
	}
	return i;
}

/* ---------------------------------------------------------------------- */
void TRJ::rd_atom_eng()
{
	int buflen;

	inatomeng.clear();
	cout<<"  Reading atom energy";
	switch(itrjformat) {
		case lmptrj:
			inatomeng.seekg(0,std::ios::end); //go to the eng
			buflen = inatomeng.tellg();
			inatomeng.seekg(0,std::ios::beg); //go to the eng
			ReadLMPAtomEng(buflen);
			break;
		case none:
		default:
			break;
	}
}

/* ---------------------------------------------------------------------- */
void TRJ::rd_frame(int iframe)
{
	int bufsize,i;

	fcount=(iframe-1);
	intrj.clear();
	switch(itrjformat) {
		case c2trj:
			intrj.seekg(location0+datasize*fcount, ios::beg);
			ReadBinFrame();
			break;
		case usrtrj:
			intrj.seekg(location0+datasize*fcount, ios::beg);
			ReadUSRFrame();
			break;
		case cpmdtrj:
			ineng.clear();
			instr.clear();
			if(iframe>(int)ascfid) {
				intrj.seekg(location0 +datasize2fid , ios::beg);
				ineng.seekg(location0e+datasize2fide, ios::beg);
				instr.seekg(location0s+datasize2fids, ios::beg);
				for(i=ascfid;i<iframe;i++) ReadCPMDFrame();
				ascfid =iframe;
				datasize2fid = intrj.tellg();
				datasize2fid-= location0;
				datasize2fide = ineng.tellg();
				datasize2fide-= location0e;
				datasize2fids = instr.tellg();
				datasize2fids-= location0s;
			} else {
				intrj.seekg(location0 , ios::beg);
				ineng.seekg(location0e, ios::beg);
				instr.seekg(location0s, ios::beg);
				for(i=0;i<iframe;i++) ReadCPMDFrame();
				ascfid =iframe;
				datasize2fid = intrj.tellg();
				datasize2fid-= location0;
				datasize2fide = ineng.tellg();
				datasize2fide-= location0e;
				datasize2fids = instr.tellg();
				datasize2fids-= location0s;
			}
			break;
		case xyztrj:
			CleanTRJcontentDouble();
			intrj.clear();
			intrj.seekg(header.byte_offset[iframe-1].xyz, ios::beg);
			ReadXYZFrame(&intrj,0);
			if(incoord.is_open()) {
				incoord.clear();
				incoord.seekg(header.byte_offset[iframe-1].coord, ios::beg);
				ReadXYZFrame(&incoord,1);
			}
			break;
		case ambertrj:
			CleanTRJcontentDouble();
			intrj.clear();
			bufsize=header.byte_offset[iframe].xyz - header.byte_offset[iframe-1].xyz;
			intrj.seekg(header.byte_offset[iframe-1].xyz, ios::beg);
			ReadAMBERFrame(&intrj,bufsize,(iframe-1),0,header.byte_offset[iframe-1].vamberbox);
			if(incoord.is_open()) {
				incoord.clear();
				bufsize=header.byte_offset[iframe].coord - header.byte_offset[iframe-1].coord;
				incoord.seekg(header.byte_offset[iframe-1].coord, ios::beg);
				ReadAMBERFrame(&incoord,bufsize,(iframe-1),1,header.byte_offset[iframe-1].camberbox);
			}
			break;
		case asctrj:
			if(iframe>(signed)ascfid) {
				intrj.seekg(location0+datasize2fid, ios::beg);
				for(i=ascfid;i<iframe;i++) ReadAscFrame();
				ascfid =iframe;
				datasize2fid  = intrj.tellg();
				datasize2fid -= location0;
			} else {
				intrj.seekg(location0, ios::beg);
				for(i=0;i<iframe;i++) ReadAscFrame();
				ascfid =iframe;
				datasize2fid = intrj.tellg();
				datasize2fid-= location0;
			}
			break;
		case lmptrj:
			instr.clear();
			bufsize=header.byte_offset[iframe].xyz - header.byte_offset[iframe-1].xyz;
			intrj.seekg(header.byte_offset[iframe-1].xyz, ios::beg);
			ReadLMPFrame(bufsize,iframe);
			if(ineng.is_open()) {
				ineng.clear();
				ineng.seekg(header.byte_offset[iframe-1].eng, ios::beg);
				ReadLMPThermo();
			}
			break;
		case none:
			default:
			break;
	}

	cell->H2others();
	prp->rho=prp->mass/(prp->V*1E-24*Na); //density in g/cc

}

