/********************************************************************************
*								 Two-Phase Thermodynamics (2PT) Program													*
*										Shiang-Tai Lin (stlin@ntw.edu.tw)														*
* Department of Chemical Engineering, National Tiawan University, Taipei, Taiwan*
*							 Prabal K. Maiti (maiti@physics.iisc.ernet.in)										*
*	Department of Physics, Indian Institute of Science, Bangalore, India, 560012	*
*									Tod A Pascal (tpascal@wag.caltech.edu)												*
*																 and																						*
*								William A Goddard III (wag@wag.caltech.edu)											*
*			 Materials and Process Simulation Center, Caltech, Pasadena, CA USA	 			*
*													Copyright (c) 2010																		*
*													All rights Reserved																		*
*										 *************************																	*
*																Cite:									 													*
* S.T. Lin, M. Blanco and W.A. Goddard, J. Chem. Phys., 2003, 119, 11792-11805	*
* S.T. Lin, P.K. Maiti and W.A. Goddard, J. Phys. Chem. B., 2010, 114, 8191-8198*
* T.A. Pascal, S.T. Lin and W.A. Goddard, PCCP, 2011, 13(1), 169-181			 			*
*********************************************************************************/

#include <sstream>
#include <iostream>
#include "model.h"
#include "utility.h"
#include "statistics.h"
#include <algorithm>

//atom/molecular info
/* ---------------------------------------------------------------------- */
MODEL::MODEL () {
	nmol=nbond=nangle=ntor=ninv=nimp=ngrp=0;
	atom=NULL;
	bond=NULL;
	mol=NULL;
	grp=NULL;
	periodic=0;
	xyzfreq=velfreq=engfreq=strfreq=0;
	timestep=0;
}

/* ---------------------------------------------------------------------- */
MODEL::~MODEL () {
	if(atom!=NULL)	delete [] atom;
	if(bond!=NULL)	delete [] bond;
	if(mol!=NULL)	 delete [] mol;
	if(grp!=NULL)	 delete [] grp;
}

/* ---------------------------------------------------------------------- */
int MODEL::rd_strt()
{
	int errcode = 0;

	strcpy(name,ctl.in_strt);	//set the name of the model to the input filename
	switch(ctl.in_strt_flag) {
		case strtlmp:
			errcode = rd_lmpdata(name);
			break;
		case strtbgf:
			errcode = rd_bgf(name);
			if(errcode) break;
			cnt2valence();//set up bonds
			element2prp(1); //assign atom mass,color,radius,crad based on element type
			break;
		case strtamber:
			errcode = rd_amberprmtop(name);
			break;
	}
	if(errcode) {cout<<"error while reading structure"<<endl; return errcode; }
	bond2mol(); //find molecules
	cal_mass(); //calculate total mass of the model
	ck_range(); //check atom range for nonperiodic systems
	cout<<natom<<" atoms "<<nbond<<" bonds "<<nmol<<" molecules...Done"<<endl;
	return 0;
}

/* ---------------------------------------------------------------------- */
int MODEL::rd_lmpdata(char *indata)
{
	int i,j,natmt;
	double * atmtmass;
	double xhi,xlo,yhi,ylo,zhi,zlo;
	double xy,xz,yz;

	char null[1024],junk[1024];

	logical hasmass = 0;

	istringstream ss;

	xy=xz=yz=0; //orthorhombic cell

	cout<<"Reading lammps data file "<<indata<<"...";
	fflush(stdout);

	ifstream inf(indata,fstream::in);

	if(!inf.is_open()) {
		cout<<" Error: Lammps data file "<<indata<<" cannot be opened."<<endl;
		return filenotopen;
	}
	

	inf.getline(null,1024);
	natmt=0;
	while (! inf.eof()) {
		inf.getline(null,1024);
		if (strspn(null," \t\n\r") == strlen(null)) continue;

		else if (strstr(null,"atoms")) { sscanf(null,"%d",&natom); init_atom(); }
		else if (strstr(null,"bonds")) { sscanf(null,"%d",&nbond); init_bond(); }

		else if (strstr(null,"atom types")) { sscanf(null,"%d",&natmt); }

		else if (strstr(null,"xlo xhi")) sscanf(null,"%lg %lg",&xlo,&xhi);
		else if (strstr(null,"ylo yhi")) sscanf(null,"%lg %lg",&ylo,&yhi);
		else if (strstr(null,"zlo zhi")) sscanf(null,"%lg %lg",&zlo,&zhi);
		else if (strstr(null,"xy xz yz")) sscanf(null,"%lg %lg %lg",&xy,&xz,&yz);
	}

	if(!natom || !natmt) {
		cerr<<"LAMMPS data file is invalid"<<endl;
		return filenotopen;
	}

	atmtmass = new double [natmt];
	periodic=3;
	cell.H[0][0]=(xhi-xlo);
	cell.H[0][1]=xy;
	cell.H[0][2]=xz;
	cell.H[1][0]=0;
	cell.H[1][1]=(yhi-ylo);
	cell.H[1][2]=yz;
	cell.H[2][0]=cell.H[2][1]=0;
	cell.H[2][2]=(zhi-zlo);
	cell.H2others();

	inf.clear();
	inf.seekg(0,ios::beg);
	while(!inf.eof()) {
		inf.getline(null,1024);
		if(strstr(null,"Masses")) {
			inf.getline(null,1024);
			for(i=0;i<natmt;i++)	{
				inf.getline(null,1024);
				ss.clear();
				ss.str(null);
				ss>>junk>>atmtmass[i];
				hasmass = 1;
			}
		}
	}
	if(! hasmass) {
		cout<<"ERROR: Masses not found in LAMMPS datafile"<<endl;
		return filenotopen;
	}

	inf.clear();
	inf.seekg(0,ios::beg);
	while(!inf.eof()) {
		inf.getline(null,1024);
		if(strstr(null,"Atoms")) {
			inf.getline(null,1024);
			for(j=0;j<natom;j++) {
				inf.getline(null,1024);
				ss.clear();
				ss.str(null);
				ss>>i;
	 			i--;
				atom[i].id = i;
				ss>>atom[i].mol>>atom[i].ffid>>atom[i].chg>>atom[i].pv[0]>>atom[i].pv[1]>>atom[i].pv[2];
				atom[i].ffid--;
				atom[i].mol=-1;
				atom[i].mass=atmtmass[atom[i].ffid];
			}

		} else if(strstr(null,"Bonds")&&nbond) {
			inf.getline(null,1024);
			for(i=0;i<nbond;i++) {
				inf>>bond[i].id>>bond[i].ffid>>bond[i].atm[0]>>bond[i].atm[1];
				bond[i].id--;
				bond[i].ffid--;
				bond[i].atm[0]--; bond[i].atm[1]--;
				bond[i].atom[0]=&atom[bond[i].atm[0]];
				bond[i].atom[1]=&atom[bond[i].atm[1]];
			}

		} else if(strstr(null,"Velocities")) {
			for(i=0;i<natom;i++) {
				inf>>null>>atom[i].vel[0]>>atom[i].vel[1]>>atom[i].vel[2];
				atom[i].vel[0] *= 1000.0;
				atom[i].vel[1] *= 1000.0;
				atom[i].vel[2] *= 1000.0;
			}
		}
	}
	delete [] atmtmass;
	return 0;
}

/* ---------------------------------------------------------------------- */
int MODEL::rd_bgf(char *indata)
{
	char null[1024];
	int i,j;	
	char *delims = {(char *)" ,\t"};
	char *charptr;
	int Symbsize=2;

	cout<<"Reading BGF file "<<indata<<"...";
	ifstream inf(indata,fstream::in);
	if(!inf.is_open()) {
		cout<<" Error: BGF file "<<indata<<" cannot be opened."<<endl;
		return filenotopen;
	}
	
	while(!inf.eof()) {
		inf.getline(null,1024);
		if (strncmp(null,"CRYSTX",6) == 0) break;
		else if ((strncmp(null,"HETATM",6) == 0) || (strncmp(null,"ATOM",4) == 0) ) break;
	}

	if(strncmp(null,"CRYSTX",6) == 0) {
		periodic = 3;
		charptr = strtok(null,delims);
		charptr = strtok(NULL, delims);
		cell.la=atof(charptr);	
		charptr = strtok(NULL, delims);
		cell.lb=atof(charptr);	
		charptr = strtok(NULL, delims);
		cell.lc=atof(charptr);	
		charptr = strtok(NULL, delims);
		cell.alpha=atof(charptr);	
		charptr = strtok(NULL, delims);
		cell.beta=atof(charptr);	
		charptr = strtok(NULL, delims);
		cell.gamma=atof(charptr);	
		while(!inf.eof()) {
			inf.getline(null,1024);
			if ((strncmp(null,"HETATM",6) == 0) || (strncmp(null,"ATOM",4) == 0) ) break;
		}
	} else {
		periodic = 0;
		cell.la=cell.lb=cell.lc=cell.alpha=cell.beta=cell.gamma=0;
	}
	cell.labc2H();
	cell.H2others();

	if((strncmp(null,"HETATM",6) == 0) || (strncmp(null,"ATOM",4) == 0)) {
		natom=1;
		while(!inf.eof()) {
			inf.getline(null,1024);
			if ((strncmp(null,"HETATM",6) == 0) || (strncmp(null,"ATOM",4) == 0) ) natom++;
			else if(strncmp(null,"CONECT",6) == 0) break;
		}
	} else {
		natom=0;
		return 1;
	}
	init_atom();
	inf.clear();
	inf.seekg(0,ios::beg);
	inf>>null;
	while((strncmp(null,"HETATM",6) != 0) && (strncmp(null,"ATOM",4) != 0)) inf>>null;
	if(strncmp(null,"ATOM",4)==0) { inf.getline(null,1024); inf>>null;}
	//get atomic info
	for(i=0;i<natom;i++) {
		atom[i].id=i;
		inf>>null>>atom[i].name>>null>>null>>null
	 		>>atom[i].pv[0]>>atom[i].pv[1]>>atom[i].pv[2]
	 		>>atom[i].fftype>>null>>null>>atom[i].chg>>null;
	//remove number in the atom name
	for(j=0;j<Symbsize;j++)	
		if(atom[i].name[j]>=48 && atom[i].name[j]<=57) 
		for(;j<Symbsize;j++) atom[i].name[j]='\0';
	}
	//get connectivity
	while(strncmp(null,"CONECT",6) != 0 && !inf.eof()) inf.getline(null,1024);
	i=0;
	while(!inf.eof()&&strncmp(null,"END",3)!=0) {
		if(strncmp(null,"CONECT",6)==0) {
			charptr = strtok(null,delims);
			charptr = strtok(NULL, delims);
			i=atoi(charptr)-1;
			atom[i].ncnt=0;
			charptr = strtok(NULL, delims);
			while(charptr!=NULL) {
				atom[i].bod[ atom[i].ncnt ] =1; //default bond order is 1
				atom[i].cnt[ atom[i].ncnt++ ]=atoi(charptr)-1;
				charptr = strtok(NULL, delims);
			}
		} else if(strncmp(null,"ORDER",5)==0) {
			charptr = strtok(null,delims);
			charptr = strtok(NULL, delims);
			if(i!=atoi(charptr)-1) {
				cout<<"Error: mismatch in rd_bgf cnect"<<endl;
				return 1;
			}
			for(j=0;j<atom[i].ncnt;j++) {
				charptr = strtok(NULL, delims);
				atom[i].bod[j]=atoi(charptr);
			}
		}
		inf.getline(null,1024);
	}
	return 0;
}

/* ---------------------------------------------------------------------- */
int MODEL::rd_amberprmtop(char *indata)
{
	char null[1024],junk[1024];
	int i;

	cout<<"Reading AMBER prmtop file "<<indata<<"...";
	ifstream inf(indata,fstream::in);
	if(!inf.is_open()) {
		cout<<" Error: AMBER prmtop file "<<indata<<" cannot be opened."<<endl;
		return filenotopen;
	}

	natom = 0;
	while(!inf.eof()) {
		inf.getline(null,1024);
		if (strstr(null,"%FLAG POINTERS")) {
			inf.getline(null,1024); //skip
			inf>>natom; //get number of atoms
			init_atom();
			inf>>junk; //skip
			inf>>nbond; //nbonds
			init_bond();
		} else if (natom && strstr(null,"%FLAG ATOM_NAME")) {
			inf.getline(null,1024); //skip
			for(i=0;i<natom;i++) {
				atom[i].id=i;
				atom[i].mol=-1;
				inf>>atom[i].name;
				if(! inf.good() || inf.eof()) return 1;
			}
		} else if (natom && strstr(null,"%FLAG CHARGE")) {
			inf.getline(null,1024); //skip
			for(i=0;i<natom;i++) { inf>>atom[i].chg; atom[i].chg /= 18.2223; if(! inf.good() || inf.eof()) return 1; }
		} else if (natom && strstr(null,"%FLAG MASS")) {
			inf.getline(null,1024); //skip
			for(i=0;i<natom;i++) { inf>>atom[i].mass; if(! inf.good() || inf.eof()) return 1; }
		} else if (natom && strstr(null,"%FLAG ATOM_TYPE_INDEX")) {
			inf.getline(null,1024); //skip
			for(i=0;i<natom;i++) { inf>>atom[i].ffid; if(! inf.good() || inf.eof()) return 1; }
		} else if (natom && strstr(null,"%FLAG AMBER_ATOM_TYPE")) {
			inf.getline(null,1024); //skip
			for(i=0;i<natom;i++) { inf>>atom[i].fftype; if(! inf.good() || inf.eof()) return 1; }
		} else if (natom && strstr(null,"%FLAG RADII")) {
			inf.getline(null,1024); //skip
			for(i=0;i<natom;i++) { inf>>atom[i].radius; if(! inf.good() || inf.eof()) return 1; }
		} else if (natom && strstr(null,"%FLAG BONDS_INC_HYDROGEN")) {
			inf.getline(null,1024); //skip
			for(i=0;i<nbond;i++) {
				bond[i].id = i;
				inf>>bond[i].atm[0]>>bond[i].atm[1]>>bond[i].ffid;
				bond[i].atm[0] /= 3; bond[i].atm[1] /= 3;
				bond[i].atom[0]=&atom[bond[i].atm[0]];
				bond[i].atom[1]=&atom[bond[i].atm[1]];
				if(! inf.good() || inf.eof()) return 1;
			}
		}
	}

	return 0;
}

/* ---------------------------------------------------------------------- */
int MODEL::bond2mol()
{
	int i,j,k,ck;
	int tota,fstatm;
	int *tmatom = new int [natom]; 

	nmol=0;
	fstatm=0;
	for(i=0;i<natom;i++) atom[i].mol=-1;
	atom[fstatm].mol=nmol;
	tmatom[nmol]=1;
	tota=0;
	while(tmatom[nmol]) {
		 MOLECULE tmol; tmol.natom=natom; tmol.init_atom();
		 tmol.atm[0]=fstatm;
		 for(i=0;i<tmatom[nmol];i++) {
			for(j=0;j<nbond;j++) {
				ck=-1;
				if(bond[j].atm[0]==tmol.atm[i]) ck=bond[j].atm[1];
			 	else if(bond[j].atm[1]==tmol.atm[i]) ck=bond[j].atm[0];
			 	if(ck>0) {
					//check if bonded atom in the list
					for(k=0;k<tmatom[nmol];k++) if(tmol.atm[k]==ck) break;
					if(k==tmatom[nmol]) { //add new atom to the list
						atom[ck].mol=nmol;
						tmol.atm[tmatom[nmol]++]=ck;
				 	}
			 	} 
			}
		 }
		 tota+=tmatom[nmol];
		 nmol++;
		 for(i=0;i<natom;i++) if(atom[i].mol<0) break;
		 if(i<natom) {
		 	tmatom[nmol]=1;
		 	fstatm=i;
		 	atom[fstatm].mol=nmol;
		 } else tmatom[nmol]=0;
	}
	if(tota!=natom) {
		cout<<" Error: bond2mol failed natom="<<natom<<" tota="<<tota<<endl;
		delete [] tmatom;
		return 1;
	}
	init_mol();
	for(i=0;i<nmol;i++) {
		 mol[i].natom=tmatom[i];
		 mol[i].init_atom();
		 mol[i].natom=0;
		 mol[i].mass=0;
	}
	for(i=0;i<natom;i++) {
		 j=atom[i].mol;
		 mol[j].atm[ mol[j].natom ]=i;
		 mol[j].atom[ mol[j].natom++ ]= & atom[i];
		 mol[j].mass+=atom[i].mass;
	}	

	if(DEBUG==2) {
		for(i=0;i<nmol;i++) {
			cout<<"mol "<<i;
			for(j=0;j<mol[i].natom;j++) cout<<" "<<mol[i].atm[j];
			cout<<endl;
		}
	}
	delete [] tmatom;
	return 0;
}

/* ---------------------------------------------------------------------- */
void MODEL::cnt2valence()
{	 //bgf2lmp

	int i,j,k,id1;

	//find bond
	nbond=0;
	for(i=0;i<natom;i++) {
		for(j=0;j<atom[i].ncnt;j++) {
			if(atom[i].id<atom[atom[i].cnt[j]].id) nbond++; 
		}
	}
	init_bond();
	k=0;
	for(i=0;i<natom;i++) {
		for(j=0;j<atom[i].ncnt;j++) {
			id1=atom[i].cnt[j];
			if(atom[i].id<atom[id1].id) {
				bond[k].id=k;
				bond[k].atm[0]=atom[i].id;
				bond[k].atm[1]=atom[id1].id;
				bond[k].atom[0]=&atom[i];
				bond[k].atom[1]=&atom[id1];
				k++;
			}
		}
	}
}

/* ---------------------------------------------------------------------- */
int MODEL::rd_grp(char *ingrpf)
{
	int i,j,id,k,atom1,atom2,is_range;
	char null[1024],dummy[1024];

	istringstream ss;
	if(strcmp(ingrpf,"none")==0) {
		i=0;
		ngrp=1;
		init_grp();
		grp[i].natom=natom;
		grp[i].init_atom();
		for(j=0;j<grp[i].natom;j++) grp[i].atm[j]=j;
		cal_grp_mass();
	} else {
		cout<<"Reading group file "<<ingrpf<<endl;
		ifstream inf(ingrpf,ios::in);
		if(!inf.is_open()) {
			cout<<"Error: group file "<<ingrpf<<" cannot be opened."<<endl;
			i=0;
			ngrp=1;
			init_grp();
			grp[i].natom=natom;
			grp[i].init_atom();
			for(j=0;j<grp[i].natom;j++) grp[i].atm[j]=j;
			cal_grp_mass();

			return filenotopen;
		}
		do { inf.getline(null,1024); } while(! strstr(null,"Total Groups:") && ! inf.eof());
		ss.clear();
		ss.str(null);
		ss>>null>>null>>null;
		ngrp=(int)atof(null);
		if(! ngrp) { cerr<<"Error: Invalid number of groups ("<<null<<") specified"<<endl; exit(1); }
		ngrp++; //last group would be the whole system
		init_grp();
		for(i=0;i<ngrp-1;i++) {
			sprintf(dummy,"Group %i",i+1);
			do { inf.getline(null,1024); } while(! strstr(null,dummy)&&!inf.eof());
			if(inf.eof()) { cerr<<"Error: Cannot find atoms for group "<<i+1<<endl; exit(1); };
			ss.clear();
			ss.str(null);
			// assuming line is "Group x Atoms y"
			ss>>null>>null;
			id=(int)atof(null);
			if(! id) { cerr<<"Error: Could not read group id. Line: "<<null<<endl; exit(1); }
			id--;
			ss>>null>>null;
			grp[id].natom=(int)atof(null);
			if(! grp[id].natom) { 
				cerr<<"Error: Cannot read number of atoms for group "<<id<<" line: "<<(int)atof(null)<<endl; 
				exit(1); 
			}
			grp[id].init_atom();
			j=atom1=is_range=0;
			while(j<grp[id].natom&&!inf.eof()) {
				inf.getline(dummy,1024);
				ss.clear();
				ss.str(dummy);
				while(ss>>null) {
					//cout<<"read: "<<null<<" j "<<j<<" natom "<<grp[id].natom<<endl;
					if(! strcmp(null,"-")) { //range specified
						if (atom1<0) {
							cerr<<"Error: Range specified but no atom!"<<endl;
							exit(1);
						}
						is_range = 1;
					} else {
						if(! is_range) {
							atom1=(int)atof(null)-1;
							if(atom1<0) { 
								cerr<<"Error: Invalid atom ("<<null<<") specified for group "<<i+1<<endl;
								exit(1); 
							}
							grp[id].atm[j]=atom1;
							atom[atom1].grp = id;
							j++;
						} else {
							atom2=(int)atof(null)-1;
							if(! atom2) { 
								cerr<<"Error: Invalid atom ("<<null<<") specified for group "<<i+1<<endl;
								exit(1); 
							}
							if(atom1>atom2) SWAP(&atom1,&atom2);
							for(k=atom1+1;k<=atom2;k++) {
								if(j==grp[id].natom) {
									cerr<<"Error: # atom specified for group "<<i+1<<" is "<<grp[id].natom<<" but read "<<j+1<<"!"<<endl;
									exit(1);
								}
								grp[id].atm[j] = k;
								//cout<<"id "<<k<<" of "<<grp[id].natom<<" in grpid "<<id<<" atom_mass: "<<atom[k].mass<<endl;
								atom[k].grp = id;
								j++;
							}
							is_range = 0;
						}
					}
				}
			}
			cout<<"    Group "<<(id+1)<<" atoms specified: "<<grp[id].natom<<" atoms read: "<<j<<endl;
			if(j != grp[id].natom) cerr<<""<<endl;
		}
		id=ngrp-1;
		grp[id].natom=natom;
		grp[id].init_atom();
		for(j=0;j<grp[id].natom;j++) grp[id].atm[j]=j;
		cal_grp_mass();
		j = 0;
		inf.clear();
		inf.seekg(ios::beg);
		do{ inf>>null; } while (strcmp(null,"Constraints") && ! inf.eof());
		if(!inf.eof()) {
			cout<<"Constraints:"<<endl;
			for(i=0;i<ngrp-1;i++) { 
				inf>>grp[i].constraint; 
				j+= grp[i].constraint; 
				cout<<"    Group "<<(i+1)<<" "<< grp[i].constraint <<endl;
			}
			grp[ngrp-1].constraint = j;
			if(inf.eof()) {
				cout<<"Error: group file Constraints format error"<<endl;
				return 1;
			}
		}

		inf.clear();
		inf.seekg(ios::beg);
		int rotsym_tot = 0;
		do{inf>>null;} while (strcmp(null,"RotationalSymmetryNumber") && ! inf.eof());
		if(!inf.eof()) {
			cout<<"Rotational Symmetry:"<<endl;
			for(i=0;i<ngrp-1;i++) {
				inf>>null; 
				grp[i].rotsym=atof(null);
				cout<<"    Group "<<(i+1)<<" "<< grp[i].rotsym<<endl;
				rotsym_tot += grp[i].rotsym;
			}
			if((rotsym_tot%(ngrp-1) == 0) && ((rotsym_tot/(ngrp-1))==grp[0].rotsym)) 
				grp[ngrp-1].rotsym = grp[0].rotsym;
			else grp[ngrp-1].rotsym = 1;
				cout<<"    Group "<<ngrp<<" "<<grp[ngrp-1].rotsym<<endl;
			if(inf.eof()) {
				cout<<"Error: group file Rotational Symmetry format error"<<endl;
				return 1;
			}
		}

		inf.clear();
		inf.seekg(ios::beg);
		do{inf>>null;} while (strcmp(null,"LinearMoleculeFlag")!=0&&!inf.eof());
		if(!inf.eof()) {
			for(i=0;i<ngrp-1;i++) inf>>grp[i].linear;
			if(inf.eof()) {
				cout<<"Error: group file linear format error"<<endl;
				return 1;
			}
		}
		grp[ngrp-1].linear = 0;

		inf.clear();
		inf.seekg(ios::beg);
		do{inf>>null;} while (strcmp(null,"GroupEnergyAvg")!=0&&!inf.eof());
		if(!inf.eof()) {
			cout<<"GroupEnergyAvg:"<<endl;
			for(i=0;i<ngrp-1;i++) { 
				inf>>grp[i].eng.sum;
				grp[i].eng.avg=grp[i].eng.sum; 
				grp[i].eng.nXY = grp[i].eng.npt = 1; 
				cout<<"    Group "<<(i+1)<<" "<< grp[i].eng.sum <<endl;
			}
			if(inf.eof()) {
				cout<<"Error: group file GroupEnergyAvg format error"<<endl;
			 	return 1;
			}
		 }

		inf.clear();
		inf.seekg(ios::beg);
		do{inf>>null;} while (strcmp(null,"GroupEnergyStd")!=0&&!inf.eof());
		if(!inf.eof()) {
			cout<<"GroupEnergyStdev:"<<endl;
			for(i=0;i<ngrp-1;i++) {
				inf>>grp[i].eng.std;
				grp[i].eng.sigma2x=grp[i].eng.std*grp[i].eng.std;
				grp[i].eng.avg=grp[i].eng.sum;
				grp[i].eng.nXY = grp[i].eng.npt = 1; 
				cout<<"    Group "<<(i+1)<<" "<< grp[i].eng.std<<endl;
			}
			if(inf.eof()) {
				cout<<"Error: group file GroupEnergyStd format error"<<endl;
				return 1;
			}
		}

		inf.clear();
		inf.seekg(ios::beg);
		do{inf>>null;} while (strcmp(null,"GroupVolume")!=0&&!inf.eof());
		if(!inf.eof()) {
			cout<<"GroupVolume:"<<endl;
			for(i=0;i<ngrp-1;i++) {
				inf>>grp[i].vol;
				cout<<"    Group "<<(i+1)<<" "<< grp[i].vol<<endl;
			}
			if(inf.eof()) {
				cout<<"Error: group file GroupVolume format error"<<endl;
				return 1;
			}
		}
		
		inf.close();
		cout<<"Finish reading group file "<<ingrpf<<endl;
	}

	return 0;
}

/* ---------------------------------------------------------------------- */
void MODEL::init_grp()
{
	if(ngrp>0) {
		if(grp!=NULL) delete [] grp;
		grp= new GROUP [ngrp];
	}
}

/* ---------------------------------------------------------------------- */
void MODEL::init_mol()
{
	if(nmol>0) {
		if(mol!=NULL) delete [] mol;
		mol= new MOLECULE [nmol];
	}
}

/* ---------------------------------------------------------------------- */
void MODEL::init_atom()
{

	if(natom>0)	{
		if(atom!=NULL) delete [] atom;
		atom = new ATOM [natom];
	}
}

/* ---------------------------------------------------------------------- */
void MODEL::init_bond()
{
	if(nbond>0) {
		if(bond!=NULL) delete [] bond;
		bond= new BOND [nbond];
	}
}

/* ---------------------------------------------------------------------- */
void MODEL::cal_mass()
{
	int i;
	prp.mass=0;
	for(i=0;i<natom;i++) {
		prp.mass += atom[i].mass;
	}
	for(i=0;i<nmol;i++) {
		mol[i].cal_mass();
	}
	cout<<"total mass "<<prp.mass<<"...";
}

/* ---------------------------------------------------------------------- */
void MODEL::cal_grp_mass()
{
	int i,j;
	for(i=0;i<ngrp;i++) {
		grp[i].mass=0;
		for(j=0;j<grp[i].natom;j++) grp[i].mass += atom[grp[i].atm[j]].mass;
	}
}

/* ---------------------------------------------------------------------- */
void MODEL::ck_range()
{
	int i,k;
	double min[3],max[3];

	if (periodic==0) {
		for(k=0;k<3;k++) { 
			cell.o[k]=0.0;
			max[k]=-1e99;
			min[k]=1e99;
		}
		for(i=0;i<natom;i++) {
			for(k=0;k<3;k++) {
				if(atom[i].pv[k]>max[k]) max[k]=atom[i].pv[k];
				else if(atom[i].pv[k]<min[k]) min[k]=atom[i].pv[k];
			}
		}
		cell.la=max[0]-min[0]; cell.alpha=90;
		cell.lb=max[1]-min[1]; cell.beta=90;
		cell.lc=max[2]-min[2]; cell.gamma=90;
		cell.labc2H();
		cell.H2others();
	}
}

/* ---------------------------------------------------------------------- */
void MODEL::cal_vcomp()
{
	int i;

	for(i=0;i<nmol;i++) mol[i].cal_vc();
}

/* ---------------------------------------------------------------------- */
void MODEL::cal_T()
{
	int i,j;
	double v2,df;

	prp.T=0;
	for(i=0;i<natom;i++) {
		v2=0;
		for(j=0;j<3;j++) v2+= (atom[i].vel[j]*atom[i].vel[j]);
		prp.T += (atom[i].mass*v2);
	}
	df=3.0*natom-(periodic>0?3:0);
	prp.T /= (df*R*0.1);	//0.1 is due to unit conversion
}
 
/* ---------------------------------------------------------------------- */
void MODEL::find_cm()
{
	int i,k;

	prp.cmpv[0]=prp.cmpv[1]=prp.cmpv[2]=0;
	for(i=0;i<natom;i++)
		for(k=0;k<3;k++) prp.cmpv[k] += atom[i].mass*atom[i].pv[k];
	
	for(k=0;k<3;k++) prp.cmpv[k]/=prp.mass;
}

/* ---------------------------------------------------------------------- */
void MODEL::find_molcm()
{
	int i;

	for(i=0;i<nmol;i++) 
		mol[i].find_cm();
}

/* ---------------------------------------------------------------------- */
void MODEL::find_molingrp()
{
	int i,j,k,l,match;

	for(i=0;i<ngrp;i++) {
		grp[i].nmol=0;
		for(j=0;j<nmol;j++) {
			match=0;
			for(k=0;k<mol[j].natom;k++) 
			for(l=0;l<grp[i].natom;l++) if(grp[i].atm[l]==mol[j].atm[k]) {match++;break;}
	
			if(match==mol[j].natom) grp[i].nmol++;
		}
		grp[i].mol=new int [grp[i].nmol];
		grp[i].nmol=0;
		for(j=0;j<nmol;j++) {
			match=0;
			for(k=0;k<mol[j].natom;k++) 
				for(l=0;l<grp[i].natom;l++) if(grp[i].atm[l]==mol[j].atm[k]) {match++;break;}
		
			if(match==mol[j].natom) grp[i].mol[ grp[i].nmol++ ]=j;
		}
	}
}

/* ---------------------------------------------------------------------- */
void MODEL::init_trj()
{
	int i;

	ctl.in_trj_flag=trj.init_trj(ctl.in_trj_n,ctl.in_trj,ctl.in_trj_flag,0,&cell,atom,&prp,natom);

	//check if the model and trj agree
	for(i=0;i<ctl.in_trj_n;i++) {
		if(natom!=trj.strj[i].header.nmovatm) {
			cout<<"Error Mismatch of atoms in model and trj # "<<i+1<<": "<<natom<<" "<<trj.strj[i].header.nmovatm<<endl;
			exit(1);
		}
	}

	if(ctl.ana_fframe==0||ctl.ana_fframe>trj.totframe) ctl.ana_fframe=trj.totframe;
	if(ctl.ana_iframe< 0||ctl.ana_iframe>trj.totframe) ctl.ana_iframe=1;
	if(ctl.ana_sframe<=0) ctl.ana_sframe=1;

	for(i=0;i<ctl.in_trj_n;i++) {
		ctl.in_trj[i].molopt = ctl.ana_vac_2pt;
		trj.strj[i].header.timestep=ctl.ana_vac_tstep;
	}
}

/* ---------------------------------------------------------------------- */
void MODEL::element2prp(int do_mass)
{
	int i,j,found;
	string atmname, fftype;
	ele *elements;
	char chars[] = "-0123456789.";

	elements = loadelements();

	for(i=0;i<natom;i++) {
		atmname.assign(toLowerCase(atom[i].name));
		fftype.assign(toLowerCase(atom[i].fftype));
		//remove _xx
		if((j=atmname.find("_",0))>-1) {atmname.erase(j,atmname.length()-j); }
		if((j=fftype.find("_",0))>-1) {fftype.erase(j,fftype.length()-j); }
		//remove all non-letters
		for(j=0;j<strlen(chars);++j) {
			atmname.erase(std::remove(atmname.begin(), atmname.end(), chars[j]), atmname.end());
			fftype.erase(std::remove(fftype.begin(), fftype.end(), chars[j]), fftype.end());
		}
		j=1;
		found = 0;
		while(j<110) {
			if(elements[j].symbol == fftype || elements[j].symbol == atmname) {
				atom[i].mass=elements[j].mass;
				found=1;
				break;
			}
			j++;
		}
		if(!found) {
			cout<<"ERROR: Cannot locate any valid element for atom #"<<i+1<<": "<<atom[i].name<<" (fftype: "<<atom[i].fftype<<")"<<endl;
			exit(1);
		}
	}
}

