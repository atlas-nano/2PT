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

//read atom trajectories

#include <iostream>
#include <vector>
#include <sstream>
#include "trj_header.h"
//#include <sys/mman.h>
//#include <sys/stat.h>

/* ---------------------------------------------------------------------- */
TRJheader::TRJheader() {

	mvatmofst=NULL;
	lmp_data_len = 3;
	byte_offset = NULL;
	xyz_flag=vel_flag=stress_flag=image_flag=force_flag=charge_flag=eng_flag=0; //dump format flags
}

/* ---------------------------------------------------------------------- */
TRJheader::~TRJheader() 
{

	if(mvatmofst!=NULL) delete [] mvatmofst;
	if(byte_offset!=NULL) delete [] byte_offset;
}

/* ---------------------------------------------------------------------- */
int TRJheader::init_header()
{
	int i;

	if(nmovatm==0) return 1;
	if(mvatmofst!=NULL) delete [] mvatmofst;
	mvatmofst=new int [nmovatm];

	for(i=0;i<nmovatm;i++) { mvatmofst[i]=i+1; }
	mvatmpfu[0]=nmovatm;   //number of movable atoms
	natmpfu[0]=nmovatm;		//total number of atoms

	return 0;
}

/* ---------------------------------------------------------------------- */
void TRJheader::CleanTRJheader(int clearatom = 1)
{
	int i,j;

	for(i=0;i<5;i++) hdr[i]='\0';  
	for(i=0;i<20;i++) icntrl[i]=0; 
	ntrjti=neexti=0;
	for(i=0;i<10;i++) { 
		for(j=0;j<80;j++) trjtic[i][j]=eextic[i][j]='\0'; 
	}
	period=molxtl=lcanon=defcel=prtthrm=lnose=lnpecan=ltmpdamp=nflusd=0;
	lmp_cell_format = 0;
	for(i=0;i<MAXFILE;i++) { mvatmpfu[i]=natmpfu[i]=0; for(j=0;j<8;j++) decusd[i][j]='\0'; }
	if (clearatom) natom=nmovatm=movatm1=movatmn=totmov=0;
	//for(i=0;i<nmovatm;i++) { mvatmofst[i]=0; }
	leexti=lparti=version=0;
	for(i=0;i<80;i++) eextit[i]=partit[i]='\0';
	xyzfreq=velfreq=engfreq=strfreq=0;
	nxyz=nvel=neng=nstr=0;
	timestep=0;
	totaccustep=0;

	strcpy(hdr,"MDTR");
	version=icntrl[0]=2010;
	ntrjti=1; strcpy(trjtic[0],"COMET trajectory");
	nflusd=1;
	strcpy(decusd[0],"TEST");
	leexti=7;
	strcpy(eextit,"NOTITLE");
	lparti=5;
	strcpy(partit,"NOPAR");
		 
}

/* ---------------------------------------------------------------------- */
int TRJheader::ReadBinHeader(ifstream *intrj)
{
	int i;
	float null;  //dummy variable, one carriage return
	double crtn;

	CleanTRJheader();
		
//  ----- Read header information -----
	intrj->read((char *)&null,sizeof(null)); 
	intrj->read((char *)&hdr,sizeof(hdr)-1);

	for(i=0;i<20;i++) intrj->read((char *)&icntrl[i],sizeof(icntrl[i])); 
	intrj->read((char *)&crtn,sizeof(crtn)); //read carriage return produced in unformatted Fortran file

	version=icntrl[0];

	if (version>=311) {
		intrj->read((char *)&ntrjti,sizeof(ntrjti));
		for (i=0;i<ntrjti;i++)		intrj->read((char *)&trjtic[i],sizeof(trjtic[i]));
	}
	intrj->read((char *)&crtn,sizeof(crtn));

	intrj->read((char *)&neexti,sizeof(neexti));
	for (i=0;i<neexti;i++) intrj->read((char *)&eextic[i],sizeof(eextic[i]));
	intrj->read((char *)&crtn,sizeof(crtn));

	if (version<=150) {
		period=false;
		molxtl=false;
		lcanon=false;
		defcel=false;
		prtthrm=false;
	} else if (version<300) {
		intrj->read((char *)&period,sizeof(period));
		intrj->read((char *)&molxtl,sizeof(molxtl));
		intrj->read((char *)&lcanon,sizeof(lcanon));
		intrj->read((char *)&defcel,sizeof(defcel));
		intrj->read((char *)&prtthrm,sizeof(prtthrm));
		intrj->read((char *)&crtn,sizeof(crtn));
		lnose =false;
		lnpecan=false;
		ltmpdamp=false;
	} else {
		intrj->read((char *)&period,sizeof(period)); //Periodicity
		intrj->read((char *)&molxtl,sizeof(molxtl)); //MolXtl
		intrj->read((char *)&lcanon,sizeof(lcanon)); //Canonical
		intrj->read((char *)&defcel,sizeof(defcel)); //DefCel
		intrj->read((char *)&prtthrm,sizeof(prtthrm)); //PertTheory
		intrj->read((char *)&lnose,sizeof(lnose)); //NoseOrHoover
		intrj->read((char *)&lnpecan,sizeof(lnpecan)); //NpTCanon
		intrj->read((char *)&ltmpdamp,sizeof(ltmpdamp)); //TempDamping
		intrj->read((char *)&crtn,sizeof(crtn));
	}

	if (version>=200) {
		intrj->read((char *)&nflusd,sizeof(nflusd));
		for(i=0;i<nflusd;i++) intrj->read((char *)&mvatmpfu[i],sizeof(mvatmpfu[i]));
		for(i=0;i<nflusd;i++) intrj->read((char *)&natmpfu[i],sizeof(natmpfu[i]));
		for(i=0;i<nflusd;i++) intrj->read((char *)&decusd[i],sizeof(decusd[i]));
		intrj->read((char *)&crtn,sizeof(crtn));
		   
		intrj->read((char *)&totmov,sizeof(totmov));
		nmovatm = totmov;
		init_header();

		for(i=0;i<totmov;i++) intrj->read((char *)&mvatmofst[i],sizeof(mvatmofst[i]));
		intrj->read((char *)&crtn,sizeof(crtn));

	} else {
		intrj->read((char *)&natom,sizeof(natom));
		 intrj->read((char *)&nmovatm,sizeof(nmovatm));
		intrj->read((char *)&movatm1,sizeof(movatm1));
		intrj->read((char *)&movatmn,sizeof(movatmn));
		intrj->read((char *)&crtn,sizeof(crtn));
	}

	intrj->read((char *)&leexti,sizeof(leexti));
	for(i=0;i<leexti;i++) intrj->read((char *)&eextit[i],sizeof(eextit[i]));
	intrj->read((char *)&crtn,sizeof(crtn));

	intrj->read((char *)&lparti,sizeof(lparti));
	for(i=0;i<lparti;i++) intrj->read((char *)&partit[i],sizeof(partit[i]));
	intrj->read((char *)&crtn,sizeof(crtn));

// ----- end of header information -----
// intrj->read((char *)&null,sizeof(null));
	return 0;
}

/* ---------------------------------------------------------------------- */
int TRJheader::ReadAscHeader(ifstream *intrj)
{
	int i; 
	char check[1024];

	CleanTRJheader();

//  ----- Read header information -----
	*intrj>>check>>hdr;
	if(strcmp(check,"MPSim-trj_ASCII_file_for_BINARY_conversion")!=0) {
		cout<<" ASCII file does not have the right format. Reading fails."<<endl;
		return filenotopen;
	}
	for(i=0;i<20;i++) *intrj>>icntrl[i];
	version=icntrl[0];

	if (version>=311) {
		*intrj>>ntrjti;
		for (i=0;i<ntrjti;i++) { 
			intrj->getline(trjtic[i],sizeof(trjtic[i]));
			intrj->clear();
		}
	}
	*intrj>>neexti;
	for (i=0;i<neexti;i++) {
		intrj->getline(eextic[i],sizeof(eextic[i]));
		intrj->clear();
	}
	if(version<=150) ;
	else if (version<300) {
		*intrj>>period>>molxtl>>lcanon>>defcel>>prtthrm;
	} else {
		*intrj>>period>>molxtl>>lcanon>>defcel>>prtthrm>>lnose>>lnpecan>>ltmpdamp;
	}
	if (version>=200) {
		*intrj>>nflusd;
		for(i=0;i<nflusd;i++) *intrj>>mvatmpfu[i];
		for(i=0;i<nflusd;i++) *intrj>>natmpfu[i];
		for(i=0;i<nflusd;i++) *intrj>>decusd[i];
		*intrj>>totmov;
		nmovatm = totmov;
		init_header();
		for(i=0;i<totmov;i++) *intrj>>mvatmofst[i];
	} else {
		*intrj>>natom>>nmovatm>>movatm1>>movatmn;
	}

	*intrj>>leexti;
	for(i=0;i<leexti;i++) *intrj>>eextit[i];
	*intrj>>lparti;
	for(i=0;i<lparti;i++) *intrj>>partit[i];

//  ----- end of header information -----
	return 0;
}

/* ---------------------------------------------------------------------- */
int TRJheader::ReadUSRHeader(ifstream *intrj)
{
	CleanTRJheader();
	intrj->read((char *)&totmov,sizeof(totmov));
	nmovatm = totmov;
	init_header();
	cout<<" usrtrj natom "<<nmovatm<<endl;

	return 0;
}

/* ---------------------------------------------------------------------- */
int TRJheader::ReadXYZHeader(ifstream *intrj, int iflag)
{
	int i,prev,fstp,totframe;
	char null[1024];

	intrj->clear();
	intrj->seekg(0,ios::beg);

	*intrj>>fstp; i=fstp;
	totmov=0;

	ntrjti=1; strcpy(trjtic[0],"XYZ trajectory");
	period=3;
	init_header();
	nmovatm = totmov = i;
	if(iflag)cout<<" xyzcoord natom "<<nmovatm<<"...";
	else cout<<" xyzvel natom "<<nmovatm<<"...";
	fflush(stdout);

	totframe=1;
	intrj->seekg(0,ios::beg);
	while(!intrj->eof()) {
		intrj->getline(null,1024);
		if(strstr(null,"i =")) totframe++;
	}
	init_header();

	if(byte_offset == NULL) byte_offset = new byteoffset[totframe];
	nxyz = totframe;
	intrj->clear();
	intrj->seekg(0, ios::beg);
	totframe = 0;
	prev = 0;
	while(!intrj->eof()) {
		intrj->getline(null,1024);
		if(strstr(null," i =")) { 
			if(!iflag) byte_offset[totframe].xyz = prev;
			else byte_offset[totframe].coord = prev;
			totframe++; 
		} else prev = intrj->tellg();
	}
	lmptmp = totframe;
	cout<<" nframe "<<totframe<<"...";

	return 0;
}

/* ---------------------------------------------------------------------- */
int TRJheader::ReadAMBERHeader(ifstream *intrj, int iflag)
{
	int i,j,k; 
	int totframe,amberbox,invalid,tot,chunk;
	long long soffset,eoffset,slen,flen;
	char null[1024],junk[8];
	float crtn;

	istringstream ss;


	amberbox = 0;
	chunk = 8 * sizeof(char);
	intrj->clear();
	intrj->seekg(0,ios::end);
	flen = intrj->tellg();
	intrj->seekg(0,ios::beg);
  
	//first get the title
	intrj->getline(null,1024);
	soffset = intrj->tellg(); //starting atom offset
	ntrjti=1; strcpy(trjtic[0],null);
	period=3;
	init_header();

	if(iflag) cout<<"coord...";
	else cout<<"vels...";

	//now get size of 1 snapshot
	i = 0;
	j = 0;
	invalid = 0;
	ss.width(8);
	while (i<totmov) {
		intrj->getline(null,1024);
		ss.clear();
		ss.str(null);
		while (ss.good() && ! ss.eof()) {
			ss.read(junk,chunk);
			j++;
			if(j % 3 == 0) i++;
			if(j % 10 ==0) ss.read((char *)&crtn,sizeof(crtn));
			if (i==totmov) break;
		}
		if(! intrj->good() || intrj->eof()) { invalid = 1; break; }
	}
	if (invalid) return invalid;

	//first check to see if we get an integer number of frames
	eoffset = intrj->tellg();
	tot = totmov * 3;
	if(ss.read(junk,chunk) && ss.good() && ! ss.eof()) { //have box info apparently, so advance past next line if necessary
		j++; k = 1;
		while (ss.good() && ! ss.eof()) {
			ss.read(junk,chunk);
			j++; k++;
			if(j % 10 ==0) ss.read((char *)&crtn,sizeof(crtn));
		}
		if(j % 10 == 0) {
			//look to see if the box info extends into the next line
			intrj->getline(null,1024);
			ss.clear();
			ss.str(null);
			while (ss.good() && ! ss.eof()) {
				ss.read(junk,chunk);
				j++; k++;
				if(j % 10 ==0) ss.read((char *)&crtn,sizeof(crtn));
			}
		}
		if(k && k<10 && k != totmov) { //probably a box
			amberbox = 1;
			eoffset = intrj->tellg();
			tot += k;
		}

	} else { //maybe not, so read next line and see if less then 10 entries
		intrj->getline(null,1024);
		ss.clear();
		ss.str(null);
		k = 0;
		while (ss.good() && ! ss.eof()) {
			ss.read(junk,chunk);
			k++;
			if(k % 10 ==0) ss.read((char *)&crtn,sizeof(crtn));
		}
		if(k && k<10 && k != totmov) { //probably a box
			eoffset = intrj->tellg();
			amberbox = 1;
			tot += k;
		}
	}
	//now determine number of frames
	slen = eoffset - soffset;
	if((flen - soffset) % slen > 0) { //trajectory is weird so quite
		cout<<endl<<"ERROR: for 1st frame: byte start "<<soffset<<" end "<<eoffset<<" len "<<slen;
		cout<<" file size: "<<flen<<" so non integer #frames"<<endl;
		return countmismatch;
	}
	if(amberbox) cout<<"box..";
	totframe = (flen - soffset)/slen;
	if(byte_offset == NULL) byte_offset = new byteoffset[totframe+1];
	intrj->clear();
	intrj->seekg(soffset,ios::beg); //back to start and store offset of 1st frame
	intrj->width(8);
	if(iflag) {
		byte_offset[0].coord = soffset;
		byte_offset[0].camberbox = amberbox;
	} else {
		byte_offset[0].xyz = soffset;
		byte_offset[0].vamberbox = amberbox;
	}
	//all other frames
	for (i=1;i<=totframe;i++) {
		if(iflag) {
			byte_offset[i].coord = soffset + i*slen;
			byte_offset[i].camberbox = amberbox;
		} else {
			byte_offset[i].xyz = soffset + i*slen;
			byte_offset[i].vamberbox = amberbox;
		}
  }

	lmptmp = totframe;
	cout<<"nframe "<<totframe<<"...";

	return 0;
}

/* ---------------------------------------------------------------------- */
int TRJheader::ReadCPMDHeader(ifstream *intrj,ifstream *ineng,ifstream *instr)
{
	int i,fstp;
	char null[1024];
	unsigned long tmp;

	CleanTRJheader();

	tmp=intrj->tellg();
	intrj->clear();
	intrj->seekg(0,ios::beg);
	//find number of atoms
	*intrj>>fstp; i=fstp;
	totmov=0;
	while(i==fstp) {
		totmov++;
		intrj->getline(null,1024);
		*intrj>>i;
	}
	nmovatm = totmov;  //number of movable atoms
	init_header();

	ntrjti=1; strcpy(trjtic[0],"CPMD trajectory");
	period=3;

	//determine xyzfreq, velfreq
	xyzfreq=velfreq=i-fstp;
	intrj->seekg(tmp, ios::beg);

	nxyz=nstr=neng=xyzfreq;
	return 0;
}

/* ---------------------------------------------------------------------- */
int TRJheader::ReadLMPHeader(ifstream *intrj, ifstream *ineng, ATOM * atom)
{
	int i,j,tstart,tend,tcurr,tprev;
	char line[1024],data[1024];
	long int tmp1,tmp2,tmp3,tmp4;
	int max;
	int tstep_atoms; 
	long long byte_per_line, offset;

	tmp1=tmp2=tmp3=tmp4=0;
	max = 10000000;
	long int *trjByte = new long int [max];
	long int *engByte = new long int [max];
	istringstream ss;

	cout<<"   Getting LAMMPS Header info...";
	fflush(stdout);


	ntrjti=1; strcpy(trjtic[0],"LAMMPS trajectory");
	period=3;
	tprev=-1;
	CleanTRJheader();

	intrj->clear();
	tmp1=intrj->tellg(); //find out where we are in the file
	intrj->seekg(0,ios::beg); //go to the beginning

	//find first timestep
	intrj->getline(line,1024);
	while(! strstr(line,"ITEM: TIMESTEP")) { intrj->getline(line,1024); };
	*intrj>>tstart;
	tend = tstart;
	trjByte[0] = 0;

	//find number of atoms
	while(! strstr(line,"ITEM: NUMBER OF ATOMS")) { intrj->getline(line,1024); };
	*intrj>>i;
	if (! i) {
		cerr<<" Error: Number of atoms cannot be determined from TRAJ/VEL file."<<endl;
		return countmismatch;
	}
	nmovatm = totmov = tstep_atoms = i;  //number of movable atoms
	init_header();

	//determine the dump format of the box
	while(! strstr(line,"ITEM: BOX BOUNDS")) { intrj->getline(line,1024); };
	if(strstr(line,"xy xz yz")) lmp_cell_format = 1;

	//determine dump format and set options
	while(! strstr(line,"ITEM: ATOMS")) { intrj->getline(line,1024); };
	ss.clear();
	ss.str(line); //copy the line into a string stream
	ss>>data>>data;
	i = 0;
	while(ss>>data) i++;
	if (! i) SetDefLMPFormat(atom,nmovatm);
	else {
		lmp_data_len = i;
		lmp_data_index = (double ***) calloc(nmovatm,sizeof(double **));
		for(i=0;i<nmovatm;i++)
			lmp_data_index[i] = (double **) calloc(lmp_data_len,sizeof(double *));
		//read the format specification line again and save pointer accordingly
		ss.clear();
		ss.str(line); //copy the line into a string stream
		ss>>data>>data;
		j=0;
		while(ss>>data) { //read arguments 1 at a time and create map
			if(!strcmp(data,"x")||!strcmp(data,"xs")||!strcmp(data,"xu")) {
				xyz_flag = 1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].pv[0];
			} else if(!strcmp(data,"y")||!strcmp(data,"ys")||!strcmp(data,"yu")) {
				xyz_flag = 1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].pv[1];
			} else if(!strcmp(data,"z")||!strcmp(data,"zs")||!strcmp(data,"zu")) {
				xyz_flag = 1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].pv[2];

			} else if(!strcmp(data,"vx")) {
				vel_flag = 1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].vel[0];
			} else if(!strcmp(data,"vy")) {
				vel_flag = 1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].vel[1];
			} else if(!strcmp(data,"vz")) {
				vel_flag = 1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].vel[2];

			} else if(!strcmp(data,"fx")) {
				force_flag = 1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].f[NET][0];
			} else if(!strcmp(data,"fy")) {
				force_flag = 1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].f[NET][1];
			} else if(!strcmp(data,"fz")) {
				force_flag = 1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].f[NET][2];

			} else if(!strcmp(data,"ix")) {
				image_flag=1;  
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].img[0];
			} else if(!strcmp(data,"iy")) {
				image_flag=1;  
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].img[1];
			} else if(!strcmp(data,"iz")) {
				image_flag=1;  
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].img[2];

			} else if(!strcmp(data,"c_stress[1]")) {
				stress_flag=1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].strs[0];
			} else if(!strcmp(data,"c_stress[2]")) {
				stress_flag=1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].strs[1];
			} else if(!strcmp(data,"c_stress[3]")) {
				stress_flag=1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].strs[2];
			} else if(!strcmp(data,"c_stress[4]")) {
				stress_flag=1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].strs[3];
			} else if(!strcmp(data,"c_stress[5]")) {
				stress_flag=1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].strs[4];
			} else if(!strcmp(data,"c_stress[6]")) {
				stress_flag=1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].strs[5];

			} else if(!strcmp(data,"q")) {
				charge_flag = 1;
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = &atom[i].chg;
			} else if(!strcmp(data,"type")) {
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = (double *)&atom[i].ffid;
			} else if(!strcmp(data,"v_atomEng")) {
				eng_flag = 1;
				for(i=0;i<nmovatm;i++) 
					lmp_data_index[i][j] = &atom[i].eng_tmp;
			} else { //key not found so assign null
				for(i=0;i<nmovatm;i++)
					lmp_data_index[i][j] = 0;
			}
			j++;
		}
	}

	//determine nxyz nvel nstr
	intrj->getline(line,1024);
	tmp3 = intrj->tellg(); //start of atom section
	while(intrj->good() && ! intrj->eof()) {
		tmp4 = intrj->tellg(); //end of atom section
		intrj->getline(line,1024);
		if(strstr(line,"ITEM: TIMESTEP")) {
			trjByte[1] = tmp4;
			*intrj>>tend;
			break;
		}
	}

	byte_per_line = ((tmp4-tmp3)/nmovatm)-10; //crude
	if(byte_per_line < 10) byte_per_line = 10;
	offset = (tstep_atoms)*byte_per_line;
	intrj->seekg(offset,ios::cur);
	offset = 0;
	i = 2;
	while(intrj->good() && ! intrj->eof()) {
		tmp4 = intrj->tellg();
		intrj->getline(line,1024);
		if(strstr(line,"ITEM: TIMESTEP")) {
			trjByte[i] = tmp4;
			i++;
			*intrj>>tend;
		} else if(strstr(line,"ITEM: NUMBER OF ATOMS")) { //this should be constant, but fix here work if it isn't
			*intrj>>tstep_atoms;
		} else if(strstr(line,"ITEM: ATOMS")) {
			offset = (tstep_atoms)*byte_per_line;
			intrj->seekg(offset,ios::cur); //navigate to previous line of new timestep
		}
	}
	if (! i) {
		cerr<<" Error: Cannot determine number of frames in TRAJ/VEL file."<<endl;
		return countmismatch;
	}
	nxyz=nvel=nstr=i;

	if(ineng->is_open()) {
		//determine neng
		tmp2=ineng->tellg();
		ineng->seekg(0,ios::beg);
		i=0;
		//navigate to the 1st timestep from xyz file
		while(! ineng->eof()) {
			ineng->getline(line,1024);
			if(strstr(line,"-- Step")) {
				sscanf(line,"%*s %*s %d",&tcurr);
				if (tcurr>=tstart && tcurr<=tend && tcurr>tprev) {
					engByte[i] = ineng->tellg();
					i++; 
					tprev=tcurr;
				}
			}
		}
		if (! i) {
			cerr<<" Error: Cannot determine number of frames in ENG file."<<endl;
			return countmismatch;
		}
		neng=i;
		if(nxyz>neng) {
			cerr<<" Error: Number of snapshots exceeds the number of data in the thermo file."<<endl;
			return countmismatch;
		}

		//determine the skips of each trj, not fully optimized yet
		xyzfreq=velfreq=engfreq=strfreq=1;
		i = (int)(nxyz/neng);
		if(i > 1) engfreq = i;
	}

	byte_offset = new byteoffset[nxyz+1];
	for(i=0;i<nxyz;i++) {
		byte_offset[i].xyz = trjByte[i];
		if (ineng->is_open()) byte_offset[i].eng = engByte[i*engfreq];
	}
	intrj->clear();
	intrj->seekg(0,std::ios::end);
	byte_offset[nxyz].xyz = intrj->tellg();
	intrj->seekg(tmp1,std::ios::beg);
	if (ineng->is_open()) {
		ineng->clear();
		ineng->seekg(0,std::ios::end);
		byte_offset[nxyz].eng = ineng->tellg();
		ineng->seekg(tmp2,std::ios::beg);
	}
	totaccustep=tend;
	lmptmp=nxyz;

	if(xyz_flag) cout<<"coords...";
	if(vel_flag) cout<<"vels...";
	if(stress_flag) cout<<"stress...";
	if(force_flag) cout<<"forces...";
	if(image_flag) cout<<"images...";
	if(charge_flag) cout<<"charge...";

	cout<<"natoms "<<nmovatm<<" nxyz "<<nxyz<<" neng "<<neng<<" ndata "<<lmp_data_len+1<<"...Done"<<endl;

	delete [] trjByte;
	delete [] engByte;

	return 0;
}

/* ---------------------------------------------------------------------- */
void TRJheader::SetDefLMPFormat(ATOM * atom, int natom) {

	int i,j;

	lmp_data_len = 4;
	lmp_data_index = (double ***) calloc(natom,sizeof(double **));
	for (i=0;i<natom;i++) {
		lmp_data_index[i] = (double **) calloc(lmp_data_len,sizeof(double *));
		lmp_data_index[i][0] = (double *)&atom[i].ffid;
		for(j=0;j<3;j++) lmp_data_index[i][j+1] = &atom[i].pv[j];
	}
}

