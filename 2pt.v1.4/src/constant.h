#ifndef CONSTANT_H
#define CONSTANT_H

#include <math.h>

using namespace std;

class trjItem {
public:
   char    vel[1024];
   char    coord[1024];
   char    thermo[1024];
   char    atom_eng[1024];
   char    stress[1024];
   int     molopt;
};

template <class SW> SW SWAP(SW *a, SW *b);
template <class SW> SW SWAP(SW *a, SW *b)
{
        SW swap=*a;
        *a=*b;
        *b=swap;
        return 0;
}

typedef int   logical;

const int    maxcnt=6;                  //maximum number of connectivity per atom
const double MEGABYTE=1024*1024;        // 
const double Na=6.0221367E23;           //Avogadro's number
const double PI=3.141592653589793116;
const double eA2debye=4.802;            // 1e x 1A = 4.802 Debye
const double kb=1.380658E-23;           // Boltzmann constant J/K
const double h=6.62606896E-34;          // Plank constant Js
const double R=kb*Na;                   //8.314472;    /gas constant J/mol K
const double vlight=2.99792458E8;        // m/s
const double Bohr2A=0.529177249;        //1 Bohr = 0.529177249 A
const double Hartree2kcalmol=627.5095;   // 1 Hartree = 627.5095 kcal/mol
const double caltoj=4.184;               // 1 cal = 4.184 J
const double au2fs=0.02418884326505;               // 1 au = 0.0242 fs

#define true          1
#define false         0
#define none          0     //none
#define filenotopen   1001  //file cannot be opened
#define unregkey      1002  //unrecognized keyword
#define countmismatch 1003  //mismatch between number of trj and coord files
#define c2trj         2001  //cerius2 trj
#define asctrj        2002  //ASCII trj
#define lmptrj        2004  //LAMMPS trj
#define usrtrj        2005  //USER  trj
#define xyztrj        2006  //XYZ trj
#define ambertrj      2007  //AMBER trj
#define charmmtrj     2008  //CHARMM trj
#define cpmdtrj       2009  //CHARMM trj
#define mol2grp       3001  //create group file for each molecule
#define strtbgf       4001  //read structure from bgf      
#define strtlmp       4002  //read structure from lammps data file
#define strtamber     4003  //read structure from AMBER prmtop file
#define strtcharmm    4004  //read structure from CHARMM psf file
#define VELTYPE                 4
#define vtotal                  VELTYPE
#define vtrans                  VELTYPE-4
#define vangul                  VELTYPE-3
#define vimvib                  VELTYPE-2
#define vrotat                  VELTYPE-1

#define FORCETYPE 9
#define NET FORCETYPE-9
#define BND FORCETYPE-8
#define ANG FORCETYPE-7
#define TOR FORCETYPE-6
#define INV FORCETYPE-5
#define VDW FORCETYPE-4
#define VDL FORCETYPE-3
#define ELS FORCETYPE-2
#define ELL FORCETYPE-1

#ifdef _WIN32
#define M_PI 3.1415926535897932384626433832795
#define M_SQRT2 1.4142135623730950488016887242097
#endif

#ifdef _WIN32
// On linux, I define this in a Makefile.
// On Windows, I like it defined by default.
#define DEBUG 1
#endif

#ifdef DEBUG
   #define ASSERT(exp) { if (!(exp)) { \
      printf( "Assertion error at line %d in %s\n", __LINE__, __FILE__ ); \
      /*sprintf(0, "This should cause a crash blah blah blah blah.");*/ \
      exit( 1 ); \
   }}
   #define ASSERT_IS_EQUAL(a,b) { \
      ASSERT( (b) - 0.0005 < (a) && (a) < (b) + 0.0005 ); \
   }
#else
   #define ASSERT(exp)
   #define ASSERT_IS_EQUAL(a,b)
#endif


inline int ROUND( float x ) {
   return x < 0 ? -(int)(-x+0.5f) : (int)(x+0.5f);
}
inline int ROUND( double x ) {
   return x < 0 ? -(int)(-x+0.5f) : (int)(x+0.5f);
}


#endif

