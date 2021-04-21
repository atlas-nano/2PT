/* Memory management module adopted from LAMMPS MD code (lammps.sandia.gov) */

#include <iostream>
#include <fstream>
#include "memory.h"

/* ----------------------------------------------------------------------
   safe free 
------------------------------------------------------------------------- */

void Memory::sfree(void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);
}

/* ----------------------------------------------------------------------
   safe malloc 
------------------------------------------------------------------------- */

void * Memory::smalloc(double n, const char *name)
{
  if (!(int)n) return NULL;
  void *ptr = malloc(n);
  if (ptr == NULL) {
    cout<<"Failed to allocate "<<n<<" bytes for array "<<name<<endl;
    exit(1);
  }
  return ptr;
}

/* ----------------------------------------------------------------------
   create a 2d int array 
------------------------------------------------------------------------- */

int **Memory::create_2d_int_array(int n1, int n2, const char *name)

{
  int *data = (int *) smalloc(n1*n2*sizeof(int),name);
  int **array = (int **) smalloc(n1*sizeof(int *),name);

  int n = 0;
  for (int i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2;
  }

  return array;
}

/* ----------------------------------------------------------------------
   free a 2d int array 
------------------------------------------------------------------------- */

void Memory::destroy_2d_int_array(double **array)

{
  if (array == NULL) return;
  sfree(array[0]);
  sfree(array);
}

/* ----------------------------------------------------------------------
   create a 2d double array 
------------------------------------------------------------------------- */

double **Memory::create_2d_double_array(int n1, int n2, const char *name)

{
  double *data = (double *) smalloc(n1*n2*sizeof(double),name);
  double **array = (double **) smalloc(n1*sizeof(double *),name);
  
  int n = 0;
  for (int i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2; 
  } 
  
  return array;
} 

/* ----------------------------------------------------------------------
   free a 2d double array 
------------------------------------------------------------------------- */

void Memory::destroy_2d_double_array(double **array)

{
  if (array == NULL) return;
  sfree(array[0]);
  sfree(array);
}

/* ----------------------------------------------------------------------
   create a 4d double array 
------------------------------------------------------------------------- */
  
float **** Memory::create_4d_float_array(int n1, int n2, int n3, int n4,
                                          const char *name)
{
  int i,j,k;
  double s1 = n1*sizeof(float ***);
  double s2 = n1*n2*sizeof(float **);
  double s3 = n1*n2*n3*sizeof(float *);
  double s4 = n1*n2*n3*n4*sizeof(float);
  float *data = (float *) smalloc(s4,name);
  float **cube = (float **) smalloc(s3,name);
  float ***plane = (float ***) smalloc(s2,name);
  float ****array = (float ****) smalloc(s1,name);

  long n = 0;
  for (i = 0; i < n1; i++) {
    array[i] = &plane[i*n2];
    for (j = 0; j < n2; j++) {
      plane[i*n2+j] = &cube[i*n2*n3+j*n3];
      for (k = 0; k < n3; k++) {
        cube[i*n2*n3+j*n3+k] = &data[n];
        n += n4;
      }
    }
  }

  return array;
}

/* ----------------------------------------------------------------------
   free a 4d float array 
------------------------------------------------------------------------- */

void Memory::destroy_4d_float_array(float ****array)
{
  if (array == NULL) return;
  sfree(array[0][0][0]);
  sfree(array[0][0]);
  sfree(array[0]);
  sfree(array);
}

/* ----------------------------------------------------------------------
   create a 3d double array 
------------------------------------------------------------------------- */

double ***Memory::create_3d_double_array(int n1, int n2, int n3,
                                         const char *name)
{
  int i,j;

  double *data = (double *) smalloc(n1*n2*n3*sizeof(double),name);
  double **plane = (double **) smalloc(n1*n2*sizeof(double *),name);
  double ***array = (double ***) smalloc(n1*sizeof(double **),name);

  int n = 0;
  for (i = 0; i < n1; i++) {
    array[i] = &plane[i*n2];
    for (j = 0; j < n2; j++) {
      plane[i*n2+j] = &data[n];
      n += n3;
    }
  }

  return array;
}

/* ----------------------------------------------------------------------
   free a 3d double array 
------------------------------------------------------------------------- */

void Memory::destroy_3d_double_array(double ***array)
{
  if (array == NULL) return;
  sfree(array[0][0]);
  sfree(array[0]);
  sfree(array);
}

/* ----------------------------------------------------------------------
   create a 4d double array 
------------------------------------------------------------------------- */

double ****Memory::create_4d_double_array(int n1, int n2, int n3, int n4,
                                          const char *name)
{
  int i,j,k;

  double *data = (double *) smalloc(n1*n2*n3*n4*sizeof(double),name);
  double **cube = (double **) smalloc(n1*n2*n3*sizeof(double *),name);
  double ***plane = (double ***) smalloc(n1*n2*sizeof(double **),name);
  double ****array = (double ****) smalloc(n1*sizeof(double ***),name);

  int n = 0;
  for (i = 0; i < n1; i++) {
    array[i] = &plane[i*n2];
    for (j = 0; j < n2; j++) {
      plane[i*n2+j] = &cube[i*n2*n3+j*n3];
      for (k = 0; k < n3; k++) {
        cube[i*n2*n3+j*n3+k] = &data[n];
        n += n4;
      }
    }
  }

  return array;
}

/* ----------------------------------------------------------------------
   free a 4d double array 
------------------------------------------------------------------------- */

void Memory::destroy_4d_double_array(double ****array)
{
  if (array == NULL) return;
  sfree(array[0][0][0]);
  sfree(array[0][0]);
  sfree(array[0]);
  sfree(array);
}

