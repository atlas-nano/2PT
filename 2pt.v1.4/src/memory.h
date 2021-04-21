#ifndef MEMORY_H
#define MEMORY_H

#include <stdlib.h>

using namespace std;

class Memory {
  public:
    void sfree(void *);
    void *smalloc(double, const char *);

    double **create_2d_double_array(int, int, const char *);
    void destroy_2d_double_array(double **);
    int **create_2d_int_array(int, int, const char *);
    void destroy_2d_int_array(double **);

    double *** create_3d_double_array(int, int, int, const char *);
    void destroy_3d_double_array(double ***);

    double **** create_4d_double_array(int, int, int, int, const char *);
    void destroy_4d_double_array(double ****);
    float ****create_4d_float_array(int, int, int, int, const char *);
    void destroy_4d_float_array(float ****);
};

#endif
