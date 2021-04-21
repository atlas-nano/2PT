#ifndef COMPUTES_H
#define COMPUTES_H

#include <stdlib.h>
#include <fftw3.h>
#include "constant.h"
#include "statistics.h"
#include "model.h"
#include <vector>

class Compute {

public:
    Compute() { }
    virtual ~Compute() {};
    virtual void init(MODEL *) = 0;
    virtual void rd_frame(MODEL *,int) = 0;
    virtual void report(ostream *, MODEL *) = 0;

};

#endif
