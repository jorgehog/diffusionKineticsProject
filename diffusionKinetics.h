#ifndef DIFFUSIONKINETICS_H
#define DIFFUSIONKINETICS_H

#include "solver.h"
#include "Schemes.h"

struct configParams {
    SCEMETYPE type;
    int N;
    int nSteps;
    double H;
    double D;
    double T;
    double dt;
};

#endif // DIFFUSIONKINETICS_H
