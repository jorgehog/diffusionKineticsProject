#ifndef EXPLICITEULER_H
#define EXPLICITEULER_H

#include "diffusionscheme.h"

class ExplicitEuler : public DiffusionScheme
{
public:
    ExplicitEuler(double dt, double dx);

    double nextTimeStep(const arma::mat &u_t, double &D, const int &i) const;

};

#endif // EXPLICITEULER_H
