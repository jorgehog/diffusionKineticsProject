#ifndef DIFFUSIONSCHEME_H
#define DIFFUSIONSCHEME_H

#include <armadillo>

class DiffusionScheme
{
public:
    DiffusionScheme(double dt, double dx);

    virtual double nextTimeStep(const arma::mat &u_t, double &D, const int &i) const = 0;

    const double dt;
    const double dx;

};

#endif // DIFFUSIONSCHEME_H
