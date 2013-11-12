#include "expliciteuler.h"

#include <assert.h>


ExplicitEuler::ExplicitEuler(double dt, double dx) :
    DiffusionScheme(dt, dx)
{
    assert(dt < dx && "invalid choise of time step. Solution will not converge.");
}

double ExplicitEuler::nextTimeStep(const arma::mat &u_t, double & D, const int & i) const
{
    return D*dt/(dx*dx)*(
                             u_t(i + 1)
                         - 2*u_t(i)
                         +   u_t(i - 1)

                        );
}


