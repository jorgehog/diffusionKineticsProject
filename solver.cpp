#include "solver.h"

#include "diffusionscheme.h"

#include <assert.h>
#include <iomanip>

Solver::Solver(DiffusionScheme *scheme, int N, double T) :
    N(N),
    scheme(scheme),
    constants(new Constants(T))
{
    assert(N > 1);
    u.set_size(N);

    CO2.set_size(N);
    c.set_size(N);

    H.set_size(N);
    OH.set_size(N);
    HCO3.set_size(N);
    CO3.set_size(N);

    Ca.set_size(N);

}

void Solver::run(int tSteps)
{
    using namespace std;

    setBoundaryAndInitialConditions();

    vec uNew(N);
    for (int t; t < tSteps; ++t){

        for (int i = 1; i < N-1; ++i) {
            uNew(i) = u(i) + scheme->nextTimeStep(u, constants->DCO2, i);
        }

        iterateKinetics();

        cout << setw(5) << t/(tSteps-1.0)*100 << " %" << endl;

        u(span(1, N-1)) = uNew(span(1, N-1));

    }

}

void Solver::iterateKinetics()
{
    vec dCO2 = -(k1() + k2()*OH)%CO2 +(ka()*H + kM2())%HCO3;

    CO2 += scheme->dt*dCO2;
    c   -= scheme->dt*dCO2;

    chargeBalance();

}


void Solver::setBoundaryAndInitialConditions()
{
    u(1) = u(0);
    u(N-1) = u(N-2); // TRUE?

    CO2.randu();
    c.randu();

    H.randu();
    OH.randu();
    HCO3.randu();
    CO3.randu();

    Ca.randu();
    Cl = as_scalar(randu<vec>(1));
    K = as_scalar(randu<vec>(1));
}



double Solver::charge(double pH)
{

    double Hp = pow(10.0, -pH);
    double OHp = KW()/Hp;
    double HCO3p = Hp*c(j)/(K2()+Hp);
    double CO3p = K2()*HCO3p/Hp;
    double q = 2*Ca(j) + Hp + K - HCO3p - 2*CO3p - OHp - Cl;

    return q;
}

double Solver::bisectRootCharge()
{
    int Nmax = 1E7;
    int N = 0;

    double a = 0;
    double b = 14;
    double c;

    double eps = 1E-16;

    while (N < Nmax){
        c = (a + b)/2;

        if (fabs(charge(c)) < eps){
            std::cout << "yay" << charge(c) << std::endl;
            return c;
        }

        if (charge(a)*charge(c) > 0) {
            a = c;
        } else {
            b = c;
        }

        N++;
    }

    std::cout << "bisection failed " << charge(c) << std::endl;

    return c;

}

void Solver::chargeBalance()
{

    double pH;

    for (j = 0; j < N; ++j) {

        pH = bisectRootCharge();

        H(j) = pow(10.0, -pH);
        OH(j)   = KW()/H(j);             // Eq. (9)
        HCO3(j) = H(j)*c(j)/(K2()+H(j)); // Eq. (10)
        CO3(j)  = K2()*HCO3(j)/H(j);     // Eq. (11)
    }



}

