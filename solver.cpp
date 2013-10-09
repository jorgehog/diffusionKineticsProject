#include "solver.h"

#include "diffusionscheme.h"

#include <assert.h>
#include <iomanip>

Solver::Solver(DiffusionScheme *scheme, int N, double T, double Patm, double Pres) :
    N(N),
    scheme(scheme),
    constants(new Constants(T)),
    Patm(Patm),
    Pres(Pres)
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

    mg_cm2_s.set_size(N);

}

void Solver::run(int tSteps)
{
    using namespace std;

    setBoundaryAndInitialConditions();

    for (int t; t < tSteps; ++t){

        diffuse(CO2, constants->DCO2);
        diffuse(Ca,  constants->DIon);
        diffuse(c,   constants->DIon);

        chargeBalance();

        iterateKinetics();

        precipitateDisolve();

        gasExchange();

        cout << setw(5) << t/(tSteps-1.0)*100 << " %" << endl;

    }

}

void Solver::diffuse(vec &u, double D)
{
    vec uNew(N-2);
    for (int i = 1; i < N-1; ++i) {
        uNew(i-1) = u(i) + scheme->nextTimeStep(u, D, i);
    }
    u(span(1, N-2)) = uNew;
}

void Solver::iterateKinetics()
{
    vec dCO2 = -(k1() + k2()*OH)%CO2 +(ka()*H + kM2())%HCO3;

    CO2 += scheme->dt*dCO2;
    c   -= scheme->dt*dCO2;

}



void Solver::setBoundaryAndInitialConditions()
{
    u(1) = u(0);
    u(N-1) = u(N-2); // TRUE?

    CO2.fill(KH()*Patm); //Henry's law
    c.zeros();

    Cl = 0.04;
    Ca.fill(0.5*Cl);
    K = 0.02;

    H.zeros();
    OH.zeros();
    HCO3.zeros();
    CO3.zeros();
    mg_cm2_s.zeros();

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

    double eps = 1E-3;

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

    //member j is needed by charge function.
    for (j = 0; j < N; ++j) {

        pH = bisectRootCharge();

        H(j) = pow(10.0, -pH);
        OH(j)   = KW()/H(j);             // Eq. (9)
        HCO3(j) = H(j)*c(j)/(K2()+H(j)); // Eq. (10)
        CO3(j)  = K2()*HCO3(j)/H(j);     // Eq. (11)
    }



}

void Solver::gasExchange()
{

    double k = 3E-4;

    double FCO2 = k*(KH()*Pres - CO2(N-1));
    CO2(N-1) = CO2(N-1) + scheme->dt*FCO2/scheme->dx;


}

void Solver::precipitateDisolve()
{

    double F, delta;
    for (int j = 1; j < N; ++j) {

        F = 1000*(9.5e-11 - (7.14e-3)*Ca(j)*CO3(j)); // Compton-Pritchard

        delta = scheme->dt*F/scheme->dx;
        Ca(j) = Ca(j) + delta;
        c(j) = c(j) + delta;

        mg_cm2_s(j) = -F*100; // Molecular weight of CaCO3

    }



}

