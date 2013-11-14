#include "solver.h"

#include "diffusionscheme.h"

#include <assert.h>
#include <iomanip>
#include <sstream>

Solver::Solver(DiffusionScheme *scheme, int N, double T, double Patm, double Pres) :
    N(N),
    scheme(scheme),
    constants(new Constants(T)),
    Patm(Patm),
    Pres(Pres)
{
    assert(N > 1);

    CO2.set_size(N);
    c.set_size(N);

    H.set_size(N);
    OH.set_size(N);
    HCO3.set_size(N);
    CO3.set_size(N);

    Ca.set_size(N);
    Cl.set_size(N);

    mg_cm2_s.set_size(N);
    accumass.set_size(N);


    pH.set_size(N);

}

void Solver::run(int tSteps)
{
    using namespace std;

    setBoundaryAndInitialConditions();

    for (int t = 0; t < tSteps; ++t){

        diffuse(CO2, constants->DCO2);
        diffuse(Ca,  constants->DIon);
        diffuse(c,   constants->DIon);
        diffuse(Cl,  constants->DIon);

        output("Done diffusion");

        chargeBalance();

        iterateKinetics();

        precipitateDisolve();

        gasExchange();




        dump(t, tSteps);


    }

    superOutput = true;
    output("Simulation finished");

}

void Solver::dump(int i, int tSteps)
{

    using namespace std;

    if (!(i%outTresh == 0)) return;

    mat A(N, 6);

    A.col(0) = c;
    A.col(1) = CO2;
    A.col(2) = Ca;
    A.col(3) = accumass;
    A.col(4) = Cl;
    A.col(5) = pH;

    std::stringstream s;
    s << "/tmp/concOut" << i << ".arma";

    A.save(s.str());

    cout.precision(3);
    cout.setf(ios::fixed);
    cout << setw(5) << i/(tSteps-1.0)*100 << " %" << endl;

}

void Solver::diffuse(vec &u, double D)
{
    vec uNew(N);
    uNew(0) = u(0);
    double fac = D*scheme->dt/(scheme->dx*scheme->dx);
    for (int i = 1; i < N-1; ++i) {
        uNew(i) = u(i) + fac*(u(i+1) -2*u(i) + u(i-1));
    }
    uNew(N-1) = u(N-1) + D*(-u(N-1) + u(N-2));
    u = uNew;
}

void Solver::iterateKinetics()
{
    vec dCO2 = -(k1() + k2()*OH)%CO2 +(ka()*H + kM2())%HCO3;

    CO2 += scheme->dt*dCO2;
    c   -= scheme->dt*dCO2;


    output("Done kinetics");


}

void Solver::output(std::string header)
{
    if (!superOutput) return;

    cout.precision(15);
    cout.setf(ios::fixed);

    int a;
    H.raw_print("H = ");
    OH.raw_print("OH = ");
    HCO3.raw_print("HCO3 = ");
    CO3.raw_print("CO3 = ");
    mg_cm2_s.raw_print("mg_cm2_s = ");
    accumass.raw_print("accumass = ");
    CO2.raw_print("CO2 = ");
    Ca.raw_print("Ca = ");
    c.raw_print("c = ");
    Cl.raw_print("Cl = ");
    std::cout << header << std::endl;
    std::cin >> a;
}



void Solver::setBoundaryAndInitialConditions()
{

    CO2.zeros();

    c.zeros();

    double ClPipe = 0.046;
    Ca.fill(0.5*ClPipe);
    Cl.fill(ClPipe);
    Cl(0) = 0;


    H.zeros();

    OH.zeros();

    OH(0) = ClPipe;

    HCO3.zeros();
    CO3.zeros();
    mg_cm2_s.zeros();
    accumass.zeros();

}



double Solver::charge(double pH)
{

    double Hp = pow(10.0, -pH);
    double OHp = KW()/Hp;
    double HCO3p = Hp*c(j)/(K2()+Hp);
    double CO3p = K2()*HCO3p/Hp;
    double q = 2*Ca(j) + Hp - HCO3p - 2*CO3p - OHp - Cl(j);

    return q;
}

double Solver::chargeDeriv(double pH)
{

    double Hp = pow(10.0, -pH);
    double HCO3p = Hp*c(j)/(K2()+Hp);


    //Setting up partial derivatives.
    double dHdP = -log(10)*Hp;

    double dHCO3dH = c(j)*(K2()*(1-Hp) + Hp)/((K2()+Hp)*(K2() + Hp));

    double dCO3dH = K2()*(dHCO3dH*Hp - HCO3p)/(Hp*Hp);

    double dOHdH = -KW()/(Hp*Hp);

//    return 2*pH;

    return (1 - dHCO3dH - 2*dCO3dH - dOHdH)*dHdP;

}

double Solver::fixPoint(double start)
{

    double eps = 1E-15;
    int NMAX = (int)1E5;

    int N = 0;


    while ((fabs(charge(start)) > eps) && (N < NMAX)) {

        start = charge(start) - start;
        N++;
    }

    return start;



}

double Solver::bisectRootCharge()
{
    int Nmax = 1E3;
    int N = 0;

    double a = 0;
    double b = 14;
    double C;

    double eps = 0.01;


    double delta = 0.001;

    while(charge(a)*charge(b) > 0) {
        b -= delta;

        if (b < a) {
            std::cout << "EXPRESSION HAS NO ROOTS!!" << std::endl;
            return 0;
        }
    }


    double CPrev = 10000;
    while (N < Nmax){
        C = (a + b)/2;

        if (fabs(CPrev - C) < eps){
            return C;
        }

        if (charge(a)*charge(C) > 0) {
            a = C;
        } else {
            b = C;
        }

        CPrev = C;

        N++;
    }

    return C;

}

double Solver::newtonMethodRootCharge(double start)
{
    double a0 = 0;
    double a1 = start;

    double eps = 1E-8;
    int NMAX = (int)1E5;

    int N = 0;

    while ((fabs(charge(a1)) > eps) && (N < NMAX)) {

        a0 = a1;

        a1 = a0 - charge(a0)/chargeDeriv(a0);

        N++;
    }



    return a1;
}

double Solver::secantMethod(double start)
{

    double a2 = 0;
    double a1 = start - start/10;
    double a0 = start + start/10;

    double eps = 1E-15;
    int NMAX;
    NMAX = (int)1000;

    int N = 0;


    while ((fabs(charge(a1)) > eps) && (N < NMAX)) {

        a2 = a0 - charge(a0)*(a1-a0)/(charge(a1)-charge(a0));

        a0 = a1;
        a1 = a2;

        N++;
    }

    return a1;
}

void Solver::chargeBalance()
{

    double pH;

    //member j is needed by charge function.
    for (j = 0; j < N; ++j) {

        pH = bisectRootCharge();
        pH = secantMethod(pH);

        this->pH(j) = pH;
        H(j) = pow(10.0, -pH);
        OH(j)   = KW()/H(j);             // Eq. (9)
        HCO3(j) = H(j)*c(j)/(K2()+H(j)); // Eq. (10)
        CO3(j)  = K2()*HCO3(j)/H(j);     // Eq. (11)
    }



    output("Done chargeBalance");
}

void Solver::gasExchange()
{

    double k = 3E-4;

    double FCO2 = k*(KH()*Pres - CO2(N-1));
    CO2(N-1) += scheme->dt*FCO2/scheme->dx;


    output("Done gas");



}

void Solver::precipitateDisolve()
{

    double dM, F, delta;
    for (int j = 1; j < N; ++j) {

        F = 1000*(9.5e-11 - (7.14e-3)*Ca(j)*CO3(j)); // Compton-Pritchard

        dM = -F*100;

        if (accumass(j) + dM < 0) {
//            std::cout << "unable to dissolve CaCO3. Nothing present" << std::endl;
            F = 0;
        }

        delta = scheme->dt*F/scheme->dx;


        Ca(j) += delta;
        c(j)  += delta;

        mg_cm2_s(j) = -F*100; // Molecular weight of CaCO3

    }

    accumass += mg_cm2_s;

    output("Done dissolve");

}

