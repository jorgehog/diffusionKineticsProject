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

//        superOutput = true;
        gasExchange();
//        superOutput = false;




        dump(t);

        cout << setw(5) << t/(tSteps-1.0)*100 << " %" << endl;

    }
    cout << c << endl;
    cout << CO2 << endl;
    cout << Ca<< endl;
    cout << H << endl;
    cout << OH << endl;
    cout << mg_cm2_s << endl;

}

void Solver::dump(int i)
{
//    c, co2, ca, H, OH, mass
    mat A(N, 8);

    A.col(0) = c;
    A.col(1) = CO2;
    A.col(2) = Ca;
    A.col(3) = H;
    A.col(4) = OH;
    A.col(5) = mg_cm2_s;
    A.col(6) = Cl;
    A.col(7) = pH;

    std::stringstream s;
    s << "/home/jorgehog/scratch/concOut" << i << ".arma";

    A.save(s.str());

}

void Solver::diffuse(vec &u, double D)
{
    vec uNew(N);
    uNew(0) = u(0);
    double fac = D*scheme->dt/(scheme->dx*scheme->dx);
    for (int i = 1; i < N-1; ++i) {
        uNew(i) = u(i) + fac*(u(i+1) -2*u(i) + u(i-1));
    }
    using namespace std;
//    std::cout << D << endl << scheme->dt << endl << scheme->dx << "\n-----"<<std::endl;
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
    CO2.raw_print("CO2 = ");
    Ca.raw_print("Ca = ");
    c.raw_print("c = ");
    Cl.raw_print("Cl = ");
    std::cout << header << std::endl;
    std::cin >> a;
}



void Solver::setBoundaryAndInitialConditions()
{

    Constants c2(373.15);
    CO2.fill(c2.KH*Patm); //Henry's law
    CO2(N-1) = KH()*Patm;
    CO2(0) = 0;
//    CO2.zeros();
    c.zeros();

    double ClPipe = 0.046;
    Ca.fill(0.5*ClPipe);
    Cl.fill(ClPipe);
    Cl(0) = 0;


    H.zeros();

    OH.zeros();

    OH(0) = ClPipe;
//    H(0) = KW()/OH(0);

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
    double q = 2*Ca(j) + Hp - HCO3p - 2*CO3p - OHp - Cl(j);

//    return pH*pH - 2;

    return q;
}

double Solver::chargeDeriv(double pH)
{

//    return (charge(pH + 0.001) - charge(pH - 0.001))/0.002;
    double Hp = pow(10.0, -pH);
    double HCO3p = Hp*c(j)/(K2()+Hp);

    //  double OHp = KW()/Hp;
    //    double CO3p = K2()*HCO3p/Hp;
    //    double dHCO3dP = dHCO3dH*dHdP;
    //    double dCO3dP = dCO3dH*dHdP;
    //    double dOHdP = dOHdH*dHdP;


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

    double eps = 0;


    double delta = 0.001;

//    double bX = b;
    while(charge(a)*charge(b) > 0) {
        b -= delta;

        if (b < a) {
            std::cout << "EXPRESSION HAS NO ROOTS!!" << std::endl;
            return 0;
        }
    }

//    if (b != bX) {
//        std::cout << "YEEHAA" << std::endl;
//        std::cout << bX << "  " << b << std::endl;
//    }
//    std::cout << charge(a)*charge(b) << std::endl;


    while (N < Nmax){
        C = (a + b)/2;

        if (fabs(charge(C)) < eps){
            return C;
        }

        if (charge(a)*charge(C) > 0) {
            a = C;
        } else {
            b = C;
        }

        N++;
    }

    return C;

}

double Solver::newtonMethodRootCharge(double start)
{
    double a0 = 0;
    double a1 = start;

    double eps = 1E-15;
    int NMAX = (int)1E5;

    int N = 0;

//    double deriv = (charge(a1 + 0.001) - charge(a1 - 0.001))/0.002;
//    std::cout << deriv << "  " << chargeDeriv(a1) << std::endl;
//    exit(1);

    while ((fabs(charge(a1)) > eps) && (N < NMAX)) {

        a0 = a1;

        a1 = a0 - charge(a0)/chargeDeriv(a0);

        N++;
//        std::cout << "new solution "<< a1<< "  " << fabs(charge(a1)) <<  std::endl;
    }

//    std::cout << N << std::endl;

    return a1;
}

double Solver::secantMethod(double start)
{

    double a2 = 0;
    double a1 = start - start/10;
    double a0 = start + start/10;

    double eps = 1E-15;
    int NMAX = (int)1E5;

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

//        pH = fixPoint(7);
        pH = bisectRootCharge();
//        std::cout << pH << "  " << charge(pH) << std::endl;
//        pH = secantMethod(pH);
//        std::cout << pH << "  " << charge(pH) << std::endl;
        pH = newtonMethodRootCharge(pH);
//        std::cout << pH << "  " << charge(pH) << std::endl;
//        std::cout << "---------------" << std::endl;

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

    double F, delta;
    for (int j = 1; j < N; ++j) {

        F = 1000*(9.5e-11 - (7.14e-3)*Ca(j)*CO3(j)); // Compton-Pritchard

        delta = scheme->dt*F/scheme->dx;
        Ca(j) += delta;
        c(j)  += delta;

        mg_cm2_s(j) = -F*100; // Molecular weight of CaCO3

    }


    output("Done dissolve");


}

