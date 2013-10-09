#ifndef SOLVER_H
#define SOLVER_H

#include "constants.h"

#include <armadillo>
using namespace arma;

#include <vector>

class DiffusionScheme;

class Solver
{
public:

    const int N;

    double Patm;
    double Pres;

    Solver(DiffusionScheme *scheme, int N, double T, double Patm, double Pres);

    DiffusionScheme* scheme;

    Constants * constants;

    vec u;

    vec CO2;
    vec c;

    vec H;
    vec OH;
    vec HCO3;
    vec CO3;

    vec Ca;
    double Cl;
    double K;

    vec mg_cm2_s;

    void run(int tSteps);

    void diffuse(vec & u, double D);

    void setBoundaryAndInitialConditions();

    void iterateKinetics();


    double charge(double pH);
    double bisectRootCharge();
    void chargeBalance();
    int j;

    void gasExchange();

    void precipitateDisolve();


    double k1() {
        return constants->k1;
    }

    double k2() {
        return constants->k2;
    }

    double ka() {
        return constants->ka;
    }

    double kM2() {
        return constants->kM2;
    }

    double KW() {
        return constants->KW;
    }

    double K2() {
        return constants->K2;
    }

    double KH() {
        return constants->KH;
    }

};

#endif // SOLVER_H
