#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <math.h>
#include <iostream>

class Constants {
public:
    double DCO2;
    double DIon;

    double k1;
    double k2;
    double ka;
    double kM2;

    double KW;
    double K2;

    Constants(double T) {
        DCO2 = 1E-5*(0.56 + 0.058*T);
        DIon = 0.7*DCO2;

        k1 = pow(10.0, 329.850-110.54*log10(T)-17265.4/T);
        k2 = pow(10.0, 13.635-2895/T);

        double k_minus1 = pow(10.0, 13.558-3617.1/T);
        double K5 = 1.707e-4;
        ka = k_minus1 / K5;
        kM2 = pow(10.0, 14.09-5308/T);

        KW = pow(10.0, 22.801 - 4787.3/T - 0.010365*T - 7.1321*log10(T));
        K2 = pow(10.0, -107.8871 + 5151.79/T - 0.03252849*T + 38.9256*log10(T) -
        563713.9/(T*T));

        std::cout << DCO2 << std::endl;
        std::cout << DIon << std::endl;
        std::cout << k1   << std::endl;
        std::cout << k2   << std::endl;
        std::cout << ka   << std::endl;
        std::cout << kM2  << std::endl;
        std::cout << KW   << std::endl;
        std::cout << K2   << std::endl;

    }

};

#endif // CONSTANTS_H

