#include "diffusionKinetics.h"

#include <iostream>
#include <libconfig.h++>
#include <assert.h>

using namespace libconfig;
using namespace std;

void parseConfigFile(const char *filename, configParams &cfgParams);
DiffusionScheme *selectScheme(SCEMETYPE type, double dx, double dt);

int main(int argc, char ** argv)
{
    //LOADING CONFIG
    assert(argc > 1 && "Config file not supplied.");
    configParams cfgParams;
    parseConfigFile(argv[1], cfgParams);


    //Calculating the cellHeigth
    double h = cfgParams.H/(cfgParams.N-1);

    //Selecting scheme based on type given in config
    DiffusionScheme *scheme = selectScheme(cfgParams.type, h, cfgParams.dt);


    //Create and run the solver
    Solver solver(scheme, cfgParams.N, cfgParams.T);
    solver.run(cfgParams.nSteps);

    return 0;
}




void parseConfigFile(const char *filename, configParams &cfgParams){

    Config cfg;

    cfg.readFile(filename);

    cfgParams.type = static_cast<SCEMETYPE>((int)cfg.lookup("scheme.type"));

    cfgParams.H  = cfg.lookup("environment.H");
    cfgParams.T  = cfg.lookup("environment.T");

    cfgParams.dt = cfg.lookup("solver.dt");
    cfgParams.N  = cfg.lookup("solver.N");
    cfgParams.nSteps = cfg.lookup("solver.steps");
}

DiffusionScheme* selectScheme(SCEMETYPE type, double dx, double dt) {

    switch (type) {

    case EXPLICIT_EULER:
        return new ExplicitEuler(dx, dt);
        break;

    default:
        return NULL;
        break;
    }

}