#include "diffusionKinetics.h"

#include <iostream>
#include <libconfig.h++>
#include <DCViz.h>
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

    DCViz viz("/tmp/concOut0.arma");
    viz.launch(true, 0, 16, 14);

    parseConfigFile(argv[1], cfgParams);


    //Calculating the cellHeigth
    double h = cfgParams.H/(cfgParams.N-1);

    //Selecting scheme based on type given in config
    DiffusionScheme *scheme = selectScheme(cfgParams.type, h, cfgParams.dt);


    //Create and run the solver
    Solver solver(scheme, cfgParams.N, cfgParams.T, cfgParams.Patm, cfgParams.Pres);
    solver.run(cfgParams.nSteps);

    return 0;
}




void parseConfigFile(const char *filename, configParams &cfgParams){

    Config cfg;

    cfg.readFile(filename);

    cfgParams.type = static_cast<SCEMETYPE>((int)cfg.lookup("scheme.type"));

    cfgParams.H  = cfg.lookup("environment.H");
    cfgParams.T  = (double)cfg.lookup("environment.T") + 273.15;

    cfgParams.dt = cfg.lookup("solver.dt");
    cfgParams.N  = cfg.lookup("solver.N");
    cfgParams.nSteps = cfg.lookup("solver.steps");

    cfgParams.Patm = cfg.lookup("environment.Patm");
    cfgParams.Pres = ((double)cfg.lookup("environment.PresOverPatm"))*cfgParams.Patm;

}

DiffusionScheme* selectScheme(SCEMETYPE type, double dx, double dt) {

    switch (type) {

    case EXPLICIT_EULER:
        return new ExplicitEuler(dt, dx);
        break;

    default:
        return NULL;
        break;
    }

}
