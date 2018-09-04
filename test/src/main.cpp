#include "mspp/msproblem.h"
#include "mspp/random.h"
//#include "mspp/cplex.h"
//#include "test/tstest.h"
#include "test/cvartest.h"
//#include "test/almtest.h"

int main(int argc, char *argv[])
{
    std::ofstream log("sddp.log");
    sys::setlog(log);
//    using O=csvlpsolver<realvar>;
    using O=cplex<realvar>;
//    twostagetest<O>();
//    cvartest<O>(0.05,0.5); // tbd vysledky ale zavisej na hodnotach
//    almtest<O>(0.05,0.5);

    psproblem rnproblem(0.5,0.05);
    gmdddistribution<double> d
           = gddistribution({0, 1.0/3.0, 2.0 / 3.0})
           * gddistribution({0.1+0, 0.1+ 1.0/3.0, 0.1+2.0 / 3.0});

    mpmcvarequivalent<psproblem> mcvproblem(rnproblem);
    sddpmethod<mpmcvarequivalent<psproblem>,gmdddistribution<double>> sm;
    sddpsolution<mpmcvarequivalent<psproblem>> sddpsol(mcvproblem);

    double sddpovalue;
    vector<const gmdddistribution<double>*> dv;
    dv.push_back(&d);
    sm.solve(mcvproblem, dv, vector<double>(0), sddpovalue, sddpsol);

    return 0;
}
