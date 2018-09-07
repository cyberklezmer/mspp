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
//    std::ofstream elog("esddp.log");
//    sys::seterr(elog);

    sys::seed(0);
//    using O=csvlpsolver<realvar>;
    using O=cplex<realvar>;
//    twostagetest<O>();
//    cvartest<O>(0.05,0.5); // tbd vysledky ale zavisej na hodnotach
//    almtest<O>(0.05,0.5);

    psproblem rnproblem(0.5,0.05);
    gmdddistribution<double> d
           = gddistribution({0, 1.0/3.0, 2.0 / 3.0})
           * gddistribution({0.1+0, 0.1+ 1.0/3.0, 0.1+2.0 / 3.0});

    vector<const gmdddistribution<double>*> dv;
    dv.push_back(&d);

    sddpmethod<psproblem,gmdddistribution<double>> dsm;
    sddpsolution<psproblem> dsddpsol(rnproblem);

    dsm.solve(rnproblem, dv, vector<double>(0), dsddpsol);

    cout << "lb: " << dsddpsol.lb() << " ubm: " << dsddpsol.ubm()
         <<     " ubb: " << dsddpsol.ubb() << endl;
    for(unsigned int i=0; i<dsddpsol.firststage().size(); i++)
        cout << "x" << i << "=" << dsddpsol.firststage()[i] << endl;

    mpmcvarequivalent<psproblem> mcvproblem(rnproblem);
    sddpmethod<mpmcvarequivalent<psproblem>,gmdddistribution<double>> sm;
    sddpsolution<mpmcvarequivalent<psproblem>> sddpsol(mcvproblem);

    for(unsigned int i = 0; i<0; i++)
    {
        sm.solve(mcvproblem, dv, vector<double>(0), sddpsol);

        cout << "lb: " << sddpsol.lb() << " ubm: " << sddpsol.ubm()
             <<     " ubb: " << sddpsol.ubb() << endl;
        for(unsigned int i=0; i<sddpsol.firststage().size(); i++)
            cout << "x" << i << "=" << sddpsol.firststage()[i] << endl;
    }
    return 0;
}
