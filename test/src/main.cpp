  #include <fstream>
#include "mspp/random.h"
#include "mspp/process.h"
#include "test/mspptest.h"

//#include "mspp/msproblem.h"
//#include "mspp/cplex.h"
//#include "test/tstest.h"
//#include "test/cvartest.h"
//#include "test/almtest.h"
#include "test/hmctest.h"
//#include "mspp/cplex.h"
#include <mcheck.h>

using namespace mspp;


void compiletest(unsigned int T)
{
    assert(T>=1);
    assert(T<=4);
    const double tol = 1e-5;

    using etapd = iidprocessdistribution<normaldistribution>;

    double etam = -0.5;
    double etasd = 0.2;

    etapd etap(0.0,normaldistribution(etam,etasd),T);

    double xim = -0.2;
    double xisd = 1;

    arnormalprocessdistribution xipd(1.0,xim,xisd,1.0,T);

    onedcovering c1({-0.2},{-1, 0.6});
    onedcovering c2({-1,0.2},{-1.5,-0.4, 0.7});
    vector<onedcovering> c({c1,c2});

    onedhmcapproximation<arnormalprocessdistribution,onedcovering>
        ha(xipd,c);
}



int main(int, char **)
{
    int res=mcheck(nullptr);
    if(res)
        cerr << "error " << res << " starting mcheck" << endl;

    std::ofstream log("sddp.log");
    sys::setlog(log);

    sys::seed(0);
    using O=csvlpsolver<realvar>;
//    using O=cplex<realvar>;
//    twostagetest<O>();
//    cvartest<O>(false,false);
//    cvartest<O>(false,true);
//    cvartest<O>(true,false);
//    cvartest<O>(true,true);
//    almtest<O>(3,1); // dává lb< ub
//    alm1test<O>(true);
//    alm1test<O>(false);


// just to check compilation


    compiletest(2);

    return 0;
}
