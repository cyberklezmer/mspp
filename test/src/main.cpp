#include <fstream>
#include "mspp/random.h"
#include "mspp/process.h"
//#include "mspp/de.h"
//#include "mspp/msproblem.h"
//#include "mspp/cplex.h"
//#include "test/tstest.h"
//#include "test/cvartest.h"
//#include "test/almtest.h"
#include <mcheck.h>

using namespace mspp;

int main(int, char **)
{
    int res=mcheck(nullptr);
    if(res)
        cerr << "error " << res << " starting mcheck" << endl;

    std::ofstream log("sddp.log");
    sys::setlog(log);

    sys::seed(0);
//    using O=csvlpsolver<realvar>;
//    using O=cplex<realvar>;
//    twostagetest<O>();
//    cvartest<O>(false,false);
//    cvartest<O>(false,true);
//    cvartest<O>(true,false);
//    cvartest<O>(true,true);
//    almtest<O>();
    return 0;
}
