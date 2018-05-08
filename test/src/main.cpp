#define MSPP_SYS

#include "mspp/cplex.h"
#include "test/tstest.h"
#include "test/cvartest.h"
#include "test/almtest.h"

#include "mspp/sddp.h"

sys* sys::fself=0;
sys gs;

int main(int argc, char *argv[])
{
    std::ofstream log("mspp.log");
    sys::setlog(log);
//    using O=csvlpsolver<realvar>;
    using O=cplex<realvar>;
    twostagetest<O>();
    cvartest<O>(0.05,0.5,false);
    cvartest<O>(0.05,0.5,true);
    almtest<O>(0.05,0.5);
}

