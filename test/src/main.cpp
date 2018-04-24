#include "mspptest.h"
#include "mspp/cplex.h"

#include "mspp/sddp.h"

int main(int argc, char *argv[])
{
//    csvlpsolver csvs;
    cplexlpsolver csvs;
    twostagetest(3,csvs);
    cvartest(0.05,0.5,csvs);
//    almtest(0.05,0.5,cps);
}

