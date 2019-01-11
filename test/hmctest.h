#ifndef HMCTEST_H
#define HMCTEST_H

#include "almproblem.h"
#include "mspp/sddp.h"
#include "mspp/de.h"
#include "mspp/cdists.h"
#include "mspp/hmcapprox.h"

namespace mspp
{

void hmctest(unsigned int T)
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

    normaldistribution nd(xim,xisd);
    using arnormaldistribution = ardistribution<normaldistribution>;

    arnormaldistribution ar(nd,1);

    using xipd = markovprocessdistribution<arnormaldistribution>;

    xipd xip(1.0,ar,T);
}

} // namespace


#endif // HMCTEST_H
