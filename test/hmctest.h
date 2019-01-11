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

    arnormalprocessdistribution xipd(1.0,xim,xisd,1.0,T);

    onedcovering c1({-0.2},{-1, 0.6});
    onedcovering c2({-1,0.2},{-1.5,-0.4, 0.7});
    vector<onedcovering> c({c1,c2});

    onedhmcapproximation<arnormalprocessdistribution,onedcovering>
        ha(xipd,c);

}

} // namespace


#endif // HMCTEST_H
