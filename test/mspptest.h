#ifndef MSPPTEST_H
#define MSPPTEST_H

#include <cmath>
#include "mspp/lpsolver.h"
#include "mspp/distribution.h"

using namespace mspp;

template <typename X>
class guiddistribution: public iddistribution<X>
{
public:
    guiddistribution(const std::vector<X>& items) : fitems(items)
    {
        assert(fitems.size());
    }

protected:
    virtual void atoms_are(std::vector<atom<X>>& a) const
    {
        unsigned int N=fitems.size();
        double p=1.0 / (double) N;
        a.resize(N);
        for(unsigned int i=0; i<N; i++)
            a[i]= {fitems[i],p};
    };
private:
    std::vector<X> fitems;
};

struct omega
{
    double o1;
    double o2;
};

extern void twostagetest(unsigned int N, const lpsolver& cps);
extern void cvartest(double alpha, double lambda, const lpsolver& cps);
extern void almtest(double alpha, double lambda, const lpsolver& cps);


#endif // MSPPTEST_H
