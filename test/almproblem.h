#ifndef ALMPROBLEM_H
#define ALMPROBLEM_H

#include "mspp/msproblem.h"

namespace mspp

{
template <bool CONSTL = true>
class almproblem: public msproblem<mpmcvar, linearfunction,
        linearmsconstraint,vector<double>,realvar,lastx>
{
    static std::vector<unsigned int> makeps(unsigned int t)
    {
        std::vector<unsigned int> res = {1};
        for(unsigned int i=1; i<=t; i++)
            res.push_back(2);
        return res;
    }
public:
    almproblem(double lambda, double alpha, unsigned int T) :
        msproblem<mpmcvar, linearfunction,
                linearmsconstraint,vector<double>,realvar,lastx>
        (makeps(T),mpmcvar(lambda,alpha))
    {}

    virtual void f_is(unsigned int k,
                      const vector<double>& zeta,
                      linearfunction& f
                      ) const
    {
        if(k)
            f.setc(0,zeta[0]);
        else
            f.setc(0,1);
    }

    virtual  void x_is(
            unsigned int k,
            const vector<double>& zeta,
            ranges<realvar>& xs,
            msconstraints<linearmsconstraint>& g
            ) const
    {
        if(k)
            xs[1].setlimits(0,1);
        else
            xs[0].setlimits(0,1);
        if(!k)
            return;
        double theta;
        if constexpr(CONSTL)
            theta = k==this->T() ? 1 : 0;
        else
            theta = zeta[1];
        if(k==1)
            g.add(linearmsconstraint({1.0,1.0,-1.0},constraint::eq, theta));
        else if(k==this->T())
        {
            g.add(linearmsconstraint({0.0,1.0,1.0,-1.0},constraint::eq, 0));
            xs[1].setlimits(theta,theta);
        }
        else
          g.add(linearmsconstraint({0.0,1.0,1.0,-1.0},constraint::eq, theta));
    }
    double minf_is(unsigned int) const
    {
        return -1000;
    }
    double maxf_is(unsigned int) const
    {
        return 1000;
    }
};

} // namespace
#endif // ALMPROBLEM_H
