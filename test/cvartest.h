#ifndef CVARTEST_H
#define CVARTEST_H

#include "mspptest.h"
#include "mspp/sddp.h"
#include "mspp/mpcproblem.h"
#include "mspp/de.h"

class cvarproblem: public msproblem<mpmcvar,
        linearfunction,
        linearmsconstraint,realvar,
        rvector<double>,allx,lastxi<vector<double>>>
{

public:
    cvarproblem(double lambda, double alpha) :
           msproblem<mpmcvar,
           linearfunction,
           linearmsconstraint,realvar,vector<double>,allx,lastxi<vector<double>>>
          (msproblemstructure({2,2}), mpmcvar(lambda,alpha))
    {}

    virtual void f_is(
                    unsigned int k,
                    const Z_t::C_t& zeta,
                    linearfunction& f) const
    {
      vector<double> xik=zeta;
      if(k)
         f = linearfunction({-xik[0],-xik[1]});
      else
         f = linearfunction({0,0});
    }



    virtual void x_is(
            unsigned int k,
            const Z_t::C_t& zeta,
            ranges<realvar>& r,
            msconstraints<linearmsconstraint>& g
            ) const
    {
        assert(k<=1);

        if(k == 0)
        {
            g.add(linearmsconstraint({1.0,1.0},
                          constraint::eq, 1.0));
            r[0].setpositive();
            r[1].setpositive();
        }
        else
        {
            r[0].setpositive();
            r[1].setpositive();
            g.add(linearmsconstraint({1.0,0,-1.0,0},
                          constraint::eq, 0));
            g.add(linearmsconstraint({0,1.0,0,-1.0},
                          constraint::eq, 0));
        }
    }
    virtual double minf_is(unsigned int) const
    {
        return -1e10;
    }
    virtual double maxf_is(unsigned int) const
    {
        return 1e10;
    }
};

template <typename O>
void cvartest(bool sddp=false, bool indirect=false)
{
    gmdddistribution<double> d
           = ldistribution({0, 1.0/3.0, 2.0 / 3.0})
           * ldistribution({0.1+0, 0.1+ 1.0/3.0, 0.1+2.0 / 3.0});

    using dist = processdistribution<gmdddistribution<double>>;
    dist pd(d,1);

    using dst = distrscenariotree<dist>;
    dst s(pd);

    cvarproblem p(0.5,0.05);

    mpmcvarequivalent<cvarproblem> ep(p);

    const double tol = 1e-5;

    double obj = -0.266666666666667;

    double sol[]={-0.1,0,1}; // u0, x0_1,x_02
    double thetas[]={0,-0.166666666666666,-0.333333333333335,
                      0,-0.166666666666669,-0.333333333333333,
                      0,-0.166666666666667,-0.333333333333332};

    // brought from testproblems.xlsx
    // x0_1,x_02,u0,theta_1,theta_2,theta_3,theta_4,theta_5,theta_6,theta_7,theta_8,theta_9,,,

    if(!sddp && !indirect)
    {
        std::cout << "CVaRtest direct by DE... ";

        stsolution<cvarproblem,dst> x(p,s);
        demethod::solve<cvarproblem,dst,O>(p, s, x);

        if(fabs(x.obj()-obj) > tol)
        {
            std::cerr << "opt="
             << obj << " expected, " << x.obj() << " achieved."
             << std::endl;

            for(unsigned int i=0; i<2+9*2; i++)
                cerr << x.x(i) << endl;
            throw;
        }

        if(fabs(x.x(0)) > tol)
        {
            std::cerr << "x[" << 0 << "]="
                 << 0 << " expected, " << x.x(0)
                 << " achieved." << std::endl;
            throw;
        }

        std::cout << "passed" << endl;
    }

    if(!sddp && indirect)
    {
        std::cout << "CVaRtest indirect by DE...";

        stsolution<mpmcvarequivalent<cvarproblem>,dst> x(ep,s);

        demethod::solve<mpmcvarequivalent<cvarproblem>,dst,O>(ep,s,x);

        if(fabs(x.obj()-obj) > tol)
        {
            std::cerr << "opt="
             << obj << " expected, " << x.obj() << " achieved."
             << std::endl;

            for(unsigned int i=0; i<3+9*3; i++)
                cerr << x.x(i) << endl;
            throw;
        }


        unsigned int k = sizeof(sol)/sizeof(sol[0]);
        for(unsigned int i=0; i<k; i++)
        {
            if(fabs(x.x(i)-sol[i]) > tol)
            {
                std::cerr << "x[" << i << "]="
                     << sol[i] << " expected, " << x.x(i)
                     << " achieved." << std::endl;
                throw;
            }
        }

        unsigned int l = sizeof(thetas)/sizeof(thetas[0]);
        unsigned int thetafactor = 3;
        for(unsigned int i=0; i<l; i++)
        {
            if(fabs(x.x(k+thetafactor*i)-thetas[i]) > tol)
            {
                std::cerr << "theta[" << i << "]="
                 << thetas[i] << " expected, " << x.x(k+2+3*i)
                 << " achieved." << std::endl;

                throw;
            }
        }
        std::cout <<  "Passed." << std::endl;
    }

    if(sddp && !indirect)
    {
        std::cout << "CVaRtest direct by SDDP... ";

        sddpsolution<cvarproblem> x(p);

        sddpmethod::solve<cvarproblem,dist>(p, pd, x);

        if(fabs(x.lb()-obj) > 0.1 || fabs(x.ubb()-obj) > 0.1)
        {
            std::cerr << "obj="
               << obj << " expected, " << x.lb()
               << "-" << x.ubb() << " achieved." << std::endl;

            for(unsigned int i=0; i<x.firststage().size(); i++)
               cerr << "x" << i << "=" << x.firststage()[i] << endl;

            throw;
        }

        std::cout << "passed" << endl;
    }

    if(sddp && indirect)
    {
        std::cout << "CVaRtest indirect by SDDP... ";

        sddpsolution<mpmcvarequivalent<cvarproblem>> x(ep);

        sddpmethod::solve<mpmcvarequivalent<cvarproblem>,dist>(ep, pd, x);

        if(fabs(x.lb()-obj) > 0.01 || fabs(x.ubb()-obj) > 0.01)
        {
            std::cerr << "obj="
               << obj << " expected, " << x.lb()
               << "-" << x.ubb() << " achieved." << std::endl;

            for(unsigned int i=0; i<x.firststage().size(); i++)
               cerr << "x" << i << "=" << x.firststage()[i] << endl;

            throw;
        }

        std::cout << "passed" << endl;
    }
}

#endif // CVARTEST_H
