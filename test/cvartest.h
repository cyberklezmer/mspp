#ifndef CVARTEST_H
#define CVARTEST_H

#include "mspptest.h"
#include "mspp/sddp.h"
#include "mspp/mpcproblem.h"
#include "mspp/de.h"

class psproblem: public msproblem< expectation,
        linearfunction,
        linearmsconstraint,realvar,
        rvector<double>,allx,lastxi<vector<double>>>
{

public:
    psproblem(double lambda, double alpha) :
           msproblem<expectation,
           linearfunction,
           linearmsconstraint,realvar,vector<double>,allx,lastxi<vector<double>>>
          (msproblemstructure({2,2})
//,mpmcvar(lambda,alpha)
           )
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
void cvartest()
{
    gmdddistribution<double> d
           = gddistribution({0, 1.0/3.0, 2.0 / 3.0})
           * gddistribution({0.1+0, 0.1+ 1.0/3.0, 0.1+2.0 / 3.0});

    using pdist = processdistribution<gmdddistribution<double>>;
    pdist pd(d,1);
    using myst = distrscenariotree<pdist>;
    myst stree(pd);

    psproblem rnproblem(0.5,0.05);

    const double tol = 1e-5;

    double cvarobj = -0.266666666666667;

    std::cout << "CVaRtest direct ";

    demethod<psproblem,myst,O> m;
    stsolution<psproblem,myst> s(rnproblem,stree);

    double o;

/*    m.solve(rnproblem, stree, o, s);

    if(fabs(o-cvarobj) > tol)
    {
        std::cerr << "cvartest direct: opt="
             << cvarobj << " expected, " << o << " achieved." << std::endl;

        for(unsigned int i=0; i<2+9*2; i++)
        {
            cerr << s.x(i) << endl;
        }
        throw;
    }

    if(fabs(s.x(0)) > tol)
    {
        std::cerr << "cvartest direct: x[" << 0 << "]="
             << 0 << " expected, " << s.x(0)
             << " achieved." << std::endl;
        throw;
    }

    if(fabs(s.x(1)-1) > tol)
    {
        std::cerr << "cvartest direct: x[" << 1 << "]="
             << 1 << " expected, " << s.x(1)
             << " achieved." << std::endl;
        throw;
    }

    std::cout << "Passed."  << std::endl;


    std::cout << "CVaRtest indirect by DE ";

    mpmcvarequivalent<psproblem> mcvproblem(rnproblem);

    demethod<mpmcvarequivalent<psproblem>,myst,O> b;

    stsolution<mpmcvarequivalent<psproblem>,myst>
            sol(mcvproblem,stree);
// jaktoze tu povolil sol(raproblem, stree)?
    double ovalue;

    b.solve(mcvproblem,stree, ovalue, sol);


    // brought from testproblems.xlsx
    // x0_1,x_02,u0,theta_1,theta_2,theta_3,theta_4,theta_5,theta_6,theta_7,theta_8,theta_9,,,

    double cvarsol[]={-0.1,0,1}; // u0, x0_1,x_02
    double cvarthetas[]={0,-0.166666666666666,-0.333333333333335,
                      0,-0.166666666666669,-0.333333333333333,
                      0,-0.166666666666667,-0.333333333333332};


    if(fabs(ovalue-cvarobj) > tol)
    {
        std::cerr << "cvartest: opt="
             << cvarobj << " expected, " << ovalue << " achieved." << std::endl;

        for(unsigned int i=0; i<3+9*3; i++)
        {
            cerr << s.x(i) << endl;
        }

        throw;
    }


    unsigned int k = sizeof(cvarsol)/sizeof(cvarsol[0]);
    for(unsigned int i=0; i<k; i++)
    {
        if(fabs(sol.x(i)-cvarsol[i]) > tol)
        {
            std::cerr << "cvartest: x[" << i << "]="
                 << cvarsol[i] << " expected, " << sol.x(i)
                 << " achieved." << std::endl;
            throw;
        }
    }

    unsigned int l = sizeof(cvarthetas)/sizeof(cvarthetas[0]);

    unsigned int thetafactor = 3;
    for(unsigned int i=0; i<l; i++)
    {
        if(fabs(sol.x(k+thetafactor*i)-cvarthetas[i]) > tol)
        {
            std::cerr << "cvartest: theta[" << i << "]="
                 << cvarthetas[i] << " expected, " << sol.x(k+2+3*i)
                 << " achieved." << std::endl;

            throw;
        }
    }
    std::cout <<  "Passed." << std::endl;

*/
    std::cout << "CVaRtest direct by SDDP ";

    sddpmethod<psproblem,pdist> dsm;
    sddpsolution<psproblem> dsddpsol(rnproblem);

    dsm.solve(rnproblem, pd,  dsddpsol);

    cout << "lb: " << dsddpsol.lb() << " ubm: " << dsddpsol.ubm()
           <<     " ubb: " << dsddpsol.ubb() << endl;
      for(unsigned int i=0; i<dsddpsol.firststage().size(); i++)
          cout << "x" << i << "=" << dsddpsol.firststage()[i] << endl;

/*
    sddpmethod<mpmcvarequivalent<psproblem>,gmdddistribution<double>> sm;


    sddpsolution<mpmcvarequivalent<psproblem>> sddpsol(mcvproblem);

    double sddpovalue;
    vector<const gmdddistribution<double>*> dv;
    dv.push_back(&d);

    sm.solve(mcvproblem,dv, vector<double>(0), sddpovalue, sddpsol);

    if(fabs(sddpovalue-cvarobj) > tol)
    {
        std::cerr << "cvartest: opt="
             << cvarobj << " expected, " << sddpovalue << " achieved." << std::endl;

        throw;
    }
*/
}

#endif // CVARTEST_H
