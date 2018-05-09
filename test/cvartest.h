#ifndef CVARTEST_H
#define CVARTEST_H

#include "mspptest.h"
#include "mspp/mpcproblem.h"
#include "mspp/de.h"

class psproblem: public msproblem< mpmcvar,
        linearfunction,
        linearmsconstraint,realvar,double>
{
    msproblemstructure initps(bool ssv)
    {
        return ssv ? msproblemstructure({2,2}) : msproblemstructure({2,0});
    }

public:
    psproblem(double lambda, double alpha, bool secstagevars=true) :
           msproblem<mpmcvar,
           linearfunction,
           linearmsconstraint,realvar,double>
          (initps(secstagevars),mpmcvar(lambda,alpha)), fssv(secstagevars)
    {}

    virtual linearfunction f_is(
                    unsigned int k,
                    const subvectors<double>& bx) const
    {
       if(fssv)
       {
          if(k)
             return linearfunction({0,0,-bx[1][0],-bx[1][1]});
          else
             return linearfunction({0,0});
       }
       else
       {
           if(k)
               return linearfunction({-bx[1][0],-bx[1][1]});
           else
               return linearfunction({0,0});
       }
    }



    virtual void x_is(
            unsigned int k,
            const subvectors<double>& bx,
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
        else if(fssv)
        {
            r[0].setpositive();
            r[1].setpositive();
            g.add(linearmsconstraint({1.0,0,-1.0,0},
                          constraint::eq, 0));
            g.add(linearmsconstraint({0,1.0,0,-1.0},
                          constraint::eq, 0));
        }
    }
private:
    bool fssv;
};

template <typename O>
void cvartest(double alpha, double lambda, bool ssv)
{

    gddistribution d = gddistribution({0, 1.0/3.0, 2.0 / 3.0})
                        * gddistribution({0.1+0, 0.1+ 1.0/3.0, 0.1+2.0 / 3.0});

    gdscenariotree<double> stree(d,1);

    psproblem rnproblem(lambda,alpha,ssv);

    const double tol = 1e-5;

    double cvarobj = -0.266666666666667;

    std::cout << "CVaRtest direct ";
    std::cout << (ssv ? "with" : "without") << " second stage vars." << std::endl;


    demethod<psproblem,gdscenariotree<double>,O> m;
    stsolution<psproblem,gdscenariotree<double>> s(rnproblem,stree);

    double o;

    m.solve(rnproblem, stree, o, s);

    if(fabs(o-cvarobj) > tol)
    {
        std::cerr << "cvartest direct: opt="
             << cvarobj << " expected, " << o << " achieved." << std::endl;

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

    std::cout << "CVaRtest indirect ";
    std::cout << (ssv ? "with" : "without") << " second stage vars." << std::endl;

    mpmcvarequivalent<psproblem> mcvproblem(rnproblem);

    demethod<mpmcvarequivalent<psproblem>,gdscenariotree<double>,O> b;

    stsolution<mpmcvarequivalent<psproblem>,gdscenariotree<double>>
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

    unsigned int thetafactor = ssv ? 3 : 1;
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
}

#endif // CVARTEST_H
