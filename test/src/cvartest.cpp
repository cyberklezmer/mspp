#include "mspptest.h"
#include "mspp/mpcproblem.h"
#include "mspp/de.h"

class psproblem: public msproblem<realvar,linearfunction,
        fulllinearconstraint,mmpcvar,fullhistory<omega>>
{

public:
    psproblem(double lambda, double alpha) : msproblem<realvar,linearfunction,
                  fulllinearconstraint,mmpcvar, fullhistory<omega>>
          ({2,2},mmpcvar(lambda,alpha))
    {}

    virtual linearfunction f_is(
                    unsigned int k,
                    const fullhistory<omega>& xi) const

    {
       if(k)
           return linearfunction({-xi[1].o1,-xi[1].o2});
       else
           return linearfunction({0,0});
    }

    virtual void xset_is(
            unsigned int k,
            const fullhistory<omega>& xi,
            vardefs<realvar>& r,
            msconstraints<fulllinearconstraint>& g
            ) const
    {
        assert(k<=1);

        r[0].set(realvar::Rplus);
        r[1].set(realvar::Rplus);

        if(k == 0)
            g.add(fulllinearconstraint({1.0,1.0},
                          constraint::eq, 1.0));
        else
        {
            g.add(fulllinearconstraint({1.0,0,-1.0,0},
                          constraint::eq, 0));
            g.add(fulllinearconstraint({0,1.0,0,-1.0},
                          constraint::eq, 0));
        }
    }
};


void cvartest(double alpha, double lambda, const lpsolver& solver)
{
    const int numleaves = 3;

    std::vector<omega> items(numleaves*numleaves);

    for(unsigned int i=0; i<numleaves; i++)
    {
        double o1 =  (double)i / numleaves;
        for(unsigned int j=0; j<numleaves; j++)
        {
           double o2 = 0.1 + (double) j / numleaves;
           items[3*i+j] = {o1, o2};
        }
    }

    guiddistribution<omega> distr(items);

    processdistribution<guiddistribution<omega>> pdistr({0,0},distr,1);

    using myscenariotree=distrscenariotree<guiddistribution<omega>>;

    myscenariotree stree(pdistr);


    psproblem raproblem(lambda,alpha);

    std::cout << "CVaRtest direct "  << std::endl;

    demethod<psproblem,myscenariotree> m;
    stsolution<psproblem,myscenariotree> s(raproblem,stree);

    double o;

    m.solve(raproblem,stree, solver, o, s);

    const double tol = 1e-5;

    double cvarobj = -0.266666666666667;

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

    std::cout << "CVaRtest direct passed."  << std::endl;

    std::cout << "CVaRtest indirect "  << std::endl;

    mmpcvarequivalent<psproblem> mcvproblem(raproblem);

    demethod<mmpcvarequivalent<psproblem>,myscenariotree> b;

    stsolution<mmpcvarequivalent<psproblem>,myscenariotree>
            sol(mcvproblem,stree);

    double ovalue;

    b.solve(mcvproblem,stree, solver, ovalue, sol);


    // brought from testproblems.xlsx
    // x0_1,x_02,u0,theta_1,theta_2,theta_3,theta_4,theta_5,theta_6,theta_7,theta_8,theta_9,,,

    double cvarsol[]={0,1,-0.1}; // x0_1,x_02,u0
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
    for(unsigned int i=0; i<l; i++)
    {
        if(fabs(sol.x(k+2+3*i)-cvarthetas[i]) > tol)
        {
            std::cerr << "cvartest: x[" << i << "]="
                 << cvarthetas[i] << " expected, " << sol.x(k+2+3*i)
                 << " achieved." << std::endl;

            throw;
        }
    }
    std::cout <<  "MPCVaRtest indirect passed." << std::endl;


}

