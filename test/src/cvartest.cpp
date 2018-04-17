#include "mspptest.h"
#include "mspp/mpcproblem.h"
#include "mspp/de.h"

class psproblem: public msproblem<realvar,linearfunction,
        fulllinearconstraint,fullhistory<omega>>
{

public:template <typename X>
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

    psproblem() : msproblem<realvar,linearfunction,
                  fulllinearconstraint,fullhistory<omega>>
          ({2,2})
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

    virtual void stageinfo_is(
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


void cvartest(double alpha, double lambda, const lpsolver& cps)
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


    guiddistribution<omega> g(items);

    processdistribution<guiddistribution<omega>> pd({0,0},g,1);

    using myscenariotree=distrscenariotree<guiddistribution<omega>>;
    myscenariotree sp(pd);


    psproblem prp;

    std::cout << "CVaRtest "  << std::endl;

    mmpcvarproblem<psproblem> cvp(prp,alpha,lambda);

    demethod<mmpcvarproblem<psproblem>,myscenariotree> b;

    stsolution<mmpcvarproblem<psproblem>,myscenariotree> sol(cvp,sp);
    double ov;

    b.solve(cvp,sp, cps, ov, sol);

    const double tol = 1e-5;

    // brought from testproblems.xlsx
    // x0_1,x_02,u0,theta_1,theta_2,theta_3,theta_4,theta_5,theta_6,theta_7,theta_8,theta_9,,,

    double cvarsol[]={0,1,-0.1}; // x0_1,x_02,u0
    double cvarthetas[]={0,-0.166666666666666,-0.333333333333335,
                      0,-0.166666666666669,-0.333333333333333,
                      0,-0.166666666666667,-0.333333333333332};

    double cvarobj = -0.266666666666667;

    if(fabs(ov-cvarobj) > tol)
    {
        std::cerr << "cvartest: opt="
             << cvarobj << " expected, " << ov << " achieved." << std::endl;

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
    std::cout <<  "MPCVaRtest passed." << std::endl;
}

