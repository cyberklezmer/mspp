#include "mspptest.h"
#include "mspp/de.h"

class tsproblem: public msproblem<realvar,linearfunction,
                fulllinearconstraint,fullhistory<omega>>
{
public:
    tsproblem() :
        msproblem<realvar,linearfunction,
                        fulllinearconstraint,fullhistory<omega>>
           (msproblemstructure({2,2}))
    {}
private:
    virtual void stageinfo_is(
            unsigned int k,
            const fullhistory<omega>& xi,
            vardefs<realvar>& xs,
            msconstraints<fulllinearconstraint>& g
            ) const
    {
        if(k==0)
        {
            xs[0].set(realvar::Rplus);
            xs[1].set(realvar::Rplus);
        }
        else if(k==1)
        {
            xs[0].set(realvar::Rplus);
            xs[1].set(realvar::Rplus);

            fulllinearconstraint& c=addg(g,k);
            c.set({xi[1].o1,1,1,0},constraint::geq, 7.0);

            fulllinearconstraint& d = addg(g,k);
            d.set({xi[1].o2,1,0,1},constraint::geq, 4.0);
        }

    }

    virtual linearfunction f_is(
                    unsigned int k,
                    const fullhistory<omega>& xi) const
    {
        return linearfunction({1.0,1.0});
    }
};


void twostagetest(unsigned int N, const lpsolver& cps)
{
    std::cout << "Twostagetest..." << std::endl;

    const unsigned int numleaves = N;
    std::vector<omega> items(numleaves*numleaves);

    for(unsigned int i=0; i<numleaves; i++)
    {
        double o1 = 1.0 + 2.0* (double)i / (numleaves-1);
        for(unsigned int j=0; j<numleaves; j++)
        {
           double o2 = 1.0/3.0 + 2.0/3.0* (double)j / (numleaves-1);
           items[3*i+j] = {o1, o2};
        }
    }

    guiddistribution<omega> g(items);

    processdistribution<guiddistribution<omega>> pd({0,0},g,1);

    using myscenariotree=distrscenariotree<guiddistribution<omega>>;
    myscenariotree sp(pd);

    tsproblem prp;

//    printvarnames(*prp);

    demethod<tsproblem,myscenariotree> b;

    stsolution<tsproblem,myscenariotree> sol(prp,sp);
    double ov;

    b.solve(prp,sp, cps, ov, sol);


    // brought from testproblems.xlsx
    // x0_0,x0_1,x1_0_0,x1_1_0,x1_0_1,x1_1_1,x1_0_2,x1_1_2,x1_0_3,x1_1_3,x1_0_4,x1_1_4,x1_0_5,x1_1_5,x1_0_6,x1_1_6,x1_0_7,x1_1_7,x1_0_8,x1_1_8,,,

    double tssol[]={2.25,2.5,2.25,0.75,2.25,0,2.25,0,3.33066907387547E-16,0.75,3.33066907387547E-16,0,0,0,0,0.75,0,0,0,0};
    double tsobj=5.75;

    const double tol = 1e-5;
    if(fabs(ov-tsobj) > tol)
    {
        std::cerr << "twostagetest: opt="
             << tsobj << " expected, " << ov << " achieved." << std::endl;
        throw;
    }

    unsigned int k = sizeof(tssol)/sizeof(tssol[0]);
    for(unsigned int i=0; i<k; i++)
    {
        if(fabs(sol.x(i)-tssol[i]) > tol)
        {
            std::cerr << "twostagetest: x[" << i << "]="
                 << tssol[i] << " expected, " << sol.x(i)
                 << " achieved." << std::endl;
            throw;
        }
    }
//    printstats(*s);
    std::cout << "Twostagetest passed." << std::endl;

/*    sddpmethod<tsproblem,myscenariotree> sm;

    sddpsolution ss;
    double sov;

    sm.solve(prp,sp, cps, sov, ss);*/

}



