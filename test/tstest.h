#ifndef TSTEST_H
#define TSTEST_H

#include "mspptest.h"
#include "mspp/de.h"

struct omega
{
    double o1;
    double o2;
};

template <typename Z>
class tsproblem: public msproblem<expectation,linearfunction,
                linearmsconstraint,realvar,omega,allx,Z>
{
public:
    tsproblem() :
        msproblem<expectation,linearfunction,
        linearmsconstraint,realvar,omega,allx,Z>
           (msproblemstructure({2,2}))
    {}
private:

    virtual void x_is(
            unsigned int k,
            const omega& o,
            ranges<realvar>& xs,
            msconstraints<linearmsconstraint>& g
            ) const
    {
        if(k==0)
        {
            xs[0].setpositive();
            xs[1].setpositive();
        }
        else if(k==1)
        {
            xs[0].setpositive();
            xs[1].setpositive();

            linearmsconstraint& c=g.add();

            c.set({o.o1,1,1,0},constraint::geq, 7.0);

            linearmsconstraint& d = g.add();
            d.set({o.o2,1,0,1},constraint::geq, 4.0);
        }

    }

    virtual void f_is(unsigned int k,
                                const omega& zeta,
                                linearfunction& f
                                ) const
    {
        f = linearfunction({1.0,1.0});
    }
};

template <typename R, typename O>
void tst()
{
    unsigned int N=3;

    const unsigned int numleaves = N;
    vector<omega> items(numleaves*numleaves);

    for(unsigned int i=0; i<numleaves; i++)
    {
        double o1 = 1.0 + 2.0* (double)i / (numleaves-1);
        for(unsigned int j=0; j<numleaves; j++)
        {
           double o2 = 1.0/3.0 + 2.0/3.0* (double)j / (numleaves-1);
           items[3*i+j] = {o1, o2};
        }
    }

    gddistribution g(items);

    processdistribution<gddistribution<omega>> pd({0,0},g,1);

//    using myscenariotree=distrscenariotree<guiddistribution<double>>;
    gdscenariotree<omega> sp(pd);

    tsproblem<R> prp;

    demethod<tsproblem<R>,gdscenariotree<omega>,O> b;

    stsolution<tsproblem<R>,gdscenariotree<omega>> sol(prp,sp);

    double ov;

    b.solve(prp,sp, ov, sol);


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

    std::cout << "Passed." << std::endl;


/*    sddpmethod<tsproblem,myscenariotree> sm;

    sddpsolution ss;
    double sov;

    sm.solve(prp,sp, cps, sov, ss);
*/
}

template <typename O>
void twostagetest()
{
    std::cout << "Twostagetest with Rxi=lastxi" << std::endl;
    tst<lastxi<omega>,O>();
}


#endif // TSTEST_H
