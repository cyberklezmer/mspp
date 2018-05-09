#ifndef TSTEST_H
#define TSTEST_H

#include "mspptest.h"
#include "mspp/de.h"


template <typename Rxi>
class tsproblem: public msproblem<expectation,linearfunction,
                linearmsconstraint,realvar,double,everything,Rxi>
{
public:
    tsproblem() :
        msproblem<expectation,linearfunction,
        linearmsconstraint,realvar,double,everything,Rxi>
           (msproblemstructure({2,2}))
    {}
private:

    virtual void x_is(
            unsigned int k,
            const subvectors<double>& xi,
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

            linearmsconstraint& c=this->addg(g,k);
            c.set({xi[1][0],1,1,0},constraint::geq, 7.0);

            linearmsconstraint& d = this->addg(g,k);
            d.set({xi[1][1],1,0,1},constraint::geq, 4.0);
        }

    }

    virtual linearfunction f_is(unsigned int k,
            const subvectors<double>& barxi) const
    {
        return k ? linearfunction({0,0,1.0,1.0})
                 : linearfunction({1,1});
    }
};

template <typename R, typename O>
void tst()
{
    unsigned int N=3;

    const unsigned int numleaves = N;
    vector<rvector<double>> items(numleaves*numleaves);

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

    processdistribution<gddistribution<double>> pd({0,0},g,1);

//    using myscenariotree=distrscenariotree<guiddistribution<double>>;
    gdscenariotree<double> sp(pd);

    tsproblem<R> prp;

    demethod<tsproblem<R>,gdscenariotree<double>,O> b;

    stsolution<tsproblem<R>,gdscenariotree<double>> sol(prp,sp);

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
    std::cout << "Twostagetest with Rxi=everything" << std::endl;
    tst<everything,O>();
    std::cout << "Twostagetest with Rxi=nothing" << std::endl;
    tst<nothing,O>();
}


#endif // TSTEST_H
