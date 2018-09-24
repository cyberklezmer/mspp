#ifndef ALMTEST_H
#define ALMTEST_H

#include "mspptest.h"
#include "mspp/sddp.h"
#include "mspp/de.h"


class almproblem: public msproblem<mpmcvar, linearfunction,
        linearmsconstraint,realvar,pair<unsigned int,rvector<double>>,lastx,lastmdxi<rvector<double>>>
{

public:
    almproblem(double lambda, double alpha) :
        msproblem<mpmcvar, linearfunction,
        linearmsconstraint,realvar,pair<unsigned int,rvector<double>>,
        lastx,lastmdxi<rvector<double>>>
        ({1,2,2,2},mpmcvar(lambda,alpha))
    {}


    virtual void f_is(unsigned int k,
                      const rvector<double>& zeta,
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
            const rvector<double>& zeta,
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
        if(k==1)
            g.add(linearmsconstraint({1.0,1.0,-1.0},constraint::eq, 0.0));
        else
            g.add(linearmsconstraint({0.0,1.0,1.0,-1.0},constraint::eq, 0.0));
        if(k==this->T())
            xs[1].setlimits(1,1);
    }

};

class almdistribution: public
        mddistribution<double, unsigned int>, public
        ddistribution<rvector<double>, unsigned int>
{
public:
    using I_t = rvector<double>;
    using C_t = unsigned int ;
    almdistribution(unsigned int n, const vector<double>& sc) :
        mddistribution<double, unsigned int>(1), fn(n), fsc(sc) {}
private:
    virtual atom<rvector<double>> atom_is(unsigned int i, const unsigned int& c) const
    {
        assert(c < fsc.size());
        return { {fsc[i] + (double) i / (double) fn}, 1.0/ (double) fn };
    }

    virtual unsigned int natoms_is(const unsigned int& c) const
    { return fn; }

    unsigned int fn;
    vector<double> fsc;
};

template <typename O>
void almtest()
{
    const double tol = 1e-5;
    double sol[]={0,0,1,0,0,0.727273,0.272727,0.272727,0,1,0,0,0.727273,0.272727,0.272727};
    double obj = 0.906063;
    unsigned int k = sizeof(sol)/sizeof(sol[0]);

    unsigned int nl = 1;

    // tbd soution nemusí vědět o scénáři?

/*    for(unsigned int i=0; i<k; i++)
    {
        if(fabs(x.x(i)-sol[i]) > tol)
        {
            std::cerr << "x[" << i << "]="
                 << sol[i] << " expected, " << x.x(i)
                 << " achieved." << std::endl;
            std::cerr << endl;
            throw;
        }
    }
    std::cout <<  "Passed." << std::endl;
*/

    using mcd = diterativedistribution<mmarkovchain,almdistribution>;
    vector<mcd> dm;
    dm.push_back(mcd(mmarkovchain({{0.5,0.5}}),
                     almdistribution(nl,{ 0.8, 1.1}))
                 );

    dm.push_back(mcd(mmarkovchain({{0.5,0.5,0},{0,0.5,0.5}}),
                     almdistribution(nl,{ 0.8*0.8, 1.1*0.8, 1.1*1.1 }))
                     );
    dm.push_back(mcd(mmarkovchain({{0.5,0.5,0,0},{0,0.5,0.5,0}
                                     ,{0,0,0.5,0.5}}),
                     almdistribution(nl,{ 0.8*0.8*0.8, 0.8*0.8*1.1,
                                          0.8*1.1*1.1, 1.1*1.1*1.1 } ))
                 );

    almproblem mp(0.5,0.05);

// tbd dávat de jen distribution
    using pdt = processdistribution<mcd,
                   laststate<typename mcd::I_t>>;

    pdt pd({0,{0}},dm);
    using mctree = distrscenariotree<pdt>;

    mctree s(pd);

    std::cout <<"MALMtest DE..." << std::endl;

    stsolution<almproblem,mctree> x(mp,s);

    demethod::solve<almproblem,mctree, O>(mp,s,x);

    if(fabs(x.obj()-obj) > tol)
    {
        std::cerr << "almtest: opt="
             << obj << " expected, " << x.obj() << " achieved." << std::endl;
        throw;
    }

    std::cout <<  "passed." << std::endl;

    return;

    /*std::cout <<"MALMtest DE..." << std::endl;

    msddpprocessdist<mmarkovchain, almdistribution> pd({0},dm);
    sddpsolution<malmproblem> ss(mp);

//sddpsolution<malmproblem> ss2 = ss;

    msddpmethod<malmproblem,msddpprocessdist<mmarkovchain, almdistribution>> ssm;

    ssm.solve(mp,pd,ss);
    if(fabs(ss.lb()-almobj) > 0.1)
    {
        std::cerr << "malmtest: opt="
             << almobj << " expected, " <<
             ss.lb() << "<" << ss.ubm() <<
             " achieved." << std::endl;
        throw;
    }
*/

}
#endif // ALMTEST_H
