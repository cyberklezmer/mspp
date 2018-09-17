#ifndef ALMTEST_H
#define ALMTEST_H

#include "mspptest.h"
#include "mspp/de.h"

class almproblem: public msproblem<mpmcvar, linearfunction,
        linearmsconstraint,realvar>
{

public:
    almproblem(double lambda, double alpha) : msproblem<mpmcvar, linearfunction,
                   linearmsconstraint,realvar>
        ({1,1,1,1},mpmcvar(lambda,alpha))
    {}


    virtual void f_is(unsigned int k,
            const vector<double>& zeta,
                      linearfunction f) const
    {
       double c=1;
       for(unsigned int i=0; i <= k; i++)
          c *= zeta[i];
       f.setc(0,c);
    }

    virtual void x_is(
            unsigned int k,
            const vector<double>& ,
            ranges<realvar>& xs,
            msconstraints<linearmsconstraint>& g) const
    {
        if(k==0)
           xs[0].setlimits(0,1);
        else if(k==1)
        {
           xs[0].setpositive();
           g.add(linearmsconstraint({1.0,1.0},constraint::geq, 0.0));
           g.add(linearmsconstraint({1.0,1.0},constraint::leq, 1.0));
        }
        else if(k==2)
        {
           xs[0].setpositive();
           g.add(linearmsconstraint({1.0,1.0,1.0},
                                    constraint::geq, 0.0));
           g.add(linearmsconstraint({1.0,1.0,1.0},
                                    constraint::leq, 1.0));
        }
        else
        {
            xs[0].setpositive();
            g.add(linearmsconstraint({1.0,1.0,1.0,1.0},
                                     constraint::eq, 1.0));
        }
    }
};

class malmproblem: public msproblem<mpmcvar, linearfunction,
        linearmsconstraint,realvar,rvector<double>,lastx,lastxi<rvector<double>>>
{

public:
    malmproblem(double lambda, double alpha) :
        msproblem<mpmcvar, linearfunction,
        linearmsconstraint,realvar,rvector<double>,
        lastx,lastxi<rvector<double>>>
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
    double almsol[]={0,0,1,0,0,0.727273,0.272727,0.272727,0,1,0,0,0.727273,0.272727,0.272727};
    double almobj = 0.906063;
    unsigned int k = sizeof(almsol)/sizeof(almsol[0]);


/*    gdscenariotree<double> st({1},gddistribution({0.8,1.1}),3);

    almproblem prp(0.5,0.02);

    std::cout <<"ALMtest..." << std::endl;

    demethod<almproblem,gdscenariotree<double>, O> b;
    stsolution<almproblem,gdscenariotree<double>> s(prp,st);

    double ov;
    b.solve(prp,st,ov,s);



    if(fabs(ov-almobj) > tol)
    {
        std::cerr << "almtest: opt="
             << almobj << " expected, " << ov << " achieved." << std::endl;
        throw;
    }

    for(unsigned int i=0; i<k; i++)
    {
        if(fabs(s.x(i)-almsol[i]) > tol)
        {
            std::cerr << "almtest: x[" << i << "]="
                 << almsol[i] << " expected, " << s.x(i)
                 << " achieved." << std::endl;
            std::cerr << "x=";

            throw;
        }
    }
    std::cout <<  "Passed." << std::endl;
*/
    unsigned int nl = 1;
// tbd udělat pro něj dobrej draw
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
    malmproblem mp(0.5,0.05);

// tbd dávat de jen distribution
    using mctree = distrscenariotree<diterativedistribution<mmarkovchain,almdistribution>>;

    mctree mt(0,dm);


    std::cout <<"MALMtest DE..." << std::endl;

    demethod<malmproblem,mctree, O> mb;
    stsolution<malmproblem,mctree> ms(mp,mt);

    double mov;
    mb.solve(mp,mt,mov,ms);

    if(fabs(mov-almobj) > tol)
    {
        std::cerr << "almtest: opt="
             << almobj << " expected, " << mov << " achieved." << std::endl;
        throw;
    }


    std::cout <<  "Passed." << std::endl;

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
