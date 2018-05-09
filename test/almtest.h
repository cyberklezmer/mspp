#ifndef ALMTEST_H
#define ALMTEST_H

#include "mspptest.h"
#include "mspp/de.h"

class almproblem: public msproblem<mpmcvar, linearfunction,
        linearmsconstraint,realvar,double>
{

public:
    almproblem(double lambda, double alpha) : msproblem<mpmcvar, linearfunction,
                   linearmsconstraint,realvar,double>
        ({1,1,1,1},mpmcvar(lambda,alpha))
    {}



    virtual linearfunction f_is(unsigned int k,
            const subvectors<double>& barxi) const
    {
       double c=barxi[0][0];
//std::cout << "c="        << c << std::endl;
       for(unsigned int i=1; i<=k; i++)
          c *= barxi[i][0];
       linearfunction r = newf(k);
       r.setc(k,c);
       return r;
    }

    virtual void x_is(
            unsigned int k,
            const subvectors<double>& xi,
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

template <typename O>
void almtest(double alpha, double lambda)
{
    const int numleaves = 2;

    gdscenariotree<double> st({1},gddistribution({0.8,1.1}),3);

    almproblem prp(lambda,lambda);

    std::cout <<"ALMtest..." << std::endl;

//    printvarnames(*cvp);
    demethod<almproblem,gdscenariotree<double>, O> b;
    stsolution<almproblem,gdscenariotree<double>> s(prp,st);

    double ov;
    b.solve(prp,st,ov,s);

    vector<double> x = s.x();
//    printstats(*s);

    const double tol = 1e-5;

//    double almsol[]={0,0,0,0,0.64,1,0,0,0,0,0,0,0.727273,0,0.264,0.272727,-0.036,0.272727,0,0,
//                     0,0.88,1,0,0,0,0,0,0,0.727273,0,0.363,0.272727,-0.0495,0.272727,0};
    double almsol[]={0,0,1,0,0,0.727273,0.272727,0.272727,0,1,0,0,0.727273,0.272727,0.272727};
    double almobj = 0.906063;

    if(fabs(ov-almobj) > tol)
    {
        std::cerr << "almtest: opt="
             << almobj << " expected, " << ov << " achieved." << std::endl;
        throw;
    }

    unsigned int k = sizeof(almsol)/sizeof(almsol[0]);
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
}

#endif // ALMTEST_H
