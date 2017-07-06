#include <iostream>
#include <cmath>
#include "probability.h"
#include "scenarios.h"
#include "lpsolvers.h"
#include "linear.h"
#include "tests.h"


struct omega
{
    double o1;
    double o2;
};

class tsproblem: public linearproblem<omega>
{

public:
    tsproblem() : linearproblem<omega>({2,2})
    {}
    virtual void f(unsigned int stage,
                   const scenario<omega>& xih,
                   linearfunction_ptr& r) const
    {
       r.reset(new linearfunction({1.0,1.0}));
    }

    virtual void constraints(
                unsigned int stage,
                const scenario<omega>& xi,
                varinfo_list_ptr& vars,
                constraint_list_ptr<linearconstraint>& constraints
                ) const
    {
        if(stage==0)
        {
            vars.reset(new varinfo_list);
            vars->push_back(varinfo(varinfo::Rplus));
            vars->push_back(varinfo(varinfo::Rplus));
        }
        else if(stage==1)
        {
            vars.reset(new varinfo_list);
            vars->push_back(varinfo(varinfo::Rplus));
            vars->push_back(varinfo(varinfo::Rplus));

            constraints.reset(new linearconstraint_list);
            constraints->push_back(linearconstraint({xi[1].o1,1,1,0},
                          linearconstraint::geq, 7.0));
            constraints->push_back(linearconstraint({xi[1].o2,1,0,1},
                          linearconstraint::geq, 4.0));
        }
    }
};


void twostagetest(unsigned int N, const lpsolver& cps)
{
    std::cout << "Twostagetest..." << std::endl;

    const unsigned int numleaves = N;
    indexedtree_ptr tp(new nxtree({1,numleaves*numleaves}));
    uniformtreeprobability_ptr pp(new uniformtreeprobability(tp));

    generaltreemapping_ptr<omega> xp(new generaltreemapping<omega>(tp) );

    for(unsigned int i=0; i<numleaves; i++)
    {
        double o1 = 1.0 + 2.0* (double)i / (numleaves-1);
        for(unsigned int j=0; j<numleaves; j++)
        {
           double o2 = 1.0/3.0 + 2.0/3.0* (double)j / (numleaves-1);
           (*xp)({0,3*i+j}) = {o1, o2};
       }
    }

    scenariotree_ptr<omega> sp(new modularscenariotree<omega>(xp,pp));
    linearproblem_ptr<omega> prp(new tsproblem());

    biglpsolution<omega> b(prp,sp);
    treesolution_ptr s;
    double ov;
    b.solve(cps,s,ov);

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
        if(fabs(s->x(i)-tssol[i]) > tol)
        {
            std::cerr << "twostagetest: x[" << i << "]="
                 << tssol[i] << " expected, " << s->x(i)
                 << " achieved." << std::endl;
            throw;
        }
    }
    printstats(*s);
    std::cout << "Twostagetest passed." << std::endl;
}


class psproblem: public linearproblem<omega>
{

public:
    psproblem() : linearproblem<omega>({2,2})
    {}
    virtual void f(unsigned int stage,
                   const scenario<omega>& xih,
                   linearfunction_ptr& r) const
    {
       r.reset(new linearfunction(2));
       if(stage)
       {
           r->coefs[0] = -xih[1].o1;
           r->coefs[1] = -xih[1].o2;
       }
       else
       {
           r->coefs[0] = 0;
           r->coefs[1] = 0;
       }
    }

    virtual void constraints(
                unsigned int stage,
                const scenario<omega>& xi,
                varinfo_list_ptr& vars,
                constraint_list_ptr<linearconstraint>& constraints
                ) const
    {
        assert(stage<=1);

        vars.reset(new varinfo_list);
        vars->push_back(varinfo(varinfo::Rplus));
        vars->push_back(varinfo(varinfo::Rplus));

        constraints.reset(new linearconstraint_list);

        if(stage == 0)
            constraints->push_back(linearconstraint({1.0,1.0},
                          linearconstraint::eq, 1.0));
        else
        {
            constraints->push_back(linearconstraint({1.0,0,-1.0,0},
                          linearconstraint::eq, 0));
            constraints->push_back(linearconstraint({0,1.0,0,-1.0},
                          linearconstraint::eq, 0));
        }
    }
};


void cvartest(double alpha, double lambda, const lpsolver& cps)
{
    std::cout << "CVaRtext..." << std::endl;
    const int numleaves = 3;
    indexedtree_ptr tp(new nxtree({1,numleaves*numleaves}));
    uniformtreeprobability_ptr pp(new uniformtreeprobability(tp));

    generaltreemapping_ptr<omega> xp(new generaltreemapping<omega>(tp) );

    for(unsigned int i=0; i<numleaves; i++)
    {
        double o1 =  (double)i / numleaves;
        for(unsigned int j=0; j<numleaves; j++)
        {
           double o2 = 0.1 + (double) j / numleaves;
           (*xp)({0,numleaves*i+j}) = {o1, o2};
        }
    }

    scenariotree_ptr<omega> sp(new modularscenariotree<omega>(xp,pp));
    linearproblem_ptr<omega> prp(new psproblem());

    linearproblem_ptr<omega> cvp(new mmpcvarproblem<omega>(prp,alpha,lambda));

    biglpsolution<omega> b(cvp ,sp);

    treesolution_ptr s;
    double ov;
    b.solve(cps,s,ov);

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
        if(fabs(s->x(i)-cvarsol[i]) > tol)
        {
            std::cerr << "cvartest: x[" << i << "]="
                 << cvarsol[i] << " expected, " << s->x(i)
                 << " achieved." << std::endl;

            for(unsigned int j=0; j<s->nx(); j++)
                std::cerr << "x[" << j << "]=" << s->x(j) << std::endl;

            throw;
        }
    }

    unsigned int l = sizeof(cvarthetas)/sizeof(cvarthetas[0]);
    for(unsigned int i=0; i<l; i++)
    {
        if(fabs(s->x(k+2+3*i)-cvarthetas[i]) > tol)
        {
            std::cerr << "cvartest: x[" << i << "]="
                 << cvarthetas[i] << " expected, " << s->x(k+2+3*i)
                 << " achieved." << std::endl;

            for(unsigned int j=0; j<s->nx(); j++)
                std::cerr << "x[" << j << "]=" << s->x(j) << std::endl;

            throw;
        }
    }
    printstats(*s);
    std::cout << "CVaRtest passed." << std::endl;
}



class almproblem: public linearproblem<double>
{

public:
    almproblem() : linearproblem<double>({1,1,1,1})
    {}
    virtual void f(unsigned int stage,
                   const scenario<double>& xi,
                   linearfunction_ptr& r) const
    {
       r.reset(new linearfunction(1));
       double x=xi[0];
       for(unsigned int i=1; i<=stage; i++)
          x *= xi[i];
       r->coefs[0] = x;
    }

    virtual void constraints(
                unsigned int stage,
                const scenario<double>& xi,
                varinfo_list_ptr& vars,
                constraint_list_ptr<linearconstraint>& constraints) const
    {
        vars.reset(new varinfo_list);
        if(stage==0)
        {
           vars->push_back(varinfo(0,1,"x"));
        }
        else if(stage==1)
        {
           vars->push_back(varinfo(varinfo::Rplus,"x"));
           constraints.reset(new linearconstraint_list);
           constraints->push_back(linearconstraint({1.0,1.0},
                                    linearconstraint::geq, 0.0));
           constraints->push_back(linearconstraint({1.0,1.0},
                                    linearconstraint::leq, 1.0));
        }
        else if(stage==2)
        {
           vars->push_back(varinfo(varinfo::Rplus,"x"));
           constraints.reset(new linearconstraint_list);
           constraints->push_back(linearconstraint({1.0,1.0,1.0},
                                    linearconstraint::geq, 0.0));
           constraints->push_back(linearconstraint({1.0,1.0,1.0},
                                    linearconstraint::leq, 1.0));
        }
        else
        {
            vars->push_back(varinfo(varinfo::Rplus,"x"));
            constraints.reset(new linearconstraint_list);
            constraints->push_back(linearconstraint({1.0,1.0,1.0,1.0},
                                     linearconstraint::eq, 1.0));
        }
    }
};


void almtest(double alpha, double lambda, const lpsolver& cps)
{
    std::cout << "ALMtest..." << std::endl;
    const int numleaves = 2;
    homogeneoustree_ptr tp(new homogeneoustree(3,numleaves));
    iidtreeprobability_ptr pp(new iidtreeprobability(tp,{0.5,0.5}));
    iidtreemapping_ptr<double> xp(new iidtreemapping<double>(tp,{0.8,1.1}));

    shiftedscenariotree_ptr<double> ss(new shiftedscenariotree<double>(xp,pp,1.0));

    linearproblem_ptr<double> prp(new almproblem());

    linearproblem_ptr<double> cvp(new mmpcvarproblem<double>(prp,alpha,lambda));

    biglpsolution<double> b(cvp,ss);

    treesolution_ptr s;
    double ov;
    b.solve(cps,s,ov);

    const double tol = 1e-5;


    double almsol[]={0,0,0,0,0.64,1,0,0,0,0,0,0,0.727273,0,0.264,0.272727,-0.036,0.272727,0,0,
                     0,0.88,1,0,0,0,0,0,0,0.727273,0,0.363,0.272727,-0.0495,0.272727,0};
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
        if(fabs(s->x(i)-almsol[i]) > tol)
        {
            std::cerr << "almtest: x[" << i << "]="
                 << almsol[i] << " expected, " << s->x(i)
                 << " achieved." << std::endl;
            std::cerr << "x=";

            for(unsigned int j=0; j<s->nx(); j++)
            {
                std::cerr << s->x(j) << ",";
                if((j+1)%20 == 0)
                     std::cerr << std::endl;
            }
            throw;
        }
    }
    printstats(*s);
    std::cout << "ALMtest passed." << std::endl;
}


int main(int argc, char *argv[])
{
  sys::set_log("mspp.log");

  cplexlpsolver cps;
  csvlpsolver csvs;

  twostagetest(3,cps);
  cvartest(0.05,0.5,cps);
  almtest(0.05,0.5,cps);
}
