#include <cmath>
#include "mspp/cplex.h"
#include "mspp/de.h"
#include "tests.h"
#include "mspp/sddp.h"
#include "mspp/mpcproblem.h"

// milp do templatu

using namespace mspp;

template <typename X>
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



struct omega
{
    double o1;
    double o2;
};

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


class psproblem: public msproblem<realvar,linearfunction,
        fulllinearconstraint,fullhistory<omega>>
{

public:
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



int main(int argc, char *argv[])
{
//    csvlpsolver csvs;
    cplexlpsolver csvs;
//    twostagetest(3,csvs);
    cvartest(0.05,0.5,csvs);
}


/*
#include <iostream>
#include "mspp/biglp.h"
#include "mspp/mmpcvar.h"
#include "mspp/cplex.h"
#include "mspp/shifted.h"
#include "tests.h"





class almproblem: public linearproblem<double>
{

public:
    almproblem() : linearproblem<double>({1,1,1,1})
    {}
    virtual void f(unsigned int stage,
                   const scenario<double>& xi,
                   linearobjective& r) const
    {
       double c=xi[0];
       for(unsigned int i=1; i<=stage; i++)
          c *= xi[i];
       r.coefs[0] = c;
    }

    virtual void msconstraints(
                unsigned int stage,
                const scenario<double>& xi,
                varranges& vars,
                constraint_list<linearconstraint>& csts) const
    {
        if(stage==0)
           vars[0]=varrange(0,1);
        else if(stage==1)
        {
           vars[0]=varrange(varrange::Rplus);
           csts.add(linearconstraint({1.0,1.0},linearconstraint::geq, 0.0));
           csts.add(linearconstraint({1.0,1.0},linearconstraint::leq, 1.0));
        }
        else if(stage==2)
        {
           vars[0]=varrange(varrange::Rplus);
           csts.add(linearconstraint({1.0,1.0,1.0},
                                    linearconstraint::geq, 0.0));
           csts.add(linearconstraint({1.0,1.0,1.0},
                                    linearconstraint::leq, 1.0));
        }
        else
        {
            vars[0]=varrange(varrange::Rplus);
            csts.add(linearconstraint({1.0,1.0,1.0,1.0},
                                     linearconstraint::eq, 1.0));
        }
    }
};

void almtest(double alpha, double lambda, const lpsolver_ptr& cps)
{
    const int numleaves = 2;
    ntree_ptr tp(new ntree(3,numleaves));
    iidprobability_xxx_ptr pp(new iidprobability_xxx(tp,{0.5,0.5}));
    iidprocess_ptr<double> xp(new iidprocess<double>(tp,{0.8,1.1}));

    shiftedscenariotree_ptr<double> ss(new shiftedscenariotree<double>(xp,pp,1.0));

    linearproblem_ptr<double> prp(new almproblem());

    std::cout <<"ALMtest..." << std::endl;

    linearproblem_ptr<double> cvp;
    cvp.reset(new mmpcvarproblem<double>(prp,alpha,lambda));

    printvarnames(*cvp);
    biglpmethod<double> b(cvp,ss,cps);

    treesolution_ptr s;
    double ov;
    b.solve(s,ov);

    printstats(*s);

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
                 <<https://slovnik.seznam.cz/en/?q=z%C3%A1stupce almsol[i] << " expected, " << s->x(i)
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
    std::cout <<  " ALMtest passed." << std::endl;
}


int main(int argc, char *argv[])
{
//  sys::set_log("mspp.log");

  lpsolver_ptr cps(new cplexlpsolver);

  csvlpsolver csvs;

  almtest(0.05,0.5,cps);
}

*/

