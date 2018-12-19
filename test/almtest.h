#ifndef ALMTEST_H
#define ALMTEST_H

#include "mspptest.h"
#include "mspp/sddp.h"
#include "mspp/de.h"


class almproblem: public msproblem<mpmcvar, linearfunction,
        linearmsconstraint,vector<double>,realvar,lastx>
{
    static std::vector<unsigned int> makeps(unsigned int t)
    {
        std::vector<unsigned int> res = {1};
        for(unsigned int i=1; i<=t; i++)
            res.push_back(2);
        return res;
    }
public:
    almproblem(double lambda, double alpha, unsigned int T) :
        msproblem<mpmcvar, linearfunction,
                linearmsconstraint,vector<double>,realvar,lastx>
        (makeps(T),mpmcvar(lambda,alpha))
    {}


    virtual void f_is(unsigned int k,
                      const vector<double>& zeta,
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
            const vector<double>& zeta,
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
    double minf_is(unsigned int) const
    {
        return -1000;
    }
    double maxf_is(unsigned int) const
    {
        return 1000;
    }
};

class almdistribution:
        public vdistribution<double, unsigned int>,
        public fdistribution<vector<double>, unsigned int>
{
public:
    using I_t = vector<double>;
    using C_t = unsigned int ;
    almdistribution(unsigned int n, const vector<double>& sc) :
        vdistribution<double, unsigned int>(1), fn(n), fsc(sc) {}
private:
    virtual atom<vector<double>> atom_is(unsigned int i, const unsigned int& c) const
    {
        assert(c < fsc.size());
        return { {fsc[c] - 0.5 +  (double) (2*i+1) / (double) (2*fn)}, 1.0/ (double) fn };
    }

    virtual unsigned int natoms_is(const unsigned int& c) const
    { return fn; }

    unsigned int fn;
    vector<double> fsc;
};

template <typename O>
void alm1test(bool equivalent)
{
    almproblem mp(0.5,0.05,1);
    lvdistribution<double> d({{0.5500},{0.8500}, {1.0500},{1.3500}});
    using dist = fdprocessdistribution<lvdistribution<double>,noxi<vector<double>>>;
    dist pd({1},d,1);

    if(!equivalent)
    {
      std::cout <<"MALMtest DE direct..." << std::endl;

      desolution<almproblem,dist,lastxi<vector<double>>,O> x(mp,pd);
      if(fabs(x.obj()-1) > 0.001)
      {
          std::cerr << "alm1test: opt="
               << 1 << " expected, " << x.obj() << " achieved." << std::endl;
          throw;
      }

      sddpsolution<almproblem,dist,lastxi<vector<double>>,O> sx(mp,pd);
      if(fabs(sx.obj().lb()-x.obj()) > 0.05)
      {
          std::cerr << "alm1test: opt="
               << x.obj() << " expected, " <<
               sx.obj().lb() << "<" << sx.obj().ubm() <<
               " achieved." << std::endl;
          throw;
      }
    }
    else
    {
          std::cout <<"MALMtest DE equivalent..." << std::endl;

          mpmcvarequivalent<almproblem> ep(mp);

          desolution<mpmcvarequivalent<almproblem>,dist,
              lastxi<vector<double>>,O> x(ep,pd);
          if(fabs(x.obj()-1) > 0.001)
          {
              std::cerr << "alm1test: opt="
                   << 1 << " expected, " << x.obj() << " achieved." << std::endl;
              throw;
          }

          sddpsolution<mpmcvarequivalent<almproblem>,dist,lastxi<vector<double>>,O> sx(ep,pd);
          if(fabs(sx.obj().lb()-x.obj()) > 0.05)
          {
              std::cerr << "alm1test: opt="
                   << x.obj() << " expected, " <<
                   sx.obj().lb() << "<" << sx.obj().ubm() <<
                   " achieved." << std::endl;
              throw;
          }
    }
    std::cout << "passed" << std::endl;
}

template <typename O>
void almtest(unsigned int T=3, unsigned int nl=1)
{
    assert(T>=1);
    const double tol = 1e-5;
    double sol[]={
            0,0,0,1,1,0,1,0,1,0,1,0,1,
            0.727273,0.727273,0.272727,1,0.272727,
            1,0.272727,1,0.272727,1,0,0,1,1,1,1,1,
            1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,1,0,1,0,
            1,0,1,0.727273,0.727273,0.272727,1,0.272727,
            1,0.272727,1,0.272727,1};
    double obj = 0.906063;
    unsigned int k = sizeof(sol)/sizeof(sol[0]);


    using mydist = fhmcdistribution<mmcdistribution,almdistribution>;
    vector<mydist> d;

    mydist d1(mmcdistribution({{0.5,0.5}}),almdistribution(nl,{ 0.8, 1.1}));

    d.push_back(d1);

    if(T>1)
        d.push_back(mydist(mmcdistribution({{0.5,0.5,0},{0,0.5,0.5}}),
            almdistribution(nl,{ 0.8*0.8, 1.1*0.8, 1.1*1.1 })));

    if(T>2)
        d.push_back(mydist(mmcdistribution({{0.5,0.5,0,0},{0,0.5,0.5,0}
                                     ,{0,0,0.5,0.5}}),
          almdistribution(nl,{ 0.8*0.8*0.8, 0.8*0.8*1.1,
                                          0.8*1.1*1.1, 1.1*1.1*1.1 } )));

    using pdt = fdprocessdistribution<mydist,laststate<typename almdistribution::I_t>>;
    pdt pd({0,{1}},d);

    almproblem mp(0.5,0.05,T);

    std::cout <<"MALMtest DE..." << std::endl;

    desolution<almproblem,pdt,lastmdxi<vector<double>>,O> x(mp,pd);

    if(T==3 && nl==1 && fabs(x.obj()-obj) > tol)
    {
        std::cerr << "almtest: opt="
             << obj << " expected, " << x.obj() << " achieved." << std::endl;
        throw;
    }

    vector<double> xs;
    x.x()->exportlinear(xs);


    if(T==3 && nl==1)
        for(unsigned int i=0; i<k; i++)
        {
            if(fabs(xs[i]-sol[i]) > tol)
            {
                std::cerr << "x[" << i << "]="
                     << sol[i] << " expected, " << xs[i]
                     << " achieved." << std::endl;
                std::cerr << endl;
                throw;
            }
        }

    std::cout <<  "obj=." <<x.obj() << std::endl;

    std::cout <<"ALMtest Markov SDDP..." << std::endl;

    mpmcvarequivalent<almproblem> ep(mp);


//    msddpsolution<mpmcvarequivalent<almproblem>,pdt,lastmdxi<vector<double>>,O> sx(ep,pd);
    msddpsolution<almproblem,pdt,lastmdxi<vector<double>>,O> sx(mp,pd);

    if(fabs(sx.obj().lb()-x.obj()) > 0.1)
    {
        std::cerr << "malmtest: opt="
             << x.obj() << " expected, " <<
             sx.obj().lb() << "<" << sx.obj().ubm() <<
             " achieved." << std::endl;
        throw;
    }

    std::cout << "almtest: lb=" <<
         sx.obj().lb() << " < " << "opt=" << x.obj()
              << " < ub=" << sx.obj().ubm() << std::endl;
    std::cout <<  "passed." << std::endl;
}

#endif // ALMTEST_H
