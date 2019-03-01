#ifndef ALMTEST_H
#define ALMTEST_H

#include "almproblem.h"
#include "mspptest.h"
#include "mspp/sddp.h"
#include "mspp/de.h"


class almdistribution:
        public vdistribution<double, unsigned int>
    ,
        public fdistribution<vector<double>, unsigned int>
{
public:
    using I_t = vector<double>;
    using C_t = unsigned int ;
    almdistribution(unsigned int n, const vector<double>& sc) :
         fn(n), fsc(sc) {}
private:
    virtual atom<vector<double>> atom_is(unsigned int i, const unsigned int& c) const
    {
        assert(c < fsc.size());
//std::cout << "c=" << c << " i=" << i << " x="
//         <<  fsc[c] - 0.5 +  (double) (2*i+1) / (double) (2*fn) << std::endl;
        return { {fsc[c] - 0.5 +  (double) (2*i+1) / (double) (2*fn)}, 1.0/ (double) fn };
    }

    virtual unsigned int natoms_is(const unsigned int& ) const
    { return fn; }

    unsigned int fn;
    vector<double> fsc;
    virtual unsigned int dim_is() const { return 1; }
};

template <typename O>
void alm1test(bool equivalent)
{
    almproblem<false> mp(0.5,0.05,1);
    lvdistribution<double>
        d({{0.5500,0.9},{0.8500,0.8}, {1.0500,1.1},{1.3500, 1.2}});
    using dist = fdprocessdistribution<lvdistribution<double>,noxi<vector<double>>>;
    dist pd({1,0},d,1);

    if(!equivalent)
    {
      std::cout <<"MALMtest DE direct..." << std::endl;

      desolution<almproblem<false>,dist,lastxi<vector<double>>,O> x(mp,pd);

      almproblem<false,nestedmcvar> nmp(0.5,0.05,1);


      sddpsolution<almproblem<false,nestedmcvar>,dist,lastxi<vector<double>>,O> sx(nmp,pd,{2});

      std::cerr << "alm1test: lb="
           <<
           sx.obj().lb() << "< opt=" << x.obj() << " <  " << sx.obj().ubm() <<
           " achieved." << std::endl;


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

          mpmcvarequivalent<almproblem<false>> ep(mp);

          desolution<mpmcvarequivalent<almproblem<false>>,dist,
              lastxi<vector<double>>,O> x(ep,pd);

          sddpsolution<mpmcvarequivalent<almproblem<false>>,dist,lastxi<vector<double>>,O> sx(ep,pd,{2});
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
    assert(T<=4);
    const double tol = 1e-5;

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

    almproblem<true> mp(0.5,0.05,T);

    std::cout <<"MALMtest DE..." << std::endl;

    desolution<almproblem<true>,pdt,hmczeta<vector<double>>,O> x(mp,pd);

    std::cout <<  "obj=." <<x.obj() << std::endl;

    vectors<unsigned int> counts;
    vectors<double> aves;
    x.x()->stats(counts, aves);

    std::cout <<"ALMtest Markov SDDP..." << std::endl;

    mpmcvarequivalent<almproblem<true>> ep(mp);


//    msddpsolution<mpmcvarequivalent<almproblem<true>>,pdt,lastmdxi<vector<double>>,O> sx(ep,pd);
    msddpsolution<almproblem<true>,pdt,hmczeta<vector<double>>,O> sx(mp,pd,vector<unsigned int>(T,1));

    cout << "Solution comparison" << endl;
    assert(aves.size()==sx.x()->vars.size());
    for(unsigned int k=0; k<aves.size(); k++)
    {
        assert(aves[k].size()==sx.x()->vars[k].size());
        cout << "Stage " << k << endl;
        for(unsigned int i=0; i<aves[k].size(); i++)
        {
          cout << aves[k][i] << "="  << sx.x()->vars[k][i]
             << "(" << sx.x()->sterrs[k][i] << ")" << endl;
        }
    }

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
