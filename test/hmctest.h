#ifndef HMCTEST_H
#define HMCTEST_H

#include "almproblem.h"
#include "mspp/sddp.h"
#include "mspp/de.h"
#include "mspp/cdists.h"
#include "mspp/hmcapprox.h"

namespace mspp
{

class expvmapping: public mapping<vector<double>,vector<double>>
{
public:

  virtual vector<double> operator() (const vector<double>& a) const
  {
    vector<double> dst;
    for(unsigned int i=0; i<a.size(); i++)
       dst.push_back(exp(a[i]));
    return dst;
  }
};

class expmapping: public mapping<double,vector<double>>
{
public:

  virtual vector<double> operator() (const double& a) const
  {
    return { exp(a) };
  }
};

/// discretized conditional distribution given a region
template <typename D, typename C, unsigned int N>
class drcdistribution : public fdistribution<double,unsigned int, false>
{
public:
    static constexpr unsigned int fn = N;
    drcdistribution(const D& d, const C& c)
        : fd(d), fc(c)
    {
        static_assert(std::is_base_of<onedcovering,C>::value);
        static_assert(std::is_base_of<qdistribution<nothing>,D>::value);
    }

    /// returns the source distribution
    const D& srcd() const { return fd; }
    const C& c() const { return fc; }
private:
    virtual unsigned int natoms_is(const unsigned int& ) const
    {
        return N;
    }
    virtual atom<double> atom_is(unsigned int i, const unsigned int& c) const
    {
        assert(c<fc.k());
        pair<double, double> r = fc.region(c);
        probability l=fd.cdf(r.first);
        probability h=fd.cdf(r.second);
        assert(l<h);

        probability cp= static_cast<double> (2*i+1) / static_cast<double>(2*N);

        probability p=l+cp*(h-l);
        assert(p<1.0);
        return {fd.quantile(p) , 1.0 / static_cast<double> (N) };
    }

    C fc;
    D fd;
};

template<typename P, typename C, unsigned int N>
class dhmcapproximation:
    public hmcapproximation<P,C,drcdistribution<typename P::M_t,C,N>>,
    public treedistribution<pair<unsigned int,double>>
{
public:
  // tohle je malinko prasárna
  using X_t = pair<unsigned int,double>;
  using R_t = drcdistribution<typename P::M_t,C,N>;
  dhmcapproximation(const P& pd, const vector<C>& c) :
    hmcapproximation<P,C,
      drcdistribution<typename P::M_t,C,N>>(pd,c),
      treedistribution<pair<unsigned int,double>>(pd.dim()-1),
    fc(c), fpd(pd) {}

private:
  virtual void branches_are(const vector<atom<X_t>>& e,
                               vector<atom<X_t>>& es) const
  {
      unsigned int t=e.size();
      if(!t)
      {
          double a = this->fpd.e().x();
          es.push_back({{0,a},1.0});
      }
      else
      {
         assert(t);
         double rep;
         if(t>1)
         {
            unsigned int i=e[t-1].x.first;
            rep = fc[t-2].e(i);
         }
         else
            rep = fpd.e().x();
          assert(t-1<fc.size());
          unsigned int k=fc[t-1].k();
          for(unsigned int i=0; i<k; i++)
          {
              const typename P::D_t& srcd = fpd.d(t);
              pair<double,double> reg = fc[t-1].region(i);
              probability rp = srcd.cdf(reg.second,rep)-srcd.cdf(reg.first,rep);
              probability p = rp / static_cast<double>(N);
              for(unsigned int j=0; j<N; j++)
              {
                  const R_t& ad=this->d(t).second();
                  double nr = ad(j,i).x;
                  es.push_back({{i,nr},p});
              }
          }
        }

  }
  vector<C> fc;
  P fpd;
  virtual unsigned int dim_is() const { return fpd.dim(); }
};


template <typename O, bool stocheta=false>
void hmctest(unsigned int T, bool cont = false)
{
    assert(T>=1);
    assert(T<=2);

    double xim = -0.3;
    double xisd = 0.5;

    arnormalprocessdistribution xipd(0,xim,xisd,1.0,T);

    onedcovering c1({0+xim},{-0.67+xim, 0.67+xim});
    onedcovering c2({-0.2+2*xim,0.2+2*xim},{-0.33+2*xim,2*xim, 0.33+ 2*xim});

    vector<onedcovering> c({c1});
    if(T==2)
        c.push_back(c2);

    using dhat
       = dhmcapproximation<arnormalprocessdistribution,onedcovering,20>;

    dhat dha(xipd,c);

    if constexpr(stocheta)
    {
       double etam = -0.5;
       double etasd = 0.2;

       using apn = discretization<normaldistribution,nothing>;

       using etapd = iidprocessdistribution<apn>;

       constexpr double  minf = -std::numeric_limits<double>::infinity();

       etapd etap(minf,apn(normaldistribution(etam,etasd),10),T);

       using xietat = xietaprocessdist<dhat,etapd>;
//       xietat xieta(dhat,etap);
//       using zetat = hmczeta<pair<double, doube>,exppmapping>;

    }
    else
    {
      using zetat = hmczeta<double,expmapping>;
  static_assert(std::is_same<zetat::R_t,vector<double>>::value);
      using apt = almproblem<true, /* nestedmcvar */ mpmcvar>;
      apt ap(0.5,0.05,T);


      desolution<apt,dhat,zetat,O> dx(ap,dha);

      std::cout <<  "obj=." << dx.obj() << std::endl;

      if(!cont)
      {
          msddpsolution<apt,dhat, zetat,O> sx(ap,dha);

          std::cout << "hmctest disc: lb=" <<
               sx.obj().lb() << " < " << "opt=" << dx.obj()
                    << " < ubm=" << sx.obj().ubm() << ", ubb="
                    << sx.obj().ubb() << std::endl;

          if(fabs(sx.obj().lb()-dx.obj()) > 0.01)
          {
              std::cerr << "malmtest: opt="
                   << dx.obj() << " expected, " <<
                   sx.obj().lb() << "<" << sx.obj().ubm() <<
                   " achieved." << std::endl;
              throw;
          }
      }
      else
      {
          using chat
             = chmcapproximation<arnormalprocessdistribution,onedcovering>;

          chat cha(xipd,c);

          msddpsolution<apt,chat, zetat,O> cx(ap,cha);

          std::cout << "hmctest cont: lb=" <<
               cx.obj().lb() << " < " << "opt=" << dx.obj()
                    << " < ubm=" << cx.obj().ubm() << ", ubb="
                    << cx.obj().ubb() << std::endl;

          if(fabs(cx.obj().lb()-dx.obj()) > 0.1)
          {
              std::cerr << "malmtest: opt="
                   << dx.obj() << " expected, " <<
                   cx.obj().lb() << "<" << cx.obj().ubm() <<
                   " achieved." << std::endl;
              throw;
          }
      }
    } // stocheta
    std::cout <<  "passed." << std::endl;
}

} // namespace


#endif // HMCTEST_H
