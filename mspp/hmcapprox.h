#ifndef HMCAPPROX_H
#define HMCAPPROX_H


#include "mspp/commons.h"
#include "mspp/process.h"
#include "mspp/sddp.h"
#include "mspp/cdists.h"

using namespace mspp;

/// Base abstract class for a finite covering of any space
/// \tparam I space element
template<typename I>
class covering : public object
{
public:
    using I_t = I;

    virtual unsigned int k() const = 0;
    virtual unsigned int i(const I& x) const = 0;
    virtual I e(unsigned int i) const = 0;
};

template <typename C>
using coverings = vector<C>;

/// \brief General Monte Carlo hidden Markov chain approximation
/// \tparam P process to be approximated (descendant of \ref markovprocessdistribution)
/// \tparam C covering (descentat of \ref covering)
/// !!! not tested yet
template<typename P, typename C>
class mchmcapproximation:
        public hmcprocessdistribution<mmcdistribution,
              altldistribution<typename P::D_t::I_t>>
{
    using D_t = altldistribution<typename P::D_t::I_t>;
    typename hmcprocessdistribution<mmcdistribution, D_t>::init
       makedists(const P& pd, const coverings<C>& c, unsigned int numsc)
    {
        using I_t = typename P::D_t::I_t;
        assert(pd.dim()==c.size()+1);
        unsigned int d=c.size();
        using pmatrix = vector<vector<probability>>;
        vector<pmatrix> p;
        const probability dp = 1.0 / static_cast<double>(numsc);

        vector<vector<vector<I_t>>> omegas;

        unsigned int rows = 1;
        for(unsigned int t=0; t<d; t++)
        {
            unsigned int cols = c[t].k();
            for(unsigned int j=0; j<rows; j++)
                p[t].push_back(vector<probability>(cols,0.0));
            omegas[t].push_back(vector<I_t>(cols));
            rows = cols;
        }

        for(unsigned int i=0; i<numsc; i++)
        {
            scenario<I_t> s = pd.draw();
            assert(s[0]==pd->x0());
            unsigned int oldi = 0;
            for(unsigned int t=0; t<d; t++)
            {
                I_t x = s[t+1];
                unsigned int newi=c[t].i();
                assert(oldi<p[t].size());
                assert(newi<p[t][i].size());
                p[t][oldi][newi] += dp;
                assert(newi < omegas[t].size());
                omegas[t][newi].push_back(s);
            }
        }
        vector<mmcdistribution> m;
        vector<D_t> xi;
        for(unsigned int t=0; t<d; t++)
        {
            m.push_back(mmcdistribution(p[t]));
            xi.push_back(D_t(omegas[t]));
        }
        return {pd.xi0(),m,xi};
    }
public:
    /// \p c is indexed from zero, i.e. covering of time \p t=1 is indexed as \p 0
     mchmcapproximation(const P& pd, coverings<C>& c, unsigned int numsc) :
       hmcprocessdistribution<mmcdistribution, D_t>(makedists(pd,c,numsc))
     {
         static_assert(std::is_base_of<markovprocessdistribution<typename P::D_t>,P>::value);
         static_assert(std::is_base_of<covering<typename C::I_t>,C>::value);
         static_assert(
             std::is_base_of<mdistribution<typename P::D_t::I_t,typename P::D_t::C_t>,
                        typename P::D_t>::value);
     }
};


///
/// \tparam Y distribution of \p xi
/// \tparam E distribution of \p eta
///
template<typename Y, typename E>
using xietadist = mijdistribution<Y,E,
       nomapping<pair<unsigned int, typename Y::I_t>>>;

/// Joint process of a hiden Markov and an I.i.d. processes
/// \tparam Y processdistribution of \p xi
/// \tparam E processdistribution of \p eta
///
template<typename Y, typename E>
class xietaprocessdist :
    public hmcprocessdistribution<typename Y::M_t,
         xietadist<typename Y::Y_t,typename E::D_t>>
{
    using Y_t = xietadist<typename Y::Y_t,typename E::D_t>;
    using M_t = typename Y::M_t;
    typename hmcprocessdistribution<M_t, Y_t>::init
       makedists(const Y& y, const E& e)
    {
        assert(y.dim()==e.dim());

        vector<M_t> m;
        vector<Y_t> xe;
        for(unsigned int i=1; i<y.dim(); i++)
        {
            m.push_back(y.d(i).first());
            xe.push_back(Y_t(y.d(i).second(),e.d(i)));
        }
        //typename Y::Y_t::I_t y0 = y.e().x().second;
        //typename E::D_t::I_t e0 = e.e().x();
        //pair<typename Y::Y_t::I_t,typename E::D_t::I_t> in(y0,e0);
        //return {in,m,xe};
        return {{y.e().x().second,e.e().x()},m,xe};
    }
public:
    xietaprocessdist(const Y& y, const E& e)
        : hmcprocessdistribution<M_t,Y_t>
            (makedists(y,e))
    {
        static_assert(std::is_base_of
            <hmcprocessdistribution<M_t,typename Y::Y_t>,Y>::value);
        static_assert(std::is_base_of
            <iidprocessdistribution<typename E::D_t>,E>::value);
    }

};


template<typename Y, typename E>
using fdxietabase
  = treedistribution<pair<unsigned int,pair<typename Y::D_t::S_t::I_t,typename E::D_t::I_t>>>;

/// Discrete version of \ref xietaprocessdist
/// \tparam Y, \tparam E see \ref xietaprocessdist
///
template<typename Y, typename E>
class dxietaprocessdist : public xietaprocessdist<Y,E>,
    public fdxietabase<Y,E>
{
public:
    using MD_t = typename Y::D_t::F_t;
    using XiD_t = typename Y::D_t::S_t;
    using EtaD_t = typename E::D_t;
    using X_t = pair<unsigned int,pair<typename Y::D_t::S_t::I_t,typename E::D_t::I_t>>;
    dxietaprocessdist(const Y& y, const E& e)
        : xietaprocessdist<Y,E>(y, e),
          treedistribution<X_t>(y.dim()-1)
    {
        static_assert(std::is_base_of
            <fdistribution<typename XiD_t::I_t, unsigned int,XiD_t::flistdef>,XiD_t>::value);
        static_assert(std::is_base_of
            <fdistribution<typename EtaD_t::I_t, nothing,EtaD_t::flistdef>,EtaD_t>::value);
    }
private:
    virtual void branches_are(const vector<atom<X_t>>& e,
         vector<atom<X_t>>& es) const
    {
        static_assert(std::is_same<pair<unsigned int,pair<
                      typename XiD_t::I_t, typename EtaD_t::I_t>>,X_t>::value);
        unsigned int k=e.size();
        if(k==0)
        {
            auto a = this->e().x();
            es.push_back({a,1.0});
        }
        else
        {
            assert(this->T());
            scenario<X_t> s;
            for(unsigned int i=0; i<k; i++)
                s.push_back(e[i].x);

            const MD_t& md = this->d(k).first();
            const XiD_t& xid = this->d(k).second().first();
            const EtaD_t& etad = this->d(k).second().second();

            vector<atom<typename EtaD_t::I_t>> ea;
            etad.atoms(ea);
            unsigned int cond = e[k-1].x.first;
            assert(cond==this->z(s));
            for(unsigned int m=0; m<md.nstates(); m++)
            {
                probability tp = md.transprob(cond,m);             ;
                vector<atom<typename XiD_t::I_t>> ya;
                xid.atoms(ya,m);

                for(unsigned int i=0; i<ya.size(); i++)
                    for(unsigned int j=0; j<ea.size(); j++)
                    {
                        es.push_back({{m,{ya[i].x,ea[j].x}},ya[i].p*ea[j].p*tp});
                    }
            }
        }
    }
    virtual unsigned int dim_is() const
       { return this->T()+1; }

};

/// \ref covering of the real line
class onedcovering : public covering<double>
{
public:
    /// the thresholds have to be ordered
    onedcovering(const vector<double>& ths, const vector<double>& es) :
      fths(ths), fes(es)
    {
        assert(fths.size()+1 == fes.size());
        for(unsigned int i=1; i<fths.size(); i++)
            assert(fths[i-1]<fths[i]);
        for(unsigned int i=1; i<fes.size(); i++)
            assert(fes[i-1]<fes[i]);
        for(unsigned int i=1; i<fes.size(); i++)
            assert(fes[i]>=fths[i-1]);
        for(unsigned int i=0; i<fes.size()-1; i++)
            assert(fes[i]<fths[i]);
    }
    virtual unsigned int k() const
    {
        return fths.size()+1;
    }
    virtual unsigned int i(const double& x) const
    {
        unsigned int i;
        for(i=0; i<fths.size(); i++)
        {
            if(x<fths[i])
                return i;
        }
        return i;
    }
    virtual double e(unsigned int i) const
    {
        assert(i<fes.size());
        return fes[i];
    }
    pair<double,double> region(unsigned int i) const
    {
        unsigned int k = fes.size();
        assert(i<k);
        double lb = i==0 ? mspp::min<double>(): fths[i-1];
        double hb = i==k-1 ? mspp::max<double>() : fths[i];
        return { lb,hb };
    }
private:
    vector<double> fths;
    vector<double> fes;
};


/// \ref Equiprobable coverings - regions are equiprobable intervals
/// with respect to the marginal distribution
/// \tparam P the proces defining the coverings
template <typename P>
class epcoverings: public coverings<onedcovering>, public object
{
    static vector<onedcovering> makec(const P& pd, const vector<unsigned int>& ns)
    {
        unsigned int d=pd.dim()-1;
        assert(ns.size()==d);
        vector<onedcovering> ret;
        for(unsigned int i=0; i<d; i++)
        {
            const typename P::M_t& distr = pd.md(i+1);
            unsigned int n=ns[i];
            assert(n);
            double p = 1.0 / static_cast<double>(n);
            double cp=0.0;
            vector<double> es, ths;

            for(unsigned int j=0;;)
            {
                double e=distr.quantile(cp+p/2.0);
                es.push_back(e);
                cp+=p;
                if(++j==n)
                    break;
                double h=distr.quantile(cp);
                ths.push_back(h);
            }
            ret.push_back(onedcovering(ths,es));
        }
        return ret;
    }
public:
    epcoverings(const P& p, const vector<unsigned int>& ns)
        : vector<onedcovering>(makec(p,ns))
    {
        static_assert(std::is_base_of
            <mmarkovprocessdistribution<typename P::D_t,typename P::M_t>,P>::value);
        static_assert(
            std::is_base_of<qdistribution<nothing>,
                       typename P::M_t>::value);
        static_assert(
            std::is_base_of<qdistribution<double>,
                       typename P::D_t>::value);
    }
};




/// conditional distribution given a region
/// \tparam D the original distribution
/// \tparam C covering
///
template <typename D, typename C>
class regconddistribution : public qdistribution<unsigned int>
{
public:
    regconddistribution(const D& d, const C& c)
        : fd(d), fc(c)
    {
        static_assert(std::is_base_of<onedcovering,C>::value);
        static_assert(std::is_base_of<qdistribution<nothing>,D>::value);
    }

    /// returns the source distribution
    const D& srcd() const { return fd; }
    const C& c() const { return fc; }
private:
    virtual double quantile_is(probability p, const unsigned int& i) const
    {
        assert(i<fc.k());
        pair<double, double> r = fc.region(i);
        double arg = p * fd.cdf(r.second) + (1-p) * fd.cdf(r.first);
        assert(arg);
        return fd.quantile(arg);
    }
    virtual probability cdf_is(double x, const unsigned int& i) const
    {
        assert(i<fc.k());
        pair<double, double> r = fc.region(i);
        probability cl = fd.cdf(r.first);
        double denom=fd.cdf(r.second)-cl;
        assert(denom>0);
        return (fd.cdf(x)-cl)/denom;
    }
    virtual double mean_is(const unsigned int& i) const
    {
        if constexpr(std::is_base_of<normaldistribution,D>::value)
        {
            stdnormaldistribution nd;
            pair<double, double> r = fc.region(i);
            double a1 = (r.first - fd.m()) / fd.sd();
            double a2 = (r.second - fd.m()) / fd.sd();
            double phi1 = pow(2*M_PI,-0.5)*exp( -a1*a1 / 2);
            double phi2 = pow(2*M_PI,-0.5)*exp(- a2*a2 / 2);
            return fd.m() - fd.sd()
                    * (phi2 - phi1) / (nd.cdf(a2)-nd.cdf(a1));
        }
        else
            return realdistribution<unsigned int>::mean_is(i);
    }


    C fc;
    D fd;
};

/// \brief \ref covering defined Hidden Markov chain approximation
/// of a 1D Markov process
/// \tparam P process to be approximated (descendant of \ref mmarkovprocessdistribution)
/// \tparam C covering (descentat of \ref onedcovering)
/// \tparam R distribution approximating region, has to have a ctr <tt>
/// R(const D& d, const C& c)</tt> where \p D is a \ref qdistribution and
/// \p D is a \ref onedcovering
///
/// Both the marginal and conditional distributions defining \p P
/// have to be descendants \ref qdistribution<double>
///

template<typename P, typename C, typename R, bool martingale, typename M=idmapping<double>>
class cdhmcapproximation:
        public hmcprocessdistribution<mmcdistribution,R>
{
    static constexpr bool mid = std::is_same<M,idmapping<double>>::value;
public:
    static void makeEe(const vector<double>& es, vector<probability>& pr, double e)
    {
        #define err "Cannot make the process of representatives martingale: "
        #define cherr " Perhaps make a better algorihtm and push "
        assert(es.size() == pr.size());
        unsigned int dim = es.size();
        if(dim==1)
        {
            if(fabs(es[0]-e) > probabilitytolerance)
                throw mspp::exception(err  "trivial covering. ");
            else
                return;
        }
        unsigned int firstright = 0;
        probability sum = 0;
        for(; firstright<dim; firstright++)
        {
            sum+=pr[firstright];
            if(es[firstright] > e)
                break;
        }
        if(firstright == 0 || firstright == dim)
            throw mspp::exception(err "all successors left or right." );
        double E=0;

        double left = 0;
        double Left = 0;
        unsigned int j=0;
        for(; j<firstright; j++)
        {
           E += pr[j]*es[j];
           Left += pr[j]*es[j];
           left += pr[j];
        }
        double right = 0;
        double Right = 0;
        for(; j<dim; j++)
        {
           E+= pr[j]*es[j];
           Right += pr[j]*es[j];
           right += pr[j];
        }
        assert(right > 0.0);
        assert(left > 0.0);
        double c = (e-Right/right) / (Left - left/right*Right);
        double d = (1.0 - c* left) / right;
        assert(c>=0);
        assert(d>=0);
        j=0;
        for(; j<firstright; j++)
        {
            pr[j] *= c;
            if(pr[j] > 1.0)
                throw mspp::exception(err cherr);
        }
        for(; j<dim; j++)
        {
            pr[j] *= d;
            if(pr[j] > 1.0)
                throw mspp::exception(err cherr);
        }
        #undef cherr
        #undef err
    }

private:
   static typename hmcprocessdistribution<mmcdistribution,R>::init
             makedists(const P& pd, const vector<C>& c,vector<double> trend)
    {
        assert(pd.dim()==c.size()+1);
        assert(trend.size()==0 || martingale);
        unsigned int d=c.size();

        vector<R> omega;
        for(unsigned int t=0; t<d; t++)
        {
            R o(pd.md(t+1),c[t]);
            omega.push_back(o);
        }

        using pmatrix = vector<vector<probability>>;
        vector<pmatrix> p(d);
        double xi0 = pd.e().x();

        unsigned int rows = 1;
        for(unsigned int t=0; t<d; t++)
        {
            unsigned int cols = c[t].k();
            for(unsigned int i=0; i<rows; i++)
            {
                p[t].push_back(vector<probability>(cols));
                double e=t==0 ? xi0 : c[t-1].e(i);
                for(unsigned int j=0; j<cols; j++)
                {
                    const cdfdistribution<double>& cd = pd.d(t+1);
                    pair<double,double> region = c[t].region(j);
                    probability pij =
                       cd.cdf(region.second,e)-cd.cdf(region.first,e);
                    p[t][i][j] = pij;
                }
                if constexpr(martingale)
                {
                    vector<probability>& pr = p[t][i];
                    vector<double> es;
                    for(unsigned int j=0; j<cols; j++)
                    {
                        if constexpr(mid)
                        {
                            es.push_back(omega[t].mean(j));
                        }
                        else
                        {
                            vector<atom<typename R::I_t>> vcs;
                            omega[t].atoms(vcs,j);
                            double sum=0;
                            for(unsigned int i=0; i<vcs.size(); i++)
                                sum += vcs[i].p * M()(vcs[i].x);
                            es.push_back(sum);
                        }
                    }
                    double tr = t < trend.size() ? trend[t] : 0.0;
                    makeEe(es,pr,M()(e)+tr);
//cout << t << ","
//     << (M()(e)+tr) << endl;
//for(unsigned int i=0; i<es.size(); i++)
//    cout << i << "," << pr[i] << "," << es[i] << endl;

                }
            }
            rows = cols;
        }

        vector<mmcdistribution> m;
        for(unsigned int t=0; t<d; t++)
            m.push_back(mmcdistribution(p[t]));
        return {xi0,m,omega};
    }
    void constexpr check() const
    {
        static_assert(std::is_base_of
            <mmarkovprocessdistribution<typename P::D_t,typename P::M_t>,P>::value);
        static_assert(std::is_base_of<onedcovering,C>::value);
        static_assert(
            std::is_base_of<qdistribution<nothing>,
                       typename P::M_t>::value);
        static_assert(
            std::is_base_of<qdistribution<double>,
                       typename P::D_t>::value);
        static_assert(
            std::is_base_of<mapping<double,double>,M>::value);
    }
public:
    /// \p c is indexed from zero, i.e. covering of time \p t=1 is indexed as \p 0
     cdhmcapproximation(const P& pd, const vector<C>& c,vector<double> trend = vector<double>(0)) :
       hmcprocessdistribution<mmcdistribution, R>(makedists(pd,c,trend))
     {
         check();
     }
     cdhmcapproximation(const P& pd, const vector<unsigned int>& ns, vector<double> trend = vector<double>(0)) :
       hmcprocessdistribution<mmcdistribution, R>
          (makedists(pd,epcoverings(pd,ns),trend))
     {
         check();
     }
};


template<typename P, typename C, bool martingale, typename M=idmapping<double>>
using chmcapproximation=cdhmcapproximation<P,C,
  regconddistribution<typename P::M_t,C>, martingale,M>;


/// Discretized conditional distribution given a region
/// \tparam D original distribution
/// \tparam C covering
/// \tparam N number of atoms
///
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

/// Discrete version of \ref hmcapproximation
/// \tparam P the process approximated
/// \tparam C covering
/// \tparam N number of atoms in each region
/// \tparam M mappaing according to which martingale is callibrated
template<typename P, typename C, unsigned int N, bool martingale, typename M=idmapping<double>>
class dhmcapproximation:
    public cdhmcapproximation<P,C,drcdistribution<typename P::M_t,C,N>, martingale, M>,
    public treedistribution<pair<unsigned int,double>>
{
public:
  // tohle je malinko prasárna (měl by se v těch templatech udělat
  // jednoznačnej pořádek
  using X_t = pair<unsigned int,double>;
  using R_t = drcdistribution<typename P::M_t,C,N>;
  dhmcapproximation(const P& pd, const vector<C>& c,
                    vector<double> trend = vector<double>(0)) :
    cdhmcapproximation<P,C,
      drcdistribution<typename P::M_t,C,N>,martingale>(pd,c,trend),
      treedistribution<pair<unsigned int,double>>(pd.dim()-1) {}

  // twice calling constructor of epcoverings. maybe covering should be
  // strored in hmcapproximation
  dhmcapproximation(const P& pd, const vector<unsigned int>& ns,
                    vector<double> trend = vector<double>(0)) :
    cdhmcapproximation<P,C,
      drcdistribution<typename P::M_t,onedcovering,N>,martingale,M>(pd,epcoverings(pd,ns),trend),
      treedistribution<pair<unsigned int,double>>(pd.dim()-1)//,
//    fc(epcoverings(pd,ns)), fpd(pd)
  {}

private:
  virtual void branches_are(const vector<atom<X_t>>& e,
                               vector<atom<X_t>>& es) const
  {
      unsigned int t=e.size();
      if(!t)
      {
          double a = this->e().x().second;
          es.push_back({{0,a},1.0});
      }
      else
      {
         unsigned int cond = e[t-1].x.first;

         vectors<probability> p = this->d(t).first().m();
         assert(cond < p.size());
         unsigned int k=p[cond].size();
         for(unsigned int i=0; i<k; i++)
         {
              probability rp = p[cond][i];
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
//  vector<C> fc;
//  P fpd;
  virtual unsigned int dim_is() const { return
        cdhmcapproximation<P,C,drcdistribution<typename P::M_t,C,N>, martingale,M>::dim_is(); }
};



/// not tested yet
/*template<typename P>
class epconthmcapproximation: public chmcapproximation<P,onedcovering>
{
public:
    epconthmcapproximation(const P& pd, const vector<unsigned int>& ns) :
        chmcapproximation<P,onedcovering>(pd, makec(pd, ns)) {}
};*/




#endif // HMCAPPROX_H
