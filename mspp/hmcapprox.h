#ifndef HMCAPPROX_H
#define HMCAPPROX_H


#include "mspp/commons.h"
#include "mspp/process.h"
#include "mspp/sddp.h"
#include "mspp/cdists.h"

using namespace mspp;

template<typename I>
class covering : public object
{
public:
    using I_t = I;

    virtual unsigned int k() const = 0;
    virtual unsigned int i(const I& x) const = 0;
    virtual I e(unsigned int i) const = 0;
};

/// \brief General Monte Carlo hidden Markov chain approximation
/// \tparam P process to be approximated (descendant of \ref markovprocessdistribution)
/// \tparam C covering (descentat of \ref covering)
template<typename P, typename C>
class mchmcapproximation:
        public hmcprocessdistribution<mmcdistribution,
              altldistribution<typename P::D_t::I_t>>
{
    using D_t = altldistribution<typename P::D_t::I_t>;
    typename hmcprocessdistribution<mmcdistribution, D_t>::init
       makedists(const P& pd, const vector<C>& c, unsigned int numsc)
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
     mchmcapproximation(const P& pd, const vector<C>& c, unsigned int numsc) :
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
/// \tparam X distribution of \p xi
/// \tparam E distribution of \p eta
///
template<typename Y, typename E>
using xietadist = mijdistribution<Y,E,
       nomapping<pair<unsigned int, typename Y::I_t>>>;

///
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


/// conditional distribution given a region
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

    C fc;
    D fd;
};

/// \brief Hidden Markov chain approximation of a 1d AR process
/// \tparam P process to be approximated (descendant of \ref mmarkovprocessdistribution)
/// \tparam C covering (descentat of \ref onecovering)
/// \tparam R region approximating distribution
///
/// Both the marginal and conditional distributions defining P
/// have to be descendants \ref qdistribution<double>

template<typename P, typename C, typename R>
class hmcapproximation:
        public hmcprocessdistribution<mmcdistribution,R>
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


   static typename hmcprocessdistribution
              <mmcdistribution,R>::init
             makedists(const P& pd, const vector<C>& c)
    {
        assert(pd.dim()==c.size()+1);
        unsigned int d=c.size();
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
                for(unsigned int j=0; j<cols; j++)
                {
                    double e=t==0 ? xi0 : c[t-1].e(i);
                    const cdfdistribution<double>& cd = pd.d(i+1);
                    pair<double,double> region = c[t].region(j);
                    probability pij =
                       cd.cdf(region.second,e)-cd.cdf(region.first,e);
                    p[t][i][j] = pij;
                }
            }
            rows = cols;
        }

        vector<mmcdistribution> m;
        vector<R> omega;
        for(unsigned int t=0; t<d; t++)
        {
            m.push_back(mmcdistribution(p[t]));

            R o(pd.md(t+1),c[t]);
            omega.push_back(o);
        }
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
    }
public:
    /// \p c is indexed from zero, i.e. covering of time \p t=1 is indexed as \p 0
     hmcapproximation(const P& pd, const vector<C>& c) :
       hmcprocessdistribution<mmcdistribution, R>(makedists(pd,c))
     {
         check();
     }
     hmcapproximation(const P& pd, const vector<unsigned int>& ns) :
       hmcprocessdistribution<mmcdistribution, R>
          (makedists(pd,makec(pd,ns)))
     {
         check();
     }
};


template<typename P, typename C>
using chmcapproximation=hmcapproximation<P,C,
  regconddistribution<typename P::M_t,C>>;

/// not tested yet
/*template<typename P>
class epconthmcapproximation: public chmcapproximation<P,onedcovering>
{
public:
    epconthmcapproximation(const P& pd, const vector<unsigned int>& ns) :
        chmcapproximation<P,onedcovering>(pd, makec(pd, ns)) {}
};*/




#endif // HMCAPPROX_H
