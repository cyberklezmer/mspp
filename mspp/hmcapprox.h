#ifndef HMCAPPROX_H
#define HMCAPPROX_H

#include "mspp/commons.h"
#include "mspp/process.h"
#include "mspp/sddp.h"

using namespace mspp;

template<typename I>
class covering
{
public:
    using I_t = I;

    virtual unsigned int dim() const = 0;
    virtual unsigned int k() const = 0;
    virtual unsigned int i(const I& x) const = 0;
    virtual I e(unsigned int i) const = 0;
};

/// \brief Hidden MC approximation
/// \tparam P process to be approximated (descendant of \ref markovprocessdistribution)
/// \tparam C covering (descentat of \ref covering)
template<typename P, typename C>
class hmcapproximation:
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
     hmcapproximation(const P& pd, const vector<C>& c, unsigned int numsc) :
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
using xietadist = ijdistribution<Y,E,
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

class onedcovering : public covering<double>
{
public:
    /// the thresholds have to be ordered
    onedcovering(vector<double>& ths, vector<double>& es) :
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
private:
    vector<double> fths;
    vector<double> fes;
};

#endif // HMCAPPROX_H
