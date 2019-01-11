#ifndef CDISTS_H
#define CDISTS_H

#include "mspp/random.h"
#include "mspp/process.h"
#include <boost/math/distributions.hpp>

namespace mspp
{

/// \addtogroup realdists Real distributions
/// \ingroup Distributions
/// @{

template <typename B>
class boostdistribution : virtual public qdistribution<novalue>
{
public:
    boostdistribution() {}
    boostdistribution(const B& d) : fd(d) {}
private:
    virtual probability cdf_is(double x, const novalue& ) const
    {
        return boost::math::cdf(fd,x);
    }
    virtual double quantile_is(probability p, const novalue& ) const
    {
        return boost::math::quantile(fd,p);
    }
    B fd;
// tbd
//      virtual double do_draw(const novalue&) const
//    {
//        generator(,fd)
//        return fd();
//    }
};

using stdnormaldistribution
  =boostdistribution<boost::math::normal_distribution<double>>;

using normaldistribution=scaleddistribution<stdnormaldistribution>;

class arnormalprocessdistribution :
        public mprocessdistribution
          <ardistribution<normaldistribution>, lastxi<double>, normaldistribution>
{
public:
    arnormalprocessdistribution
         (double xi0, double m, double sd, double ar, unsigned int T):
      mprocessdistribution<ardistribution<normaldistribution>,
                               lastxi<double>, normaldistribution>
          (xi0,
             ardistribution<normaldistribution>(normaldistribution(m,sd),ar),
             T),
      vdistribution<double,novalue>(T+1)
    { assert(T>=1); }

private:
    virtual normaldistribution md_is(unsigned int i) const
    {
        assert(i>0);
        const auto& ard=this->d(1);
        const auto& nd = ard.srcd();
        const auto& dd = this->e();
        double a = ard.a();
        double m = nd.m();
        double sd = nd.sd();
        double v = sd*sd;
        double mu, var;
        if(a==1)
        {
            mu = dd.x() + i*m;
            var = i*v;
        }
        else
        {
           mu =  dd.x()*pow(a,i) + m*(1.0-pow(a,i))/(1.0-a);
           var = v * (1.0-pow(a,2*i))/(1.0-a);
        }
        return normaldistribution(mu,sqrt(var));
    }
};

/// @} // realdists

/// \addtogroup realprocs Real distributions
/// \ingroup processes
/// @{

/// @} // realprocs

} // namespace


#endif // CDISTS_H
