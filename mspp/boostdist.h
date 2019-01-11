#ifndef BOOSTDIST_H
#define BOOSTDIST_H

#include "mspp/random.h"
#include <boost/math/distributions.hpp>

namespace mspp
{

/// \addtogroup realdists Real distributions
/// \ingroup Distributions
/// @{

template <typename B>
class boostdistribution : virtual public qdistribution<nothing>
{
public:
    boostdistribution() {}
    boostdistribution(const B& d) : fd(d) {}
    virtual probability cdf(double x, const nothing& ) const
    {
        return boost::math::cdf(fd,x);
    }
    virtual double quantile(probability p, const nothing& ) const
    {
        return boost::math::quantile(fd,p);
    }
private:
    B fd;
    virtual double do_draw(const nothing&) const
    {
        return fd();
    }
};

using stdnormaldistribution
  =boostdistribution<boost::math::normal_distribution<double>>;

class normaldistribution : public

/// @}

} // namespace

#endif // BOOSTDIST_H

