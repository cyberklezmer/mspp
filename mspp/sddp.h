#ifndef SDDP_H
#define SDDP_H

#include "mspp/msproblem.h"

namespace mspp
{

/// \addtogroup sddp SDDP
/// \ingroup sms
/// @{

template <typename P, typename Z>
using sddpsolution = variables<typename P::V_t>;

template <typename X=double>
class indexer
{
public:
    virtual unsigned int index(const subvectors<X>& s) const = 0;
};

template <typename X=double>
class trivialindexer : public indexer<X>
{
public:
    virtual unsigned int index(const subvectors<X>& s) const
    {
        return 0;
    }
};


template <typename P, typename Z, typename O, typename I=trivialindexer<typename P::X_t>>
class sddpmethod : public object
{
public:
    bool solve(
             const P& p,
             const Z& z,
             double& optimal,
             sddpsolution<P,Z>& sol)
    {
// TO BE DONE
    }
};


/// @}


} // namespace

#endif // SDDP_H
