#ifndef SDDP_H
#define SDDP_H

#include "mspp/msproblem.h"

namespace mspp
{

/// \addtogroup sddp SDDP
/// \ingroup sms
/// @{


using sddpsolution = variables;

template <typename C>
class indexer
{
public:
    virtual unsigned int index(const C& c) const = 0;
};

template <typename C>
class trivialindexer
{
public:
    virtual unsigned int index(const C&) const
    {
        static_assert(std::is_same<C, emptycondition<typename C::X_t>>::value);
        return 0;
    }
};


template <typename P, typename Z>
class sddpmethod : public object
{
public:
    bool solve(
             const P& p,
             const Z& z,
             const lpsolver& lps,
             double& optimal,
             sddpsolution& sol)
    {
//        if constexpr (std::is_same<typename P::R_t,nestedmcvar>::value)
//           throw("do it");
//        else
//           throw("emulate nested mcvar");
        static_assert( std::is_same<typename P::C_t,
                         emptycondition<typename P::X_t>>::value );
//        return solve(p,z,lps, optimal, sol, trivialindexer<typename P::C_t>);
    }
};


/// @}


} // namespace

#endif // SDDP_H
