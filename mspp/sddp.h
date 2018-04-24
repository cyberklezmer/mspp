#ifndef SDDP_H
#define SDDP_H

#include "mspp/msproblem.h"

namespace mspp
{

/// \addtogroup sddp SDDP
/// \ingroup sms
/// @{

template <typename V, typename F, typename R, typename X>
using swimsproblem
   =msproblem<V, F,interstagelinearconstraint,R,
               emptycondition<X>>;

template <typename X>
using sddppdistribution=processdistribution<imcdistribution<X>>;

using sddpsolution = std::vector<variable>;

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
        static_assert(
            std::is_same<typename P::G_t,interstagelinearconstraint>::value);
        static_assert(
            std::is_same<typename P::C_t,emptycondition<typename P::C_t::X_t>>::value);
        if constexpr (std::is_same<typename P::R_t,nestedmcvar>::value)
           throw("do it");
        else
           throw("emulate nested mcvar");
    }
};


/// @}


} // namespace

#endif // SDDP_H
