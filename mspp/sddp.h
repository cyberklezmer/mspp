#ifndef SDDP_H
#define SDDP_H


namespace mspp
{

/// \addtogroup sddp SDDP
/// \ingroup sms
/// @{

template <typename V, typename F, typename X>
using swimsproblem
   =msproblem<V, F,interstagelinearconstraint,
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
        throw("to be done");
    }
};


/// @}


} // namespace

#endif // SDDP_H
