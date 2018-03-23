#ifndef CPLEX_H
#define CPLEX_H

#include "mspp/solver.h"

namespace mspp
{

/// \addtogroup solvers Solvers
/// @{




class cplexlpsolver : public lpsolver
{
public:
    virtual void solve(const vardefs<realvar>& vars,
            const std::vector<double>& f,
            const std::vector<sparselinearconstraint_ptr>& constraints,
            const std::vector<std::string>& varnames,
            std::vector<double>& sol,
            double& objvalue) const;
};

/// @}

} // namespace


#endif // CPLEX_H
