#ifndef CPLEX_H
#define CPLEX_H

#include "mspp/lpsolver.h"

namespace mspp
{

class cplexlpsolver : public lpsolver
{
public:
    virtual void solve(const varrange_list& vars,
            const linearfunction& objective,
            const std::vector<sparselinearconstraint_ptr>& constraints,
            const std::vector<std::string>& varnames,
            std::vector<double>& sol,
            double& objvalue) const;
};

}


#endif // CPLEX_H
