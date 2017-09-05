#ifndef LINEAR_H
#define LINEAR_H

#include "mspp/commons.h"

namespace mspp
{

class linearfunction: public object
{
public:
    linearfunction(const std::vector<double>& acoefs)
       : coefs(acoefs) {}

    linearfunction(unsigned int dim=0)
             : coefs(dim,0.0) {}
    std::vector<double> coefs;
    unsigned int dim() const
       { return coefs.size(); }
};


using linearfunction_ptr=std::shared_ptr<linearfunction>;

class linearconstraint
{
public:
    enum type {eq, geq, leq};
    linearconstraint(const std::vector<double>& alhs,
                       type at, double arhs) :
         lhs(alhs), t(at), rhs(arhs) {}

    linearconstraint(unsigned int lhsdim, type at = eq, double arhs=0.0 )
        : lhs(lhsdim,0), t(at), rhs(arhs) {}
    std::vector<double> lhs;

    type t;
    double rhs;
    unsigned int dim() { return lhs.size(); }
};

using linearconstraint_ptr = std::shared_ptr<linearconstraint>;

using linearconstraint_list = std::vector<linearconstraint>;
using linearconstraint_list_ptr = std::shared_ptr<linearconstraint_list>;

template<class Xi>
using linearproblem = problem<linearfunction,linearconstraint,Xi>;

template<class Xi>
using linearproblem_ptr  = std::shared_ptr<linearproblem<Xi>>;

}

#endif // LINEAR_H
