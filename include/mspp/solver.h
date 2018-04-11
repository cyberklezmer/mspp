#ifndef SOLVER_H
#define SOLVER_H

#include "mspp/problem.h"

namespace mspp
{

/// \addtogroup solvers Solvers
/// @{

class sparselinearconstraint // tbd sjednotit s linearconstraintem
{
public:
//    enum type {eq, geq, leq};

    struct iitem
    {
        unsigned int index;
        double value;
    };

    sparselinearconstraint(
       linearconstraint::type at = linearconstraint::eq, double arhs=0.0 )
        : flhsdim(0), t(at), rhs(arhs) {}

    linearconstraint::type t;
    double rhs;
    unsigned int numnz() const { return fnzs.size(); }
    const iitem& nz(unsigned int i) const { return fnzs[i]; }
    unsigned int lhsdim() { return flhsdim; }

    std::shared_ptr<std::vector<double>> lhs() const
    {
        std::shared_ptr<std::vector<double>> r(new std::vector<double>(flhsdim,0));
        for(unsigned int i=0; i < fnzs.size(); i++)
        {
            const iitem& it = fnzs[i];
            (*r)[it.index] = it.value;
        }
        return r;
    }


    void set_lhs(unsigned int i, double v)
    {
        fnzs.push_back({i,v}); /*tbd sort*/
        if(i+1 > flhsdim)
            flhsdim = i+1;
    }
private:
    std::vector<iitem> fnzs;
    unsigned int flhsdim;
};

using sparselinearconstraint_ptr = std::shared_ptr<sparselinearconstraint>;


class lpsolver : public object
{
public:
    virtual void solve(const vardefs<realvar>& vrs,
            const std::vector<double>& f,
            const std::vector<sparselinearconstraint_ptr>& msconstraints, // tbd predef constraintslist
            const std::vector<std::string>& varnames,
            std::vector<double>& sol,
            double& objvalue) const = 0;
};

/// @}

} // mspp

#endif // SOLVER_H
