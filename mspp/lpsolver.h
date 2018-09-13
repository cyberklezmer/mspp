#ifndef LPSOLVER_H
#define LPSOLVER_H

#include "mspp/msproblem.h"

namespace mspp
{

/// \addtogroup solvers Solvers
/// @{

class sparselinearconstraint: public constraint
{
public:
    struct iitem
    {
        unsigned int index;
        double value;
    };

    sparselinearconstraint( constraint::type at = constraint::eq, double arhs=0.0 )
        : constraint(at), flhsdim(0), rhs(arhs) {}

    double rhs;
    unsigned int numnz() const { return fnzs.size(); }
    const iitem& nz(unsigned int i) const { return fnzs[i]; }

    std::shared_ptr<vector<double>> lhs() const
    {
        std::shared_ptr<vector<double>> r(new vector<double>(flhsdim,0));
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
    virtual unsigned int xdim() const  { return flhsdim; }
private:
    vector<iitem> fnzs;
    unsigned int flhsdim;
};

template <typename V>
class lpproblem : public object
{
public:
    vector<range<V>> vars;
    vector<double> f;
    vector<sparselinearconstraint> constraints;
    vector<std::string> varnames;
    bool check() const
    {
        unsigned int d=vars.size();
        if(d==0)
            return false;
        if(f.size()!=d)
            return false;
        if(varnames.size()!=d)
            return false;
        for(unsigned int i=0; i<constraints.size(); i++)
            if(constraints[i].xdim()>d)
               return false;
        return true;
    }
};

template <typename V>
class lpsolver : public object
{
public:
    virtual void solve(
            const lpproblem<V>& lp,
            vector<double>& sol,
            double& objvalue) const = 0;
};

/// @}

} // mspp


#endif // LPSOLVER_H
