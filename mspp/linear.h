#ifndef LINEAR_H
#define LINEAR_H

#include "mspp/commons.h"

namespace mspp
{

class linearconstraint : virtual public constraint
{
public:
    linearconstraint(unsigned int xdim) :
        constraint(xdim), flhs(xdim), frhs(0) {}
    linearconstraint(const vector<double>& lhs, double rhs=0,
            constraint::type t=constraint::eq) :
        constraint(lhs.size(),t), flhs(lhs), frhs(rhs) {}
    double lhs(unsigned int i) const
    {
        assert(i<flhs.size());
        return flhs[i];
    }
    void setlhs(unsigned int i, double v)
    {
        assert(i<flhs.size());
        flhs[i]=v;
    }
    void setlhs(const vector<double>& lhs)
    {
        assert(lhs.size()==xdim());
        flhs = lhs;
    }
    void setrhs(double rhs ) { frhs = rhs; }
    double rhs() const { return frhs; }
private:
    vector<double> flhs;
    double frhs;
};

} // namespace

#endif // LINEAR_H
