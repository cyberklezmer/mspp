#ifndef SOLUTION_H
#define SOLUTION_H

#include "msproblem.h"
#include "process.h"

namespace mspp
{

template <typename P, typename D, typename Z, typename X, typename O>
class solution : public object
{
public:
    using P_t = P;
    using D_t = D;
    using Z_t = Z;

protected:
    solution(const pair<ptr<X>,O>& s) :
        fx(s.first), fo(s.second)
    {
        static_assert(
             std::is_base_of<processdistribution<typename D::D_t,typename D::Z_t>,D>::value,
                 "D has to be a descendant of processdistribuition");
        static_assert(
             std::is_base_of<zeta<typename Z::X_t,typename Z::R_t>,Z>::value,
                 "Z has to be a descendant of zeta");

        static_assert(
             std::is_convertible<typename Z::R_t,typename P::Y_t>::value,
                        "Z and P are incompatible");
        static_assert(
             std::is_convertible<typename Z::X_t,typename D::X_t>::value,
                        "Z and D are incompatible");

    }
public:
    const O& obj() const { return fo; }
    ptr<X> x() const { return fx; }
private:
    ptr<X> fx;
    O fo;
};

} // namespace

#endif // SOLUTION_H
