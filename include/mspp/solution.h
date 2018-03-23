#ifndef SOLUTION_H
#define SOLUTION_H

#include "mspp/problem.h"
#include "mspp/distributions.h"

namespace mspp
{
/// \addtogroup sms Solution Methods
/// @{


template<typename P, typename Z, typename X, typename S, typename R>
class solutionmethod : public object
{
public:
    solutionmethod(const P& problem,
                   const Z& process,
                   const S& solver) :
       fp(problem), fz(process), fs(solver)
    {
    }
    virtual R solve(double& optvalue, X& optsolution) = 0;
protected:
    const P& fp;
    const Z& fz;
    const S& fs;
};


/// \addtogroup solutions Optimal solutions
/// \ingroup Solution
/// @{


template <typename X>
class sctreesolution : public object, public sctreecallback<std::vector<variable>>
{
public:
    sctreesolution(const problemstructure& ps,
               const scenariotree<X>& st) :
               fps(ps), fst(st)
    {}
public:
    void gothrough()
    {
        fst->foreachnode(this);
    }
    virtual void assign(const std::vector<variable>& x)
    {
        fx = x;
    }
    std::vector<variable> &x() {return fx; }
    void callback(const indexedhistory<X>& s)
    {
    }

private:
    problemstructure fps;
    scenariotree<X> fst;
    std::vector<variable> fx;
};

/*
template<typename P, typename X, typename S, typename R>
class sctreesolutionmethod
   : public solutionmethod<P,scenariotree<X>,sctreesolution<X>,S, R>
{
public:
    sctreesolutionmethod(const ptr<P>& problem,
                   const ptr<scenariotree<X>>& process,
                   const ptr<S>& solver) :
       solutionmethod<P,scenariotree<X>,
           sctreesolution<X>,S, R> (problem,process,solver)
    {}
};
*/

/// @}

/// @}



} // namespace

#endif // SOLUTION_H
