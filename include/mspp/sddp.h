#ifndef SDDP_H
#define SDDP_H

#include "mspp/solution.h"

/// \addtogroup sddp SDDP
/// \ingroup sms
/// @{


namespace mspp
{

using sddpsolution = std::vector<variable>;


template <typename Xi>
class markovindexer : public object
{
public:
    markovindexer(const scenario<Xi>& s): fs(s)
    {}
    unsigned int i() const { return index_is(fs); }
    operator unsigned int() const { return i(); }
private:
    virtual unsigned index_is(const scenario<X>& s) const = 0;
    scenario<X> fs;
};

template <typename Xi>
class trivialindexer : public markovindexer<Xi>
{
public:
    trivialindexer(const scenario<Xi>& s)
        : markovindexer<Xi>(s) {}
private:
    unsigned int index_is() const { return 0; }
};

template<typename F, typename C>
using sddpproblem=problem<realvar, F, interstagelinearconstraint, C>;

template <typename Xi>
using sddppdistribution=processdistribution<Xi,mcdistribution<Xi,nocondition<Xi>>>;

template <typename Xi, typename F, typename M>
class markovsddpmethod : public solutionmethod
         <sddpproblem<F,lastvalue<Xi>>,
          sddpdistribution<Xi>,
          sddpsolution, lpsolver, bool>
{
public:
    markovsddpmethod(const sddpproblem<F,C>& p,
               const sddpdistribution<Xi>& z,
               const markovcondition& m,
               const lpsolver& lps )
     : solutionmethod
        <sddpproblem<F,C>,sddpdistribution<Xi>,sddpsolution,lpsolver,bool>
         (p,z,lps)
    {}
private:
    markovcondition& fm;
};

template <typename Xi, typename F>
using sddpmethod=markovsddpmethod<Xi,F,trivialindexer<Xi>>;

/// @}


} // namespace

#endif // SDDP_H
