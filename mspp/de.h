#ifndef DE_H
#define DE_H

#include "mspp/msproblem.h"
#include "mspp/distribution.h"
#include "mspp/lpsolver.h"
#include "mspp/stsolution.h"

namespace mspp
{

/// \addtogroup de Deterministic equivalent
/// \ingroup sms
/// @{



template <typename P, typename S>
using desolution=stsolution<P,S>;


template<typename P,typename S>
class demethod : public object, public sctreecallback<typename S::X_t>
{
public:
    bool solve(
             const P& p,
             const S& z,
             const lpsolver& lps,
             double& optimal,
             desolution<P,S>& sol)
        {
            fp = &p;
            fz = &z;
            foffsets.resize(this->fp->T()+1);
            fdim = 0;
            fobj.clear();
            fvars.clear();
            fconstraints.clear();
            fvarnames.clear();
            assert(this->fz->T() == this->fp->T());

            this->fz->foreachnode(this);
            std::vector<double> x(fdim);
            lps.solve(fvars,fobj,fconstraints,fvarnames,x,optimal);
            sol.set(x);
            return true;
        }

private:
    const P* fp;
    const S* fz;
    unsigned int fdim;
    std::vector<double> fobj;
    vardefs<realvar> fvars;
    std::vector<sparselinearconstraint> fconstraints;
    std::vector<std::string> fvarnames;

public:
    virtual void callback(const indexedhistory<typename S::X_t>& a)
    {
        unsigned int stage = a.size()-1;
        probability up = a.uncprob();
        scenario<typename S::X_t> s = a.s();

        unsigned int thisstagedim = this->fp->d[stage];

        // calling original problems \p stageinfo

        vardefs<realvar> vars;
        msconstraints<typename P::G_t> constraints;

        demethod::fp->stageinfo(stage,s,vars,constraints);

        linearfunction f = this->fp->f(stage,s);

        foffsets[stage] = fdim;
        fobj.resize(fdim + thisstagedim);

        // this stage variables and f

        for(unsigned int i=0; i<thisstagedim; i++)
        {
            fobj[fdim] = up * f.c(i);
            const realvar& v = vars[i];
            fvars.push_back(v);
            std::ostringstream s;
            s << this->fp->varname(stage,i) << "@";
            for(unsigned int j=0;;)
            {
                s << a[j].i;
                if(++j==a.size())
                    break;
                s << "-";
            }
            fvarnames.push_back(s.str());
            fdim++;
        }

        // this stage msconstraints
        if(constraints.size())
            for(unsigned int j=0; j<constraints.size(); j++)
            {
                const typename P::G_t& s = constraints[j];

                sparselinearconstraint d;
                unsigned int src=0;
                unsigned int i=0;
                for(; i<=stage; i++)
                {
                    unsigned int dst=foffsets[i];
                    for(unsigned int r=0;
                          r<this->fp->d[i] && src<s.lhssize();
                          r++)
                    {
                        double v = s.lhs(src++);
                        if(v)
                           d.set_lhs(dst,v);
                        dst++;
                    }
                    d.rhs = s.rhs();
                    d.t = s.t();
                }
                fconstraints.push_back(d);
            }
    }

private:
    std::vector<unsigned int> foffsets;
};

/// @}


} // namespace

#endif // DE_H
