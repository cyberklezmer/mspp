#ifndef DE_H
#define DE_H

#include "mspp/msproblem.h"
#include "mspp/random.h"
#include "mspp/lpsolver.h"
#include "mspp/mpcproblem.h"
#include "mspp/stsolution.h"

namespace mspp
{

/// \addtogroup de Deterministic equivalent
/// \ingroup sms
/// @{

template<typename P,typename Z>
struct demethodstate
{
    const P& p;
    const Z& z;

    unsigned int dim;
    ptr<lpproblem<typename P::V_t>> lp;
    vector<unsigned int> offsets;
};

template<typename P,typename Z, typename O>
class demethod : public object, public sctreecallback<typename Z::X_t, demethodstate<P,Z>>
{
public:
    using V_t = typename P::V_t;
    demethod()
    {
        static_assert(
             std::is_same<typename P::C_t,expectation>::value
               || std::is_same<typename P::C_t,mpmcvar>::value
               || std::is_same<typename P::C_t,nestedmcvar>::value
                    );
    }

    bool solve(
             P& p,
             Z& z,
             double& optimal,
             stsolution<P,Z>& sol)
        {
            if constexpr(std::is_same<typename P::C_t,expectation>::value)
            {
                assert(p.T() == z.T());

                demethodstate<P,Z> s={p,z};
                s.lp.reset(new lpproblem<V_t>);
                s.offsets.resize(p.T()+1);
                s.dim = 0;

                z.foreachnode(this, &s);

                vector<double> x(s.dim);
                O solver;
//                lpsolver<typename P::V_t> lps ;
                solver.solve(*(s.lp),x,optimal);
                sol.set(x);
                return true;
            }
            else if constexpr(std::is_same<typename P::C_t,mpmcvar>::value)
            {
                mpmcvarequivalent<P> e(p);
                stsolution<mpmcvarequivalent<P>,Z> es(e,z);
                demethod<mpmcvarequivalent<P>,Z,O> m;
                m.solve(e, z, optimal, es);

                stsolreducer<stsolution<mpmcvarequivalent<P>,Z>,
                             stsolution<P,Z>> r;
                r.convert(es,sol);
                return true;
            }
            return false;
        }
public:
    virtual void callback(const indexedpath<rvector<typename Z::X_t>>& a,
                          demethodstate<P,Z>* state) const
    {
        unsigned int stage = a.size()-1;

        probability up = a.uncprob();
        scenario<typename Z::X_t> s = a.pth();

        unsigned int thisstagedim = state->p.d[stage];

        // calling original problems \p x

        ranges<V_t> vars;
        msconstraints<typename P::G_t> constraints;

        state->p.x(s,vars,constraints);

        linearfunction f = state->p.f(s);

        state->offsets[stage] = state->dim;

        lpproblem<V_t>& lp = *(state->lp);

        lp.f.resize(state->dim + thisstagedim);

        // this stage variables and initialization of the objective

        for(unsigned int i=0; i<thisstagedim; i++)
        {
            lp.f[state->dim] = 0;
            const range<V_t>& v = vars[i];
            lp.vars.push_back(v);
            std::ostringstream s;
            s << state->p.varname(stage,i) << "@";
            for(unsigned int j=0;;)
            {
                s << a[j].i;
                if(++j==a.size())
                    break;
                s << "-";
            }
            lp.varnames.push_back(s.str());
            state->dim++;
        }

        // this stage objective

        unsigned int src=0;

        for(unsigned int i=0; i<=stage; i++)
        {
            unsigned int dst=state->offsets[i];
            for(unsigned int r=0; r<state->p.d[i]; r++)
            {
                if(state->p.includedinbarx(i,r,stage))
// std::cout << stage << ":" << up << ":" <<  f.c(src) << std::endl;
                    lp.f[dst] += up * f.c(src++);
                dst++;
            }
        }

        // this stage msconstraints
        for(unsigned int j=0; j<constraints.size(); j++)
        {
            const linearmsconstraint& s = constraints[j];

            sparselinearconstraint d;
            unsigned int src=0;

            for(unsigned int i=0; i<=stage; i++)
            {
                unsigned int dst=state->offsets[i];
                for(unsigned int r=0;
                      r<state->p.d[i];
                      r++)
                {
                    if(state->p.includedinbarx(i,r,stage))
                    {
                        double v = s.lhs(src++);
                        if(v)
                           d.set_lhs(dst,v);
                    }
                    dst++;
                }
            }
            d.rhs = s.rhs();
            d.settype(s.t());
            lp.constraints.push_back(d);
        }
    }
};

/// @}



} // namespace

#endif // DE_H
