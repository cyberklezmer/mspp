#ifndef DE_H
#define DE_H

#include <utility>
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

template<typename P,typename S>
struct demethodstate
{
    const P& p;
    const S& x;

    unsigned int dim;
    ptr<lpproblem<typename P::V_t>> lp;
    vector<unsigned int> offsets;
};

class demethod : public object
{
    template <typename P, typename S, typename L>
    class decallback: public sctreecallback<typename S::I_t, demethodstate<P,S>>
    {
    public:
        using V_t = typename P::V_t;
        void solve(
                 const P& p,
                 const S& xi,
                 stsolution<P,S>& sol)
        {
            static_assert(
                 std::is_same<typename P::O_t,expectation>::value
                   || std::is_same<typename P::O_t,mpmcvar>::value
                   || std::is_same<typename P::O_t,nestedmcvar>::value
                  );
            static_assert(
                 std::is_same<typename P::I_t, typename S::I_t>::value
                        );

            if constexpr(std::is_same<typename P::O_t,expectation>::value)
            {
                assert(p.T() == xi.T());

                demethodstate<P,S> s={p,xi};
                s.lp.reset(new lpproblem<V_t>);
                s.offsets.resize(p.T()+1);
                s.dim = 0;

                xi.foreachnode(this, &s);

                vector<double> x(s.dim);
                L solver;
                double optimal;
                solver.solve(*(s.lp),x,optimal);
                sol.set(x,optimal);
            }
            else if constexpr(std::is_same<typename P::O_t,mpmcvar>::value)
            {
                mpmcvarequivalent<P> e(p);
                stsolution<mpmcvarequivalent<P>,S> es(e,xi);
                demethod::solve<mpmcvarequivalent<P>,S,L>(e, xi, es);

                stsolreducer<stsolution<mpmcvarequivalent<P>,S>,
                             stsolution<P,S>> r;
                r.convert(es,sol);
            }
        }
        virtual void callback(const indexedpath<typename S::I_t>& a,
                              demethodstate<P,S>* state) const
        {
            unsigned int stage = a.size()-1;

            probability up = a.uncprob();
            typename S::I_t xi = a.pth();

            unsigned int thisstagedim = state->p.d[stage];

            // calling original problems \p x

            vector<range<V_t>> vars;
            vector<typename P::G_t> constraints;

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

            unsigned int dst=state->offsets[stage];
            for(unsigned int src=0; src<state->p.d[stage];)
                lp.f[dst++] += up * f.c(src++);

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
public:
    template <typename P,typename S, typename L>
    static void solve(
             const P& p,
             const S& xi,
             stsolution<P,S>& sol)
    {
        static_assert
         (std::is_convertible<S&,scenariotree<typename S::I_t>&>::value,
          "S has to be a descendant of scenariotree");
        decallback<P,S,L> ds;
        ds.solve(p,xi,sol);
    }
};

/// @}



} // namespace

#endif // DE_H
