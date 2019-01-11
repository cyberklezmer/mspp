#ifndef DE_H
#define DE_H

#include <utility>
#include "mspp/solution.h"
#include "mspp/process.h"
#include "mspp/lpsolver.h"
#include "mspp/mpcproblem.h"

namespace mspp
{

/// \addtogroup de Deterministic equivalent
/// \ingroup sms
/// @{

template<typename V, typename R>
struct dexitem
{
    R zeta;
    variables<V> x;
};

/// \brief Solution of \ref msproblem with \ref treedistribution
/// as stochastic parameter
/// \tparam V type of the problems's variables
/// \tparam D the distribution
/// \tparam Z a mapping from scenarios into the problem's random parameters
///
/// \p D has to be a descendant of \ref treedistribution,
/// \p Z has to be a descendant of \ref zeta
template <typename V, typename D, typename Z>
class dex : public
     treedistribution<dexitem<V,typename Z::R_t>,
        typename D::E_t, unsigned int>
{
public:
    using E_t = typename D::E_t;
    using X_t = dexitem<V,typename Z::R_t>;

    class exportcb : public tdcallback<X_t,variables<V>>
    {
        virtual void callback(const scenario<X_t>& xi,
                              probability up,
                              vector<V>* state) const
        {
            assert(xi.size());
            unsigned int k=xi.size()-1;

            variables<V> src = xi[k].x;
            for(unsigned int i=0; i<src.size(); i++)
                state->push_back(src[i]);
        }
    };

    dex(const D& d,const msproblemstructure& aps,
        const vectors<bool> excl = 0):
        treedistribution
          <dexitem<V,typename Z::R_t>,typename D::E_t,unsigned int>(d.T()),
        vdistribution<dexitem<V,typename Z::R_t>,nothing>(d.T()+1),
        fxi(ptr<D>(new D(d))), fps(aps), fexcl(excl)
       {
           static_assert(std::is_base_of<
              treedistribution<typename D::X_t,
                               typename D::E_t,
                               typename D::A_t>,D>::value);
        }
/*    desolution(const P& p, const ptr<D> d) :
        fxi(d), fps(p.d)
       {}*/
    void set(const variables<V>& x)
    {
        fx=ptr<variables<V>>(new variables<V>(x));
    }
    void set(const ptr<variables<V>> x)
    {
        fx=x;
    }

    void exportlinear(variables<V>& v) const
    {
        if(fexcl.size()==0)
            v = *fx;
        else
        {
            exportcb cb;
            this->foreachnode(&cb,&v);
        }
    }
    void set(const vectors<bool> excl)
    {
        assert(fps.size()>=excl.size());
        fexcl=excl;
    }
    const msproblemstructure& ps() const { return fps; }
private:

    virtual void branches_are(const vector<E_t>& e,
         vector<E_t>& es) const
    {
        fxi->branches(e,es);
    }

    virtual void beforefe(unsigned int& s) const
    {
        s=0;
    }
    virtual void afterfe(unsigned int& a) const
    {
        assert(a==fx->size());
    }


    virtual void index2sinfo(const vector<E_t>& e,
                                X_t& x,
                                probability& up,
                                unsigned int& s) const
    {
        assert(e.size());
        typename D::A_t a;
        scenario<typename D::X_t> sc;

        up = 1;
        for(unsigned int i=0; i<e.size(); i++)
        {
            up *= e[i].p;
            sc.push_back(e[i].x);
        }

        unsigned int k=e.size()-1;
        vector<V> v;
        for(unsigned int i=0;i<fps[k]; i++)
        {
           assert(s < fx->size());

           if(k>=fexcl.size() || i >= fexcl[k].size() || !fexcl[k][i])
              v.push_back((*fx)[s]);
           s++;
        }

        x.zeta = Z()(sc);
        x.x = v;
    }

    ptr<variables<V>> fx;
    vectors<bool> fexcl;
    ptr<D> fxi;
    msproblemstructure fps;
};

template <typename P, typename D>
struct destate
{
    const P& p;
    const D& x;

    unsigned int dim;
    ptr<lpproblem<typename P::V_t>> lp;
    vector<unsigned int> offsets;
};

template <typename P, typename D, typename Z>
class decallback: public tdcallback<typename D::X_t,destate<P,D>>
{
public:
    using X_t = typename D::X_t;
    using V_t = typename P::V_t;

    virtual void callback(const scenario<X_t>& xi,
                          probability up,
                          destate<P,D>* state) const
    {
        unsigned int stage = xi.size()-1;

        unsigned int thisstagedim = state->p.d()[stage];

        // calling original problems \p x

        vector<range<V_t>> vars;
        vector<typename P::G_t> constraints;

        typename Z::R_t z = Z()(xi);
        state->p.x(stage,z,vars,constraints);

        linearfunction f = state->p.f(stage,z);

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
            s<< i;
//                for(unsigned int j=0;;)
//                {
//                    s << i;
//                    if(++j==a.size())
//                        break;
//                    s << "-";
//                }
            lp.varnames.push_back(s.str());
            state->dim++;
        }

        // this stage objective

        unsigned int dst=state->offsets[stage];
        for(unsigned int src=0; src<state->p.d()[stage];)
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
                      r<state->p.d()[i];
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

/// \brief Optimal solution computed by DE method
/// \tparam P the problem solution of which the class is
/// \tparam D distribution
/// \tparam Z mapping from scenarios to random parameters of \p P
/// \tparam X optimal decision variables
/// \tparam O optimal values
/// \tparam L linear solver
template <typename P, typename D, typename Z, typename L>
class desolution : public solution<P,D,Z,dex<typename P::V_t,D,Z>,double>
{
    using P_t = P;
    using D_t = D;
    using Z_t = Z;
    using X_t = typename D::X_t;
    using V_t = typename P::V_t;

public:
    desolution(const P& p, const D& xi) :
        solution<P,D,Z,dex<V_t,D,Z>,double>(solve(p,xi))
    {
    }

    static pair<ptr<dex<V_t,D,Z>>,double> solve(
             const P& p,
             const D& xi)
    {
/*        static_assert
         (std::is_convertible<S&,scenariotree<typename S::I_t>&>::value,
          "S has to be a descendant of scenariotree");*/

        static_assert(
             std::is_same<typename P::O_t,expectation>::value
               || std::is_same<typename P::O_t,mpmcvar>::value
               || std::is_same<typename P::O_t,nestedmcvar>::value
              );
        static_assert(
             std::is_same<typename P::Y_t, typename Z::R_t>::value
                    );

        if constexpr(std::is_same<typename P::O_t,expectation>::value)
        {
            assert(p.T() == xi.T());

            destate<P,D> s={p,xi};
            s.lp.reset(new lpproblem<V_t>);
            s.offsets.resize(p.T()+1);
            s.dim = 0;

            decallback<P,D,Z> ds;

            using dst = destate<P,D>;

            xi.foreachnode(&ds, &s);

            vector<double> x(s.dim);

            double optimal;

            L solver;
            solver.solve(*(s.lp),x,optimal);
            ptr<dex<V_t,D,Z>> sol(new dex<V_t,D,Z>(xi,p.d()));
            sol->set(x);
            return {sol,optimal};
        }
        else if constexpr(std::is_same<typename P::O_t,mpmcvar>::value)
        {
            mpmcvarequivalent<P> e(p);
            desolution<mpmcvarequivalent<P>,D,Z,L> es(e,xi);
            ptr<dex<typename P::V_t,D,Z>> sol = es.x();
            vectors<bool> excl(p.T()+1);
            for(unsigned int i=0; i<=p.T(); i++)
                excl[i].push_back(true);
            for(unsigned int i=1; i<p.T(); i++)
                excl[i].push_back(true);
            sol->set(excl);
            return { sol, es.obj() };
        }
    }
};

/// @}



} // namespace


#endif // DE_H
