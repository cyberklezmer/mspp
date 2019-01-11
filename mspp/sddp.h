#ifndef SDDP_H
#define SDDP_H

#include "ScenarioModel.h"
#include "DiscreteDistribution.h"
#include "SddpSolver.h"
#include "mspp/msproblem.h"
#include "mspp/process.h"
#include "mspp/markov.h"
#include "mspp/mpcproblem.h"
#include "mspp/solution.h"
#include "mspp/cplex.h"

#include <mcheck.h>

// zdalipak zobrazuju zetou?

namespace mspp
{

/// \addtogroup sddp SDDP
/// \ingroup sms
/// @{

template <typename X>
class laststate : public zeta<pair<unsigned int, X>,unsigned int>
{
public:
    virtual unsigned int operator() (const scenario<pair<unsigned int, X>>& s) const
    {
        assert(s.size());
        return s[s.size()-1].first;
    }
};

class onlylaststate : public mapping<pair<unsigned int, unsigned int>, unsigned int>
{
public:
    virtual unsigned int operator() (const pair<unsigned int, unsigned int>& d) const
    {
        return d.second;
    }
};

/// \brief Hidden Markov chain distribution
/// \tparam M the Markov chain
/// \tparam Y distribution conditioned by <tt>unsigned int</tt>
template<typename M, typename Y>
using hmcdistribution = mijdistribution<M,Y,onlylaststate>;

/// Finite support hidden Markov chain distribution
template<typename M, typename D >
using fhmcdistribution = fijdistribution<M,D,onlylaststate>;

/// \brief Hidden Markov Chain process
/// \tparam M Markov Chain
/// \tparam Y distribution, conditioned by <tt>unsigned int</tt>
///
template<typename M, typename Y>
class hmcprocessdistribution: public processdistribution<
        hmcdistribution<M,Y>,
        laststate<typename Y::I_t>>
{
public:
    using M_t = M;
    using Y_t = Y;
    hmcprocessdistribution(
                  typename Y::I_t xi0,
                  const M& md, const Y& xid, unsigned int T) :
        processdistribution<hmcdistribution<M,Y>,
         laststate<typename Y::I_t>>({0,xi0},hmcdistribution<M,Y>(md,xid))
    {}

    static vector<hmcdistribution<M,Y>>
       makexi(const vector<M>& m, const vector<Y>& xi )
    {
        vector<hmcdistribution<M,Y>> ret;
        assert(xi.size()==m.size());
        for(unsigned int i=0; i<xi.size(); i++)
            ret.push_back(hmcdistribution<M,Y>(m[i],xi[i]));
        return ret;
    }
    hmcprocessdistribution(
           typename Y::I_t xi0,
           const vector<M>& m,
           const vector<Y>& xi
            ) :
        processdistribution<hmcdistribution<M,Y>,
           laststate<typename Y::I_t>>({0,xi0},makexi(m,xi))
    {}

    struct init
    {
        typename Y::I_t xi0;
        const vector<M>& m;
        const vector<Y>& xi;
    };

    hmcprocessdistribution(const init& i) :
        processdistribution<hmcdistribution<M,Y>,
           laststate<typename Y::I_t>>
             (diracdistribution<pair<unsigned int, typename Y::I_t>>({0,i.xi0}),
              makexi(i.m,i.xi))
    {}
//
};


template <typename X>
class onlyx : public mapping<pair<unsigned int, X>, X>
{
public:
    virtual X operator() (const pair<unsigned int, X>& d) const
    {
        return d.second;
    }
};


/// Scenario restriction
template <typename X>
using lastmdxi = lastxi<pair<unsigned int, X>,onlyx<X>>;


template <typename V>
using sddpx = variables<V>;

//         assert(fsv.size()==fv.size());
// for(unsigned int i=0; i<fv.size(); i++)
//    fv[i] = fsv[i];

/// SDDP objective value
struct sddpobj
{
    double lb() const { return flb; }
    double ubm() const { return fubm; }
    double ubb() const { return fubb; }
    time_t time() const { return ftime;}
    void set(double lb, double ubm, double ubb, time_t t)
    {
        flb=lb;
        fubm=ubm;
        fubb=ubb;
        ftime = t;
    }
private:
    double flb;
    double fubm;
    double fubb;
    time_t ftime;
};


template<typename P, typename D, typename M, bool markov>
class msppsddpmodel : public ScenarioModel
{
    static constexpr bool fdebug = true;
    static constexpr bool fmarkov = markov;
    static unsigned int countDebugModel()
    {
        static unsigned int cnt=0;
        return ++cnt;
    }

    using C_t = typename D::C_t;
    class msspDistribution : public DiscreteDistribution
    {
        static constexpr bool fdebug = true;

    public:
        msspDistribution(const D* d, C_t c) :
            DiscreteDistribution(mat()),
            fd(d), fc(c) {}
        // only because I have a temporarymember
        virtual void GenerateAnalytic(mat &sample, int count)
        {
            if(fdebug)
            {
                sys::log() << "GenerateAnalytic called" << endl
                           << "(count=" << count << ")" << endl;
            }

            sample.set_size(count, fd->dim());
            for(unsigned int i=0; i<count; i++)
            {
                if constexpr(fmarkov)
                {
                    vector<double> x = M()({fc,fd->draw(fc)});
                    sample.row(i) = rowvec(x);
                }
                else
                {
                    vector<double> x = M()(fd->draw(fc));
                    sample.row(i) = rowvec(x);
                }
            }
            if(fdebug)
                sys::log() << sample <<  " returned" << endl;
        }

        virtual double DistributionFunction(double x) const
        {
            assert(0);
        }
        virtual double InverseDistributionFunction(double y) const
        {
            assert(0);
        }

    private:
        const D* fd;
        C_t fc;
    };

    static constexpr bool fnestedcvar =
            std::is_same<typename P::O_t,nestedmcvar>::value
           || std::is_same<typename P::O_t,mpmcvar>::value;

    double lambda;
    double alpha;

    void initandcheck()
    {
        assert(fxi0.size());
        assert(fp.T()>0);
        assert(fdistrs.size()==fp.T());
        static_assert(std::is_same<typename P::V_t,realvar>::value);
        static_assert(std::is_same<typename P::Y_t,vector<double>>::value);

        for(unsigned int i=0; i<=fp.T(); i++)
        {
            assert(fp.maxf(i) <  max<double>());
            assert(fp.minf(i) >  min<double>());
        }
        // mp cvar comes reformulated
        static_assert(
             std::is_same<typename P::O_t,expectation>::value
               || fnestedcvar
              );
        if constexpr(fnestedcvar && std::is_same<typename P::O_t,mpmcvar>::value)
           assert(fp.T() < 2);
        if constexpr(fnestedcvar)
        {
            alpha = fp.rho.alpha;
            assert(alpha>0);
            lambda = fp.rho.lambda;
        }
        else
            alpha = lambda = 0;
        // check if the nested problem depend only on previous decision
        // variables
        typename P::R_t r;
        for(unsigned int k=1; k<=fp.T(); k++)
            for(unsigned int i=0; i<k-1; i++)
                 assert(!r.included(i,k));
    }

public:
       msppsddpmodel(const P& p,
                     const vector<const D*>& distrs,
                     const vector<double> &xi0,
                     const vector<const mcdistribution*>& chains)
           : ScenarioModel(p.d().size()), fp(p), fdistrs(distrs),
                            fxi0(xi0), fchains(chains)
       {
           if(fdebug)
               sys::log() << "markov ctr called" << endl;
           if(fdebug)
               sys::log() << "stages = " << p.d().size() << endl;

           assert(fmarkov);
           static_assert(std::is_same<typename D::C_t,unsigned int>::value);
           initandcheck();
       }

       msppsddpmodel(const P& p,
                     const vector<const D*>& distrs,
                     const vector<double>& xi0 = 0)
           : ScenarioModel(p.d().size()), fp(p), fxi0(xi0), fdistrs(distrs)
       {
           if(fdebug)
               sys::log() << "swi ctr called" << endl;
           if(fdebug)
               sys::log() << "stages = " << p.d().size() << endl;

           assert(!fmarkov);
           static_assert(std::is_same<typename D::C_t,nothing>::value);
           initandcheck();
       }

       virtual ~msppsddpmodel() {}

       // would be much easier if I provided only the sample
       virtual Distribution* GetDistribution(unsigned int stage, unsigned int state) const
       {
           if(fdebug)
               sys::log() << "GetDistribution called" << endl
                          << "(" << stage << "," << state << ")" << endl;

            assert(stage > 0);
            assert(stage-1 <= fdistrs.size());
            if(stage==1)
            {
                mat sample(1,fxi0.size());
                sample.row(0)= rowvec(fxi0);
                DiscreteDistribution* d = new DiscreteDistribution(sample);
                d->UnbiasedEstimate(); // to initialize mu_, which is later used in
                             // determining the dimension of rv.'s
                return d;
            }
            if constexpr(fmarkov)
                return new msspDistribution(fdistrs[stage-2],state-1);
            else
            {
                nothing n;
                return new msspDistribution(fdistrs[stage-2],n);
            }
       }

        virtual mat GetTransitionProbabilities(unsigned int stage) const
        {
           if(fdebug)
               sys::log() << "GetTransitionProbabilities called" << endl;

           assert(stage > 0);
           assert(stage-1 <= fchains.size());
           if(fmarkov)
           {
               const mcdistribution& ch = *fchains[stage-1];
               unsigned int dstn=ch.nstates();
               if(stage == 1)
               {
                   mat prob(1,dstn);
                   for(unsigned int i=0; i<dstn; i++)
                   {
                       atom<unsigned int> a = ch(i,0);
                       prob(0,i) = a.p;
                   }
                   return prob;
               }
               else
               {
                   unsigned int srcn = fchains[stage-2]->nstates();

                   mat prob(srcn,dstn);
                   for(unsigned int i=0; i<srcn; i++)
                      for(unsigned int j=0; j<dstn; j++)
                      {
                           atom<unsigned int> a = ch(j,i);
                           prob(i,j) = a.p;
                       }
                   return prob;
               }
           }
           else
               return ones(1);
        }
        virtual unsigned int GetStatesCountStage(unsigned int stage) const
        {
           if(fdebug)
               sys::log() << "GetStatesCountStage called" << endl
                          << "(" << stage << ")" << endl;

           assert(stage > 0);
           if(fmarkov && stage>1)
           {
               assert(stage-1 <= fchains.size());
               if(fdebug)
                   sys::log() << fchains[stage-2]->nstates() <<" returned" << endl;
               return fchains[stage-2]->nstates();
           }
           else
           {
               if(fdebug)
                   sys::log() << "1 returned" << endl;
               return 1;
           }
        }

        virtual StageDependence GetStageDependence() const
        {
           if(fdebug)
               sys::log() << "GetStageDependence called" << endl;

           return fmarkov ? MARKOV : STAGE_INDEPENDENT;
        }

        virtual unsigned int GetDecisionSize(unsigned int stage) const
        {
           if(0 && fdebug)
               sys::log() << "GetDecisionSize called (may be my call, too)" << endl
                          << "(" << stage << ")" << endl;

            unsigned int nx = fp.d()[stage-1];
            if constexpr(fnestedcvar)
            {
                unsigned int rv;
                if(stage-1<fp.T())
                    rv = nx+1;
                else
                    rv = nx;
                if(0 && fdebug)
                    sys::log() << rv << " returned" << endl;
                return rv;
            }
            else
            {
                if(0 && fdebug)
                    sys::log() << nx << " returned" << endl;
                return nx;
            }
        }
        virtual double GetRecourseLowerBound(unsigned int stage) const
       {
           if(fdebug)
               sys::log() << "GetRecourseLowerBound called" << endl;
           return fp.minf(stage-1);
       }
        virtual double GetRecourseUpperBound(unsigned int stage) const
       {
           if(fdebug)
                sys::log() << "GetRecourseUpperBound called" << endl;
           return fp.maxf(stage-1);
       }
        virtual void FillSubradient(unsigned int stage,
                                    const double *prev_decisions,
                                    const double *scenario,
                                    double objective,
                                    const double *duals,
                                    double &recourse,
                                    double *subgradient) const
        {
           if(fdebug)
              sys::log() << "FillSubradient called";

           assert(stage >= 1);
           assert(stage <= fp.T()+1);
           if(stage == 1)
           {
               recourse = fnestedcvar ? (1-lambda) * objective : objective;
               if(fdebug)
                  sys::log() << "(1)" << endl << "recourse=" << recourse << " returned" << endl;

               return;
           }
           unsigned int k= stage - 1;
           vector<range<typename P::V_t>> vars;
           vector<typename P::G_t> constraints;
           vector<double> zeta(scenario, scenario+ fdistrs[k-1]->dim());
           if(!k)
                for(unsigned int i=0; i<fxi0.size(); i++)
                    assert(fxi0[i]==scenario[i]);
           fp.x(k,zeta,vars,constraints);
           if(fdebug)
           {
               sys::log() << "(stage=" << stage << ",pd=(";
               for(int i=0; i<GetDecisionSize(stage-1);i++)
                   sys::log() << prev_decisions[i] << ",";
               sys::log() << "),sc=(";
               for(int i=0; i< fdistrs[stage-2]->dim(); i++)
                   sys::log() << scenario[i] << ",";
               sys::log() << "), o="<< objective << ",duals=(";
               for(int i=0; i< constraints.size(); i++)
                   sys::log() << duals[i] << ",";
               sys::log() << "))" << endl;
           }

           if(fdebug)
           {
               sys::log() << "fp.x called, " << vars.size() << " vars, "
                           << constraints.size() << " constraints returned"
                           << endl;
               for(unsigned int i=0; i<constraints.size(); i++)
               {
                   for(unsigned int j=0; j<constraints[i].xdim(); j++)
                        sys::log() << constraints[i].lhs(j) << ", ";
                   sys::log() << endl;
               }
           }

           unsigned int nv = GetDecisionSize(stage-1);
           if(fnestedcvar)
           {
               double parentu = prev_decisions[nv-1];
               recourse = (1-lambda) * objective;
               if (objective > parentu)
                   recourse += (lambda / alpha) * (objective - parentu);
           }
           else
               recourse = objective;

           double* sgptr = subgradient;
           unsigned int vind = 0;

           for(unsigned int i=0; i<fp.d()[k-1]; i++,sgptr++)
           {
               double g = 0;
               if(fp.includedinbarx(k-1,i,k))
               {
                  for(int k=0; k< constraints.size(); k++)
                       g -= duals[k] * constraints[k].lhs(vind);
                  vind++;
               }
               if(fnestedcvar)
               {
                   double parentu = prev_decisions[nv-1];

                   *sgptr = (1-lambda) * g +
                        (objective > parentu ? lambda / alpha * g : 0);
               }
               else
                   *sgptr = g;
           }
           if(fnestedcvar) // sg for u
               *sgptr++ = objective > prev_decisions[nv-1]
                           ? - lambda / alpha : 0;
           assert(sgptr-subgradient == GetDecisionSize(stage-1));
           if(fdebug)
           {
               sys::log() << "g=(";
               for(unsigned int i=0; i<nv; i++)
                    sys::log() << subgradient[i] << ",";
               sys::log() << "), recourse=" << recourse << " returned" << endl << endl;
               sys::log().flush();
           }
       }

        virtual double CalculateUpperBound(unsigned int stage,
                                           const double *prev_decisions,
                                           const double *decisions,
                                           const double *scenario,
                                           double recourse_estimate,
                                           bool cut_tail) const
        {
           if(fdebug)
           {
              sys::log() << "CalculateUpporBound called(stage=" << stage << ")" << endl;
           }
           assert(!cut_tail); // still not ready for this

           assert(stage <= fp.T()+1);
           assert(stage > 0);
           unsigned int k = stage - 1;

           bool root_node = (stage == 1);
           bool last_stage = (stage == GetStagesCount());

           vector<typename P::V_t> vars;
//          vector<typename P::G_t> constraints;
           vector<double> zeta =
                k ? vector<double>(scenario, scenario+ fdistrs[k-1]->dim())
                  : fxi0;
           if(!k)
                for(unsigned int i=0; i<fxi0.size(); i++)
                    assert(fxi0[i]==scenario[i]);

           if(fdebug)
           {
               sys::log() << "(" << stage << ",pd=(";
               if(k)
                   for(int i=0; i<GetDecisionSize(stage-1);i++)
                       sys::log() << prev_decisions[i] << ",";
               sys::log() << "),decisions=(";
               for(int i=0; i<GetDecisionSize(stage);i++)
                   sys::log() << decisions[i] << ",";
               sys::log() << "),sc=(";
               for(int i=0; i< zeta.size(); i++)
                   sys::log() << scenario[i] << ",";
               sys::log() << "), re="<< recourse_estimate << ")" << endl;
           }

           typename P::F_t f = fp.f(k,zeta);

           vector<double> x(decisions,decisions+fp.d()[k]);
           double node_value = f(x);

           if(fdebug)
               sys::log() << "node value = " << node_value << endl;

           //coefficients from current stage
           if(!fnestedcvar)
           {
               if(fdebug)
                   sys::log() << "returned: "
                      << node_value+ recourse_estimate << endl;

               return node_value + recourse_estimate;
           }
           else
           {
               double u = last_stage ? 0 : decisions[fp.d()[k]];

               double res;
               if(root_node)
                   res =  recourse_estimate + node_value + u;
               else
               {
                   double prevu = prev_decisions[fp.d()[k-1]];

                   double rv;
                   if(last_stage)
                       rv = node_value;
                   else
                       rv = recourse_estimate + node_value + lambda * u;
                   res = (1-lambda) * rv;
                  if(rv > prevu)
                       res += lambda / alpha * (rv - prevu);
               }
               if(fdebug)
                   sys::log() << res << " returned." << endl;
               return res;
           }
        }
        virtual void BuildCoinModel(CoinModelWrapper * coin_model, unsigned int stage, const double *prev_decisions, const double *scenario, vector<string> &decision_vars, vector<string> &dual_constraints) const {}
        virtual void BuildCplexModel(
               IloEnv &env,
               IloModel &model,
               IloExpr &objective,
               unsigned int stage,
               const double *prev_decisions,
               const double *scenario,
               vector<IloNumVar> &decision_vars,
               vector<IloRange> &dual_constraints) const
        {
           if(fdebug)
               sys::log() << "BuildCplexModel called" << endl;

           assert(stage <= fp.T()+1);
           assert(stage > 0);
           unsigned int k = stage - 1;

           if(fdebug)
           {
               sys::log() << "(" << stage << ",(";
               if(stage > 1)
                  for(int i=0; i<GetDecisionSize(stage-1);i++)
                     sys::log() << prev_decisions[i] << ",";
               sys::log() << "),(";
               unsigned int d = !k ? fxi0.size() : fdistrs[stage-2]->dim();
               for(int i=0; i<d ; i++)
                   sys::log() << scenario[i] << ",";
               sys::log() << "))" << endl;
           }

           bool first_stage = !k;
           bool last_stage = k == fp.T();
           vector<range<typename P::V_t>> vars;
           vector<typename P::G_t> constraints;
           vector<double> zeta=
                   k ? vector<double>(scenario, scenario+ fdistrs[k-1]->dim())
                     : fxi0;
           if(!k)
                for(unsigned int i=0; i<fxi0.size(); i++)
                    assert(fxi0[i]==scenario[i]);

           fp.x(k,zeta,vars,constraints);

           typename P::F_t f = fp.f(k,zeta);

           IloNumVarArray x(env,vars.size());
           assert(vars.size() == fp.xdim(k));
           for(int i=0; i<vars.size(); i++)
           {
               const range<realvar>& v = vars[i];

               IloNum l = v.islinf() ? -IloInfinity : v.l();
               IloNum h = v.ishinf() ? IloInfinity : v.h();

               IloNumVar::Type t;
               switch(v.t())
               {
                   case range<realvar>::realt:
                      t = IloNumVar::Float;
                   break;
                   case range<realvar>::intt:
                      t = IloNumVar::Int;
                   break;
                   case range<realvar>::bint:
                      t = IloNumVar::Bool;
                   break;
               default: assert(0);
               }

               std::string vn = fp.varname(stage-1,i);

               x[i]=IloNumVar(env,l,h,t,vn.c_str());
               double c = f.c(i);
               if(fdebug)
                   sys::log() << "Adding "<< c << "*" << vn << " to objective"  << endl;
               objective += c * x[i];
               decision_vars.push_back(x[i]);
           }
           model.add(x);

           //objective according to the risk aversion lambda
           if (!last_stage && fnestedcvar)
           {
              //last stage has no recourse
              //var var = variable to calculate CVaR = VaR level
              IloNumVar u(env, GetRecourseLowerBound(stage), GetRecourseUpperBound(stage),"u");
              objective += lambda * u;
              model.add(u);
              decision_vars.push_back(u);
           }

           for(int h=0; h< constraints.size(); h++)
           {
               IloExpr v(env);

               const linearconstraint& c = constraints[h];
               double rhs=c.rhs();
               unsigned int srccoef=0;
               if(!first_stage)
               {
                   // thx to the check from ctr, we are sure that no < stage-1 vars are involved in the constraint

                   const double* srcv = prev_decisions;
                   for(unsigned int i=0; i<fp.d()[k-1];i++)
                   {
                        if(fp.includedinbarx(k-1,i,k))
                            rhs -= c.lhs(srccoef++) * *srcv;
                        srcv++;
                   }
                }
                unsigned int srcvar = 0;
                for(unsigned int i=0; i<fp.xdim(k); i++)
                {
                   if(fp.includedinbarx(k,i,k))
                       v+= c.lhs(srccoef++)*x[i];
                   srcvar++;
               }
               IloRange r =
                   c.t() == constraint::eq ? IloRange(env,rhs,v,rhs)
                   : (c.t() == constraint::leq ? IloRange(env,v,rhs)
                       : IloRange(env,rhs,v));
               v.end();
               model.add(r);
               dual_constraints.push_back(r);
           }
           if(0 && fdebug)
           {
               std::ostringstream s;
               s << "model" << countDebugModel() << ".lp";

               IloCplex cplex(model);

               cplex.setOut(env.getNullStream());
               sys::log() << "Partial model exported to " << s.str() << endl;
               cplex.exportModel(s.str().c_str());
           }
        }

   pair<ptr<sddpx<typename P::V_t>>,sddpobj> solve()
   {
       const unsigned int SEED = 350916502;

       const unsigned int DESCENDANTS = 100; //100
       const unsigned int REDUCED_DESCENDANTS = 0; //0 == no reduction
       const unsigned int DERIVATIVE_ITERATIONS = 10;
       //vectors for collecting stats and printout


       sys::log() << "Starting solving:" << endl;

       //set the seed
       if (SEED > 0) {
               sys::seed(SEED + 1); //fixes the tree for each iteration
       }

       //vector for statistics and output
       rowvec weights_o;
       double lb_o;
       double ub_o_m;
       double ub_o_b;
       time_t time_o;

       time_o = time(NULL);

       //configure the SddpSolver
       SddpSolverConfig config;
       config.debug_solver = 1;
       config.samples_per_stage = DESCENDANTS;
       config.reduced_samples_per_stage = REDUCED_DESCENDANTS;
       config.solver_strategy = STRATEGY_DEFAULT; // STRATEGY_CONDITIONAL
       // config.cut_nodes_not_tail = true;
       //there are more settings, for instace:
       //config.debug_solver = true;
       SddpSolver solver(this, config);

       solver.Solve(weights_o, lb_o, ub_o_m, ub_o_b);
       time_o = time(NULL) - time_o;

       //printout
       sys::log() << "Original problem: " << weights_o << endl;
       sys::log() << "Original problem @ " << time_o << "s, lb: " << lb_o << ", ub_mean:" << ub_o_m << ", ub_bound:" << ub_o_b << endl;

       ptr<sddpx<typename P::V_t>> vars(new sddpx<typename P::V_t>);
       unsigned int sn;
       if constexpr(std::is_same<typename P::O_t,nestedmcvar>::value ||
                    std::is_same<typename P::O_t,mpmcvar>::value)
           sn = weights_o.n_cols - 1;
       else
           sn = weights_o.n_cols;
       for(unsigned int i=0; i<sn; i++)
           vars->push_back(weights_o(i));
       sddpobj obj;
       obj.set(lb_o, ub_o_m, ub_o_b, time_o);
       return { vars, obj};
   }

private:
//        static constexpr double EPSILON = 0.000001;

        const P& fp;
        vector<const D*> fdistrs;
        vector<double> fxi0;
        vector<const mcdistribution*> fchains;
public:
    static pair<ptr<sddpx<typename P::V_t>>,sddpobj> solve
        (const P& p,const vector<const D*>& ds, const vector<double>& xi0)
    {
       if constexpr(std::is_same<typename P::O_t,mpmcvar>::value)
       {
            if(p.T()>=2)
            {
                mpmcvarequivalent<P> e(p);
                msppsddpmodel<mpmcvarequivalent<P>,D,M,false> mm(e,ds,xi0);
                return mm.solve();
            }
        }
       msppsddpmodel<P,D,M,false> mm(p,ds,xi0);
       return mm.solve();

     }
    static pair<ptr<sddpx<typename P::V_t>>,sddpobj> solve
       (const P& p,const vector<const D*>& ds,
            const vector<const mcdistribution*>& ms, const vector<double>& xi0 )
    {
        if constexpr(std::is_same<typename P::O_t,mpmcvar>::value)
        {
            if(p.T()>=2)
            {
                mpmcvarequivalent<P> e(p);
                msppsddpmodel<mpmcvarequivalent<P>,D,M,true> mm(e,ds,xi0,ms);
                return mm.solve();
            }
        }
        msppsddpmodel<P,D,M,true> mm(p,ds,xi0,ms);
        return mm.solve();
    }
};

/// \brief Base class of SDDP solutions
/// \tparam P the problem solved
/// \tparam D the distribution (descendant of \ref processdistribution)
/// \tparam Z the restriction (has to be \p llastxi)
/// \tparam O the solver (has to be \p cplex<realvar>)
///
template <typename P, typename D, typename Z, typename O>
class sddpsolbase : public solution<P,D,Z,sddpx<typename P::V_t>,sddpobj>
{
protected:
    static constexpr void check()
    {
        static_assert(std::is_same<O,cplex<realvar>>::value);
        static_assert(std::is_same
              <Z,lastxi<vector<double>, typename Z::M_t>>::value);
        static_assert(
             std::is_same<typename P::O_t,expectation>::value
               || std::is_same<typename P::O_t,mpmcvar>::value
               || std::is_same<typename P::O_t,nestedmcvar>::value,
                 "only problems with expectation, mp- or nested- cvar can be solved"
              );
        static_assert(
             std::is_base_of<processdistribution
                <typename D::D_t, typename D::Z_t>,D>::value);
        /*        static_assert(
           (markov && std::is_same<typename D::D_t::S_t::I_t,vector<double>>::value)
             || std::is_same<typename D::X_t,vector<double>>::value,
                 "only distributions with X=vector<double> may be used"   );*/
    }

public:
    using P_t = P;
    using D_t = D;
    using Z_t = Z;
    using X_t = typename D::X_t;
    using V_t = typename P::V_t;

    sddpsolbase(const pair<ptr<sddpx<typename P::V_t>>,sddpobj>& p) :
        solution<P,D,Z,sddpx<typename P::V_t>,sddpobj>(p)
    {}
};

/// \brief Solution by plain SDDP
/// \tparam D the process distribution (has to be stage-wise independent)
///
/// For other template parameters, see \ref sddpsolbase
template <typename P, typename D, typename Z, typename O>
class sddpsolution: public sddpsolbase<P,D,Z,O>
{
public:
    sddpsolution(const P& p, const D& d)
        : sddpsolbase<P,D,Z,O>(solve(p,d))
    {

        static_assert(std::is_same
           <typename D::Z_t, noxi<typename D::X_t>>::value,
                      "D has to ve stagewise independent");
    }
private:
    static pair<ptr<sddpx<typename P::V_t>>,sddpobj> solve
        (const P& p,const D& xi)
    {
        vector<const typename D::D_t*> ds;
        for(unsigned int i=1; i<=xi.T(); i++)
            ds.push_back(&(xi.d(i)));
        return
         msppsddpmodel<P, typename D::D_t, typename Z::M_t, false>
                ::solve(p,ds,xi.e().x());
    }
};

/// \brief Solution by Markov SDDP
///
/// The distribution of componentso of \p D has to be
/// a joint distribution (\ref jdistribution) of \ref fmcdistribution and a descendant
/// of \ref mdistribution conditioned by <tt>unsigned int</tt>.
/// For the other template parameters, see \ref sddpsolbase
///
template <typename P, typename D, typename Z, typename O>
class msddpsolution: public sddpsolbase<P,D,Z,O>
{
public:
    using D_t = typename D::D_t;
    using F_t = typename D_t::F_t;
    using S_t = typename D_t::S_t;

    msddpsolution(const P& p, const D& d)
        : sddpsolbase<P,D,Z,O>(solve(p,d))
    {
        static_assert(std::is_base_of<
            ijdistribution
               <typename D_t::F_t,typename D_t::S_t,typename D_t::M_t>, D_t>::value);

        static_assert(std::is_same<typename D_t::M_t, onlylaststate>::value);

        static_assert(std::is_base_of<mcdistribution,F_t>::value);

        static_assert(std::is_base_of<
             mdistribution<typename S_t::I_t, unsigned int>,
             S_t>::value);

        static_assert(std::is_same
           <typename D::Z_t, laststate<vector<double>>>::value);
    }
private:
    static pair<ptr<sddpx<typename P::V_t>>,sddpobj> solve
        (const P& p,const D& xi)
    {
        vector<const typename D_t::S_t*> ds;
        vector<const mcdistribution*> ms;
        for(unsigned int i=1; i<=xi.T(); i++)
        {
            ds.push_back(&(xi.d(i).second()));
            ms.push_back(&(xi.d(i).first()));
        }
        return
         msppsddpmodel<P, typename D_t::S_t, typename Z::M_t, true>
                ::solve(p,ds,ms,xi.e().x().second);
    }
};


/*
template <typename P, typename D, typename Z, typename O>
class msddpsolution : public sddpsolbase<P,D,Z,O,true>
{
public:
    using D_t = typename D::D_t;
    using F_t = typename D_t::F_t;
    using S_t = typename D_t::S_t;
    static void constexpr check()
    {
    }

    msddpsolution(const P& p, const D& d) : sddpsolbase<P,D,Z,O,true>(p,d)
    {
        check();

    }
    solve()
    {
    }
};
*/
/*
template <typename P, typename D, typename Z, typename O>
class msddpsolution : public solution<P,D,Z,sddpx<typename P::V_t>,sddpobj>
{
public:
    msddpmethod()
    {
        static_assert(
             std::is_same<typename P::C_t,expectation>::value
               || std::is_same<typename P::C_t,mpmcvar>::value
               || std::is_same<typename P::C_t,nestedmcvar>::value
              );
        static_assert(
             std::is_same<typename D::X_t::I_t,vector<double>>::value);
    }

    static void solve(
             const P& p,
             const D& xi,
             sddpsolution<P>& sol)
    {
        for(unsigned int i=1; i<=xi.T(); i++)
        {
            ds.push_back(&(xi.d(i).second()));
            ms.push_back(&(xi.d(i).first()));
        }

    }


    static void solve(
             const P& p,
             const vector<const typename D::X_t*>& xi,
             const typename D::X_t::I_t& xi0,
             const vector<const mcdistribution*>& m,
             sddpsolution<P>& sol)
    {
        if constexpr(std::is_same<typename P::C_t,mpmcvar>::value)
        {
            if(p.T()>=2)
            {
                mpmcvarequivalent<P> e(p);
                sddpsolution<mpmcvarequivalent<P>> es(e);
                msddpmethod<mpmcvarequivalent<P>,D> sm;
                sm.solve(e, xi, xi0, m, es);

                variables<typename P::V_t> v(sol.firststage().size());
                for(unsigned int i=0; i<v.size(); i++)
                    v[i]=es.firststage()[i];
                sol.set(v,es.lb(),es.ubm(),es.ubb(),es.time());
                return;
            }
        }
        msppsddpmodel<P,typename D::X_t> mm(p,xi,xi0,m);
        mm.solve(sol);
    }
};

*/

/// @}





} // namespace

#endif // SDDP_H
