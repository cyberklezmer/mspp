#ifndef SDDP_H
#define SDDP_H

#include "ScenarioModel.h"
#include "DiscreteDistribution.h"
#include "SddpSolver.h"
#include "mspp/msproblem.h"
#include "mspp/process.h"
#include "mspp/mpcproblem.h"
#include "mspp/cplex.h"



namespace mspp
{

/// \addtogroup sddp SDDP
/// \ingroup sms
/// @{



class markovchain : public ddistribution<unsigned int,unsigned int>
{
public:
    unsigned int natoms() const
    {
        return nstates();
    }
private:
    virtual probability transprob(unsigned int from, unsigned int to) const = 0;
    virtual unsigned int nstates() const = 0;
    virtual atom<unsigned int>
         atom_is(unsigned int i,const unsigned int& c) const
    {
        return {i,this->transprob(c,i)};
    }
    virtual unsigned int natoms_is(const unsigned int& ) const
    {
         return nstates();
    }
};






template<typename P, typename D>
class msppsddpmodel : public ScenarioModel
{
    static constexpr bool fdebug = true;


    using C_t = typename D::C_t;
    class msspDistribution : public DiscreteDistribution
    {
        static constexpr bool fdebug = true;

    public:
        msspDistribution(const D* d, C_t& c) :
            DiscreteDistribution(mat()),
            fd(d), fc(c) {}
        // only because I have a te,úprarymember
        virtual void GenerateAnalytic(mat &sample, int count)
        {
            if(fdebug)
            {
                sys::log() << "GenerateAnalytic called" << endl
                           << "(count=" << count << ")" << endl;
            }

            sample.set_size(count, fd->dim);
            for(unsigned int i=0; i<count; i++)
            {
                rvector<double> x = fd->draw(fc);
                sample.row(i) = rowvec(x);
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
            std::is_same<typename P::C_t,nestedmcvar>::value;
    static bool constexpr fmarkov = std::is_same<typename D::C_t,unsigned int>::value;

    double lambda;
    double alpha;

    void initandcheck()
    {
        assert(fp.T()>0);
        assert(fdistrs.size()==fp.T());
        static_assert(std::is_same<typename P::V_t,realvar>::value);
        static_assert(std::is_same<typename P::I_t,rvector<double>>::value);

        // mp cvar comes reformulated
        static_assert(
             std::is_same<typename P::C_t,expectation>::value
               || fnestedcvar
              );
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
                     const rvector<double> &xi0,
                     const vector<const markovchain*>& chains)
           : ScenarioModel(p.d.size()), fp(p), fdistrs(distrs),
                            fxi0(xi0), fchains(chains)
       {
           if(fdebug)
               sys::log() << "markov ctr called" << endl;
           if(fdebug)
               sys::log() << "stages = " << p.d.size() << endl;

           assert(fmarkov);
           static_assert(std::is_same<typename D::C_t,unsigned int>::value);
           initandcheck();
       }

       msppsddpmodel(const P& p,
                     const vector<const D*>& distrs,
                     const rvector<double>& xi0 = 0)
           : ScenarioModel(p.d.size()), fp(p), fxi0(xi0), fdistrs(distrs)
       {
cout << "stages:" << p.d.size()<<endl;

           if(fdebug)
               sys::log() << "swi ctr called" << endl;
           if(fdebug)
               sys::log() << "stages = " << p.d.size() << endl;

           assert(!fmarkov);
           static_assert(std::is_same<typename D::C_t,nocondition>::value);
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
                return new DiscreteDistribution(sample);
            }
            if constexpr(fmarkov)
                return new msspDistribution(fdistrs[stage-2],state);
            else
            {
                nocondition n;
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
               const markovchain& ch = *fchains[stage-1];
               unsigned int dstn=ch.natoms();
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
                   unsigned int srcn = fchains[stage-1]->natoms();

                   mat prob(srcn,dstn);
                   for(unsigned int i=0; i<srcn; i++)
                      for(unsigned int j=0; j<dstn; j++)
                      {
                           atom<unsigned int> a = ch(j,i);
                           prob(j,i) = a.p;
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
                   sys::log() << fchains[stage-2]->natoms() <<" returned" << endl;
               return fchains[stage-2]->natoms();
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
           if(fdebug)
               sys::log() << "GetDecisionSize called (may be my call, too)" << endl
                          << "(" << stage << ")" << endl;

            unsigned int nx = fp.d[stage-1];
            if constexpr(std::is_same<typename P::C_t,nestedmcvar>::value)
            {
                if(fdebug)
                    sys::log() << nx+1 << " returned" << endl;
                return nx+1;
            }
            else
            {
                if(fdebug)
                    sys::log() << nx << " returned" << endl;
                return nx;
            }
        }
        virtual double GetRecourseLowerBound(unsigned int stage) const
       {
           if(fdebug)
               sys::log() << "GetRecourseLowerBound called" << endl;
           return std::numeric_limits<double>::min();
       }
        virtual double GetRecourseUpperBound(unsigned int stage) const
       {
           if(fdebug)
                sys::log() << "GetRecourseUpperBound called" << endl;
           return std::numeric_limits<double>::max();
       }
        virtual void FillSubradient(unsigned int stage, const double *prev_decisions, const double *scenario, double objective, const double *duals, double &recourse, double *subgradient) const
        {
           if(fdebug)
              sys::log() << "FillSubradient called" << endl;

           assert(stage >= 1);
           assert(stage <= fp.T()+1);
           if(stage == 1)
           {
               recourse = fnestedcvar ? (1-lambda) * objective : objective;


               return;
           }


           unsigned int k= stage - 1;
           bool last_stage = k == fp.T();
           ranges<typename P::V_t> vars;
           msconstraints<typename P::G_t> constraints;
           rvector<double> zeta = rvector<double>(scenario, scenario+ fdistrs[k-1]->dim);
           if(!k)
                for(unsigned int i=0; i<fxi0.size(); i++)
                    assert(fxi0[i]==scenario[i]);
           fp.x(k,zeta,vars,constraints);

           if(fdebug)
           {
               sys::log() << "(" << stage << ",(";
               for(int i=0; i<GetDecisionSize(stage-2);i++)
                   sys::log() << prev_decisions[i] << ",";
               sys::log() << "),(";
               for(int i=0; i< fdistrs[stage-1]->dim; i++)
                   sys::log() << scenario[i] << ",";
               sys::log() << "),"<< objective << ",(";
               for(int i=0; i< constraints.size(); i++)
                   sys::log() << duals[i] << ",";
               sys::log() << "))" << endl;
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

           double g = 0;
           double* sgptr = subgradient;
           unsigned int vind = 0;

           for(unsigned int i=0; i<fp.d[k-1]; i++,sgptr++)
           {
               if(fp.includedinbarx(k-1,i,k))
               {
                  for(int k=0; k< constraints.size(); k++)
                       g -= duals[k] * constraints[k].lhs(vind);
                  vind++;
               }
               if(fnestedcvar)
               {
                   double parentu = prev_decisions[nv-1];

                   *sgptr = (1-lambda) * g -
                        objective > parentu ? lambda / alpha * g : 0;
               }
               else
                   *sgptr = g;
           }
           if(fnestedcvar) // sg for u
               *sgptr = objective > prev_decisions[nv-1]
                           ? - lambda / alpha : 0;
        }

        virtual double CalculateUpperBound(unsigned int stage, const double *prev_decisions, const double *decisions, const double *scenario, double recourse_estimate, bool cut_tail) const
        {
           if(fdebug)
           {
              sys::log() << "CalculateUpporBound called" << endl;
           }
           assert(!cut_tail); // still not ready for this

           assert(stage <= fp.T()+1);
           assert(stage > 0);
           unsigned int k= stage - 1;

           bool root_node = (stage == 1);
           bool last_stage = (stage == GetStagesCount());

           ranges<typename P::V_t> vars;
           msconstraints<typename P::G_t> constraints;
           rvector<double> zeta =
                k ? rvector<double>(scenario, scenario+ fdistrs[k-1]->dim)
                  : fxi0;
           if(!k)
                for(unsigned int i=0; i<fxi0.size(); i++)
                    assert(fxi0[i]==scenario[i]);

           typename P::F_t f = fp.f(k,zeta);

           vector<double> x(decisions,decisions+fp.d[k]);
           double node_value = f(x);

           //coefficients from current stage
           if(!fnestedcvar)
               return node_value + recourse_estimate;
           else
           {
               double u = last_stage ? 0 : decisions[fp.d[k]];

               if(root_node)
               {
                   return recourse_estimate + (1-lambda) * node_value
                            + lambda * u;
               }
               else
               {
                   double prevu = decisions[fp.d[k-1]];
                   double rv = node_value + lambda * u;
                   double tv = (1-lambda) * rv;
                   if(tv > prevu)
                       tv += lambda / alpha * (tv > prevu);
                   return tv;
               }
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
           unsigned int k= stage - 1;

           if(fdebug)
           {
               sys::log() << "(" << stage << ",(";
               if(stage > 1)
                  for(int i=0; i<GetDecisionSize(stage-2);i++)
                     sys::log() << prev_decisions[i] << ",";
               sys::log() << "),(";
               for(int i=0; i< !k ? fxi0.size() : fdistrs[stage-1]->dim; i++)
                   sys::log() << scenario[i] << ",";
               sys::log() << "))" << endl;
           }

           bool first_stage = !k;
           bool last_stage = k == fp.T();
           ranges<typename P::V_t> vars;
           msconstraints<typename P::G_t> constraints;
           rvector<double> zeta=
                   k ? rvector<double>(scenario, scenario+ fdistrs[k-1]->dim)
                     : fxi0;
           if(!k)
                for(unsigned int i=0; i<fxi0.size(); i++)
                    assert(fxi0[i]==scenario[i]);

           fp.x(k,zeta,vars,constraints);

           typename P::F_t f = fp.f(k,zeta);

           IloNumVarArray x(env,vars.size());

           assert(vars.size() == fp.barxdim(k));
           for(int i=0; i<vars.size(); i++)
           {
               const range<realvar>& v = vars[i];

               const IloNum l = v.islinf() ? -IloInfinity : v.l();
               const IloNum h = v.ishinf() ? IloInfinity : v.h();
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
                   for(unsigned int i=0; i<fp.d[k-1];i++)
                   {
                        if(fp.includedinbarx(k-1,i,k))
                            rhs -= c.lhs(srccoef++) * *srcv;
                        srcv++;
                   }
                }
                unsigned int srcvar = 0;
                for(unsigned int i=0; i<fp.barxdim(k); i++)
                {
                   if(fp.includedinbarx(k,i,k))
                       v+= c.lhs(srccoef++)*x[i];
                   srcvar++;
               }
               IloRange r =
                   c.t() == constraint::eq ? IloRange(env,c.rhs(),v,c.rhs())
                   : (c.t() == constraint::leq ? IloRange(env,v,c.rhs())
                       : IloRange(env,c.rhs(),v));
               v.end();
               model.add(r);
               dual_constraints.push_back(r);
           }
           if(fdebug)
           {
               std::ostringstream s;
               s << "model" << (rand() % 100) << ".lp";

               IloCplex cplex(model);
               sys::log() << "Partial model exported to " << s.str() << endl;
               cplex.exportModel(s.str().c_str());
           }
        }
private:
//        static constexpr double EPSILON = 0.000001;

        const P& fp;
        vector<const D*> fdistrs;
        rvector<double> fxi0;
        vector<const markovchain*> fchains;
};


template <typename P>
class sddpsolution : public variables<typename P::V_t>
{
public:
    sddpsolution(const P& p) :
        variables<typename P::V_t>(p.d[0]) {}
    void set(const variables<typename P::V_t>& v)
    {
        assert(v.size()==this->size());
        for(unsigned int i=0; i<v.size(); i++)
            (*this)[i] = v[i];
    }
};


template <typename P, typename D /* typename O*/>
class sddpmethod : public object
{
public:
    sddpmethod()
    {
        static_assert(
             std::is_same<typename P::C_t,expectation>::value
               || std::is_same<typename P::C_t,mpmcvar>::value
               || std::is_same<typename P::C_t,nestedmcvar>::value
              );
    }

    static void solve(
             const P& p,
             processdistribution<D>& xi,
             double& optimal,
             sddpsolution<P>& sol)
    {
        assert(xi.xi0.size()==0);
        vector<const D*> ds;
        for(unsigned int i=0; i<xi.size(); i++)
            ds.push_back(&xi(i));
        return solve(p,ds,xi.xi0(),optimal,sol);
    }

    static void solve(
             const P& p,
             const vector<const D*>& xi,
             const rvector<double> xi0,
             double& optimal,
             sddpsolution<P>& sol)
    {
        if constexpr(std::is_same<typename P::C_t,mpmcvar>::value)
        {
            mpmcvarequivalent<P> e(p);
            sddpsolution<mpmcvarequivalent<P>> es;
            sddpmethod<mpmcvarequivalent<P>,D> sm;
            sm.solve(e, xi, xi0, optimal, es);
            for(unsigned int i=0; i<sol.size(); i++)
                sol[i]=es[i];
        }
        else
        {
            msppsddpmodel<P,D> mm(p,xi,xi0);
            do_solve(&mm,optimal,sol);
        }

    }

    static void solve(
             const P& p,
             const vector<const D*>& xi,
             const rvector<double> xi0,
             const vector<const markovchain*>& m,
             double& optimal,
             sddpsolution<P>& sol)
    {
        if constexpr(std::is_same<typename P::C_t,mpmcvar>::value)
        {
            mpmcvarequivalent<P> e(p);
            sddpsolution<mpmcvarequivalent<P>> es;
            sddpmethod<mpmcvarequivalent<P>,D> sm;
            sm.solve(e, xi, xi0, m, optimal, es);
            for(unsigned int i=0; i<sol.size(); i++)
                sol[i]=es[i];
        }
        else
        {
            msppsddpmodel<P,D> mm(p,xi,xi0,m);
            do_solve(&mm,optimal,sol);
        }
    }
private:
    static void do_solve(
             msppsddpmodel<P,D> *m,
             double& optimal,
             sddpsolution<P>& sol)
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
        config.samples_per_stage = DESCENDANTS;
        config.reduced_samples_per_stage = REDUCED_DESCENDANTS;
        config.solver_strategy = STRATEGY_DEFAULT; // STRATEGY_CONDITIONAL
        // config.cut_nodes_not_tail = true;
        //there are more settings, for instace:
        //config.debug_solver = true;
        SddpSolver solver(m, config);

        solver.Solve(weights_o, lb_o, ub_o_m, ub_o_b);
        time_o = time(NULL) - time_o;

        //printout
        sys::log() << "Original problem: " << weights_o << endl;
        sys::log() << "Original problem @ " << time_o << "s, lb: " << lb_o << ", ub_mean:" << ub_o_m << ", ub_bound:" << ub_o_b << endl;

        variables<typename P::V_t> vars;
        unsigned int sn;
        if constexpr(std::is_same<typename P::C_t,nestedmcvar>::value)
            sn = weights_o.n_cols;
        else
            sn = weights_o.n_cols - 1;
        for(unsigned int i=0; i<sn; i++)
            vars.push_back(weights_o(i));
        sol.set(vars);
    }
};




/// @}


template <typename I=double>
struct dsitem
{
    unsigned int state;
    I x;
};


template <typename I=double>
class laststate : public zeta<unsigned int>
{
public:
    laststate(const scenario<dsitem<I>>& s)
    {
        assert(s.size());
        fstate = s[s.size()-1].state;
    }
    operator unsigned int () const { return fstate; }
private:
    unsigned int fstate;
};

/*template <typename I=double>
class lastxi : public zeta<I>
{
public:
    lastxi(const scenario<dsitem<I>>& s)
    {
        assert(s.size());
        fx = s[s.size()-1].x;
    }
    operator I () const { return fx; }
private:
    I fx;
};*/



template <typename D, typename M>
class statedrivendistribution :
        public distribution<dsitem<typename D::I_t>,laststate<typename D::I_t>>
{
    void checktypes() const
    {
        static_assert(
                    std::is_same<typename D::C_t,nocondition>::value
                    || std::is_same<typename D::C_t,unsigned int>::value);
    }
public:

    statedrivendistribution(const ptr<D>& d, const ptr<M>& m):
      fd(d),fm(m)
    {
        checktypes();
    }
    statedrivendistribution(const D& d, const M& m):
      fd(new D(d)),fm(new M(m))
    {
        checktypes();
    }

    D& d() const { return fd; }
    M& m() const { return fm; }
private:
    ptr<D> fd;
    ptr<M> fm;
};


} // namespace

#endif // SDDP_H
