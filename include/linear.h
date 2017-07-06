#ifndef LINEAR_H
#define LINEAR_H

#include <sstream>
#include <iostream>
#include <iterator>
#include "problem.h"
#include "probability.h"
#include "solution.h"

class linearfunction: public object
{
public:
    linearfunction(const std::vector<double>& acoefs)
       : coefs(acoefs) {}

    linearfunction(unsigned int dim=0)
             : coefs(dim,0.0) {}
    std::vector<double> coefs;
    unsigned int dim() const
       { return coefs.size(); }
};


using linearfunction_ptr=std::shared_ptr<linearfunction>;

class linearconstraint
{
public:
    enum type {eq, geq, leq};
    linearconstraint(const std::vector<double>& alhs,
                       type at, double arhs) :
         lhs(alhs), t(at), rhs(arhs) {}

    linearconstraint(unsigned int lhsdim, type at = eq, double arhs=0.0 )
        : lhs(lhsdim,0), t(at), rhs(arhs) {}
    std::vector<double> lhs;

    type t;
    double rhs;
    unsigned int dim() { return lhs.size(); }
};

using linearconstraint_ptr = std::shared_ptr<linearconstraint>;

using linearconstraint_list = std::vector<linearconstraint>;
using linearconstraint_list_ptr = std::shared_ptr<linearconstraint_list>;


template<class Xi>
using linearproblem = problem<linearfunction,linearconstraint,Xi>;

template<class Xi>
using linearproblem_ptr  = std::shared_ptr<linearproblem<Xi>>;

class lpsolver
{
public:
    virtual void solve(const varinfo_list& vars,
            const linearfunction& objective,
            const constraint_list<linearconstraint_ptr>& constraints,
            std::vector<double>& sol,
            double& objvalue) const = 0;
};

template<class Xi>
class biglpsolution : public treecallback
{
public:
    biglpsolution(const linearproblem_ptr<Xi>& pp,
                   const scenariotree_ptr<Xi>& sp)
     : fp(pp), fs(sp), foffsets(pp->T()+1)
    {
    }

private:

    std::vector<unsigned int> foffsets;

    //solve

    unsigned int fdim;
    linearfunction fobj;
    varinfo_list fvars;
    std::vector<linearconstraint_ptr> fconstraints;

    //list
    unsigned int fvindex;

public:
    virtual void callback(const path& p)
    {
        unsigned int stage = p.size()-1;
        unsigned int thisstagedim = fp->stagedim(stage);

        scenario<Xi> xi = fs->s(p);

        varinfo_list_ptr vars;
        constraint_list_ptr<linearconstraint> constraints;

        fp->get_constraints(stage,xi,vars,constraints);

        linearfunction_ptr objective;
        fp->get_f(stage,xi,objective);

        foffsets[stage] = fdim;
        fobj.coefs.resize(fdim + thisstagedim);
        prob up = fs->up(p);

        for(unsigned int i=0; i<thisstagedim; i++)
        {
            fobj.coefs[fdim] = up * objective->coefs[i];
            const varinfo& v = (*vars)[i];
            fvars.push_back(v);
            std::ostringstream s;
            if(v.n == "")
                s << "x" << i;
            else
                s << v.n;
            s << "@";
            for(unsigned int j=0;;)
            {
                s << p[j];
                if(++j==p.size())
                    break;
                s << "-";
            }
            fvars.rbegin()->n = s.str();
            fdim++;
        }
        if(constraints)
            for(unsigned int j=0; j<constraints->size(); j++)
            {
                linearconstraint& s = (*constraints)[j];

                linearconstraint_ptr d(new linearconstraint(fdim));

                unsigned int src=0;
                for(unsigned int i=0; i<=stage; i++)
                {
                    unsigned int dst=foffsets[i];
                    for(unsigned int r=0;
                          r<fp->stagedim(i) && src<s.lhs.size();
                          r++)
                        d->lhs[dst++] = s.lhs[src++];
                    d->rhs = s.rhs;
                    d->t = s.t;
                }
                fconstraints.push_back(d);
            }
    }

    void solve(const lpsolver& lps, treesolution_ptr& sol, double& optimal)
    {
        foffsets.resize(fp->T()+1);
        fdim = 0;
        fobj.coefs.clear();
        fvars.clear();
        fconstraints.clear();
        assert(fs->depth() == fp->T()+1);

        fs->t()->foreachnode(this);
        std::vector<double> x(fdim);

        lps.solve(fvars,fobj,fconstraints,x,optimal);
        sol.reset(new treesolution(fp->stagedims(),fs->tp(),x));
    }

private:
    linearproblem_ptr<Xi> fp;
    scenariotree_ptr<Xi> fs;
};

#endif // LINEAR_H