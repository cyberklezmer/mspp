#ifndef BIGLP_H
#define BIGLP_H

#include "mspp/linear.h"
#include "mspp/lpsolver.h"

namespace mspp
{

template<class Xi>
class biglpmethod :
        public solutionmethod<linearproblem<Xi>,scenariotree<Xi>,
                                      treesolution>,
        public treecallback
{
public:
    biglpmethod(const linearproblem_ptr<Xi>& pp,
                   const scenariotree_ptr<Xi>& sp,
                 const lpsolver_ptr& lps )
     : solutionmethod<linearproblem<Xi>,scenariotree<Xi>,treesolution>(pp,sp), flps(lps),
       foffsets(pp->T()+1)
    {
    }

private:


    unsigned int fdim;
    linearfunction fobj;
    varinfo_list fvars;
    std::vector<sparselinearconstraint_ptr> fconstraints;
    std::vector<std::string> fvarnames;

    // stored expectation constraints
    std::vector<linearconstraint> fsrcexpconstraints;
    std::vector<sparselinearconstraint*> fdstexpconstraints;
    std::vector<prob> fups;


    // list
    unsigned int fvindex;

public:
    virtual void callback(const path& p)
    {
        unsigned int stage = p.size()-1;
        prob up = this->fd->up(p);
        scenario<Xi> xi = this->fd->s(p);

        unsigned int thisstagedim = this->fp->stagedim(stage);

        // calling original problems \p constraints
        varinfo_list_ptr vars;
        constraint_list_ptr<linearconstraint> constraints;

        this->fp->constraints(stage,xi,vars,constraints);

        linearfunction_ptr objective;
        this->fp->f(stage,xi,objective);

        foffsets[stage] = fdim;
        fobj.coefs.resize(fdim + thisstagedim);

        // this stage variables and objective

        for(unsigned int i=0; i<thisstagedim; i++)
        {
            fobj.coefs[fdim] = up * objective->coefs[i];
            const varinfo& v = (*vars)[i];
            fvars.push_back(v);
            std::ostringstream s;
            s << this->fp->varname(stage,i) << "@";
            for(unsigned int j=0;;)
            {
                s << p[j];
                if(++j==p.size())
                    break;
                s << "-";
            }
            fvarnames.push_back(s.str());
            fdim++;
        }

        // stored constraints with expectations
        for(unsigned int i=0; i<fsrcexpconstraints.size(); i++)
        {
            assert(i<fdstexpconstraints.size());

            for(unsigned int src=this->fp->stageoffset(stage), dst=foffsets[stage];
                src<this->fp->dimupto(stage); src++,dst++)
            {
                if(src < fsrcexpconstraints[i].lhs.size())
                {
                     double coef = fsrcexpconstraints[i].lhs[src];
                     if(coef)
                        fdstexpconstraints[i]->set_lhs(dst,coef*up/fups[i]);
                }
            }

        }

        // this stage constraints
        if(constraints)
            for(unsigned int j=0; j<constraints->size(); j++)
            {
                linearconstraint& s = (*constraints)[j];

                sparselinearconstraint_ptr d(new sparselinearconstraint());
                unsigned int src=0;
                unsigned int i=0;
                for(; i<=stage; i++)
                {
                    unsigned int dst=foffsets[i];
                    for(unsigned int r=0;
                          r<this->fp->stagedim(i) && src<s.lhs.size();
                          r++)
                    {
                        double v = s.lhs[src++];
                        if(v)
                           d->set_lhs(dst,v);
                        dst++;
                    }
                    d->rhs = s.rhs;
                    d->t = s.t;
                }
                if(s.lhs.size() > this->fp->dimupto(stage))
                {
                    assert(stage < this->fp->T());
                    fsrcexpconstraints.push_back(s);
                    fdstexpconstraints.push_back(d.get());
                    fups.push_back(up);
                }
                fconstraints.push_back(d);
            }
    }

    void solve( treesolution_ptr& sol, double& optimal)
    {
        foffsets.resize(this->fp->T()+1);
        fdim = 0;
        fobj.coefs.clear();
        fvars.clear();
        fconstraints.clear();
        fvarnames.clear();
        assert(this->fd->depth() == this->fp->T()+1);

        fsrcexpconstraints.clear();
        fdstexpconstraints.clear();
        fups.clear();


        this->fd->t()->foreachnode(this);
        std::vector<double> x(fdim);

        flps->solve(fvars,fobj,fconstraints,fvarnames,x,optimal);
        sol.reset(new treesolution(this->fp->stagedims(),this->fd->tp(),x));
    }

private:
    lpsolver_ptr flps;
    std::vector<unsigned int> foffsets;

};

}

#endif // BIGLP_H
