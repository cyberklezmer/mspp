#ifndef CVAR_H
#define CVAR_H

#include "linear.h"


template<class Xi>
class mcvarproblem: public linearproblem<Xi>
{
protected:
    static std::vector<unsigned int>
       dims(const std::vector<unsigned int>& d)
    {
        assert(d.size()>0);
        std::vector<unsigned int> r(d.size());
        r[0]=d[0]+1;
        int i=1;
        for(; i< d.size()-1; i++)
            r[i]=d[i]+2;
        r[i] = d[i]+1;
        return r;
    }

    mcvarproblem(const linearproblem_ptr<Xi>& lp, double alpha, double lambda) :
        linearproblem<Xi>(dims(lp->stagedims())),
        flp(lp), falpha(alpha), flambda(lambda)
    {
        assert(flp->T()>0);
    }
public:
    virtual std::string varname(unsigned int stage, unsigned int i) const
    {
        if(i<flp->stagedim(stage))
            return flp->varname(stage,i);
        if(stage == 0)
            return "u";
        if(stage == this->T())
            return "theta";
        else if(i==flp->stagedim(stage))
            return "theta";
        else
            return "u";
    }

protected:
    double mu() const { return 1.0 - flambda; }
    double nu() const { return 1.0 - flambda + flambda / falpha; }

    void origvarsandconstraints(
            unsigned int stage,
            const scenario<Xi>& xih,
            varinfo_list_ptr& vars,
            linearconstraint_list_ptr& constraints
            ) const
    {
        bool laststage = stage == this->T();

        unsigned int newsize = this->dimupto(stage);

        flp->get_constraints(stage,xih,vars,constraints);

        if(stage >0)
            vars->push_back(varinfo(varinfo::R));
        if(!laststage)
            vars->push_back(varinfo(varinfo::R));

        if(constraints)
        {
            for(linearconstraint_list::iterator p = constraints->begin();
                 p != constraints->end(); p++)
            {
                assert(p->lhs.size()==flp->dimupto(stage));
                p->lhs.resize(newsize);

                unsigned int src = flp->dimupto(stage)-1;
                unsigned int dst = newsize-1;

                for(unsigned int i=stage; i>0; i--)
                {
                    p->lhs[dst--] = 0;
                    if(i<this->T())
                        p->lhs[dst--] = 0;
                    for(int j=0; j<flp->stagedim(i); j++)
                        p->lhs[dst--] = p->lhs[src--];
                }
                p->lhs[dst] = 0; // the rest of the constraint is in the right place
            }
        }
    }

    linearproblem_ptr<Xi> flp;
    double falpha;
    double flambda;
};

template<class Xi>
class mmpcvarproblem: public mcvarproblem<Xi>
{
public:
    mmpcvarproblem(const linearproblem_ptr<Xi>& lp, double alpha, double lambda) :
        mcvarproblem<Xi>(lp, alpha, lambda)
    {
    }

    virtual void f(
            unsigned int stage,
            const scenario<Xi>& xih,
            linearfunction_ptr& f) const
    {
        unsigned int k = this->stagedim(stage);
        f.reset(new linearfunction(k));
        if(!stage)
        {
            linearfunction_ptr orig;
            this->flp->get_f(0,xih,orig);
            for(unsigned int i=0; i<this->flp->stagedim(0); i++)
                f->coefs[i]=orig->coefs[i];
        }
        f->coefs[k-1]=1;
        if(stage > 0 && stage < this->T())
            f->coefs[k-2]=1;
    }

    virtual void constraints(
            unsigned int stage,
            const scenario<Xi>& xih,
            varinfo_list_ptr& vars,
            linearconstraint_list_ptr& constraints
            ) const
    {
        bool laststage = stage == this->T();

        unsigned int newsize = this->dimupto(stage);

        this->origvarsandconstraints(stage,xih, vars,constraints);

        if(stage)
        {
            if(!constraints)
                 constraints.reset(new linearconstraint_list);

            linearfunction_ptr c;
            this->flp->get_f(stage,xih,c);

            constraints->push_back(linearconstraint(newsize,linearconstraint::geq));
            constraints->push_back(linearconstraint(newsize,linearconstraint::geq));
            linearconstraint_list::reverse_iterator muc = constraints->rbegin();
            linearconstraint_list::reverse_iterator nuc = muc + 1;

            unsigned int i = newsize-(laststage ? 1 : 2);
            muc->lhs[i] = 1;
            nuc->lhs[i--] = 1;

            for(int j=this->flp->stagedim(stage)-1; j>=0; j--)
            {
                muc->lhs[i] = -this->mu() * c->coefs[j];
                nuc->lhs[i--] = -this->nu() * c->coefs[j];
            }

            muc->lhs[i] = this->mu();
            nuc->lhs[i] = this->nu();
        }
    }
};

template<class Xi>
using mmpcvarproblem_ptr = std::shared_ptr<mmpcvarproblem<Xi>>;


template<class Xi>
class mncvarproblem: public mcvarproblem<Xi>
{

public:
    mncvarproblem(const linearproblem_ptr<Xi>& lp, double alpha, double lambda) :
         mcvarproblem<Xi>(lp, alpha, lambda)
    {
    }


    virtual void f(
            unsigned int stage,
            const scenario<Xi>& xih,
            linearfunction_ptr& f) const
    {
        unsigned int k = this->stagedim(stage);
        f.reset(new linearfunction(k));
        if(stage == 0)
        {
            linearfunction_ptr orig;
            this->flp->get_f(0,xih,orig);
            unsigned int i=0;
            for(; i<this->flp->stagedim(0); i++)
                f->coefs[i]=orig->coefs[i];
            f->coefs[i] = 1;
        }
        else if(stage==1)
            f->coefs[k-2] = 1;
    }

    virtual void constraints(
            unsigned int stage,
            const scenario<Xi>& xih,
            varinfo_list_ptr& vars,
            linearconstraint_list_ptr& constraints
            ) const
    {
        bool laststage = stage == this->T();

        unsigned int newsize = this->dimupto(laststage ? stage : stage+1);

        origvarsandconstraints(stage, xih, vars, constraints);

        if(stage)
        {
            if(!constraints)
                 constraints.reset(new linearconstraint_list);

            constraints->push_back(linearconstraint(newsize,linearconstraint::geq));
            constraints->push_back(linearconstraint(newsize,linearconstraint::geq));

            linearconstraint_list::reverse_iterator muc = constraints->rbegin();
            linearconstraint_list::reverse_iterator nuc = muc + 1;

            unsigned int i = this->stageoffest(stage)-1;
            muc->lhs[i] = this->mu();
            nuc->lhs[i++] = this->nu();

            linearfunction_ptr c;
            this->flp->get_f(stage,xih,c);

            for(unsigned int j = 0; j<this->flp->stagedim(stage); i++,j++)
            {
                muc->lhs[i] = -this->mu()*c->coefs[j];
                nuc->lhs[i] = -this->nu()*c->coefs[j];
            }

            muc->lhs[i] = 1; // theta
            nuc->lhs[i++] = 1;
            if(!laststage)
            {
                muc->lhs[i] = -this->mu(); // u_{k+1}
                nuc->lhs[i] = -this->nu();

                unsigned int nextthetapos = this->dimupto(stage+1)
                                              - (stage+1 == this->T() ? 1 : 2);
                muc->lhs[nextthetapos] = -this->mu();
                nuc->lhs[nextthetapos] = -this->nu();
            }
        }
    }
};

template<class Xi>
using mncvarproblem_ptr = std::shared_ptr<mncvarproblem<Xi>>;



#endif // CVAR_H
