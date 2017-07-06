#ifndef MMPCVAR_H
#define MMPCVAR_H

#include "linear.h"

template<class Xi>
class mmpcvarproblem: public linearproblem<Xi>
{
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

public:
    mmpcvarproblem(const linearproblem_ptr<Xi>& lp, double alpha, double lambda) :
        linearproblem<Xi>(dims(lp->stagedims())), // tbd ordering
        flp(lp), falpha(alpha), flambda(lambda)
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
            flp->get_f(0,xih,orig);
            for(unsigned int i=0; i<k-1; i++)
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

        flp->get_constraints(stage,xih,vars,constraints);

        if(stage >0)
            vars->push_back(varinfo(varinfo::R,"theta"));
        if(!laststage)
            vars->push_back(varinfo(varinfo::R,"u"));

        if(constraints)
        {
            for(linearconstraint_list::iterator p = constraints->begin();
                 p != constraints->end(); p++)
            {
//                    std::cout << "sd: " << flp->dimupto(stage) << std::endl;

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

        if(stage)
        {
            if(!constraints)
                 constraints.reset(new linearconstraint_list);
            linearfunction_ptr c;
            flp->get_f(stage,xih,c);

            constraints->push_back(linearconstraint(newsize,linearconstraint::geq));
            constraints->push_back(linearconstraint(newsize,linearconstraint::geq));
            linearconstraint_list::reverse_iterator muc = constraints->rbegin();
            linearconstraint_list::reverse_iterator nuc = muc + 1;

            unsigned int i = newsize-(laststage ? 1 : 2);
            muc->lhs[i] = 1;
            nuc->lhs[i--] = 1;

            for(int j=flp->stagedim(stage)-1; j>=0; j--)
            {
                muc->lhs[i] = -mu() * c->coefs[j];
                nuc->lhs[i--] = -nu() * c->coefs[j];
            }

            muc->lhs[i] = mu();
            nuc->lhs[i] = nu();
        }
    }
    double mu() const { return 1.0 - flambda; }
    double nu() const { return 1.0 - flambda + flambda / falpha; }
public:
    linearproblem_ptr<Xi> flp;
    double falpha;
    double flambda;
};

template<class Xi>
using mmpcvarproblem_ptr = std::shared_ptr<mmpcvarproblem<Xi>>;


#endif // MMPCVAR_H