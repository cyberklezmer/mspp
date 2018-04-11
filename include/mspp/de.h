#ifndef DE_H
#define DE_H

#include "mspp/problem.h"
#include "mspp/distribution.h"
#include "mspp/solver.h"

namespace mspp
{

/// \addtogroup de Deterministic equivalent
/// \ingroup sms
/// @{

template <typename S>
class desolution: public object, public sctreecallback<variables>
{    
public:
    desolution(const S& s, const msproblemstructure& ps) :
        fs(ptr<S>(new S(s))), fps(ps)
       {}
    desolution(const ptr<S> s, const msproblemstructure& ps) :
        fs(s), fps(ps)
       {}
    void set(const variables& x) { fx=ptr<variables>(new variables); }
    void set(const ptr<variables> x) { fx=x; }
    void foreach()
    {
        fs->foreach(this);
    }
    virtual void callback(const indexedhistory<variables>&) {};
private:
    ptr<variables> fx;
    ptr<S> fs;
    msproblemstructure fps;
};

template<typename P,typename S,typename Xi>
class demethod : public object, public sctreecallback<Xi>
{
public:
    static bool solve(
             const P& p,
             const S& z,
             const lpsolver& lps,
             double& optimal,
             desolution<S>& sol)
            {/*
                foffsets.resize(this->fp->T()+1);
                fdim = 0;
                fobj.clear();
                fvars.clear();
                fconstraints.clear();
                fvarnames.clear();
                assert(this->fz->T() == this->fp->T());

                fsrcexpconstraints.clear();
                fdstexpconstraints.clear();
                fups.clear();


                this->fz->foreachnode(this);
                std::vector<double> x(fdim);

                fs.solve(fvars,fobj,fconstraints,fvarnames,x,optimal);
                sol->assign(x);*/
                return true;
            }

private:
    unsigned int fdim;
    std::vector<double> fobj;
    vardefs<realvar> fvars;
    std::vector<sparselinearconstraint_ptr> fconstraints;
    std::vector<std::string> fvarnames;

    // stored expectation msconstraints
    std::vector<linearconstraint> fsrcexpconstraints;
    std::vector<sparselinearconstraint*> fdstexpconstraints;
    std::vector<probability> fups;


    // list
    unsigned int fvindex;

public:
    virtual void callback(const indexedhistory<Xi>& a)
    {
/*        unsigned int stage = a.size()-1;
        probability up = a.uncprob();
        scenario<Xi> xi = a.s();

        unsigned int thisstagedim = this->fp->d[stage];

        // calling original problems \p msconstraints

        ptr<vardefs<realvar>> vars;
        ptr<msconstraints<linearconstraint>> msconstraints;
        demethod::fp->msconstraints(stage,xi,vars,msconstraints);

        ptr<linearobjective> f;
        this->fp->f(stage,xi,f);

        foffsets[stage] = fdim;
        fobj.resize(fdim + thisstagedim);

        // this stage variables and f

        for(unsigned int i=0; i<thisstagedim; i++)
        {
            fobj[fdim] = up * (*f)[i];
            const realvar& v = (*vars)[i];
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

        // stored msconstraints with expectations
        for(unsigned int i=0; i<fsrcexpconstraints.size(); i++)
        {
            assert(i<fdstexpconstraints.size());

            for(unsigned int src=this->fp->d.upto(stage), dst=foffsets[stage];
                src<this->fp->d.sum(stage); src++,dst++)
            {
                if(src < fsrcexpconstraints[i].lhssize())
                {
                     double coef = fsrcexpconstraints[i].lhs(src);
                     if(coef)
                        fdstexpconstraints[i]->set_lhs(dst,coef*up/fups[i]);
                }
            }

        }

        // this stage msconstraints
        if(msconstraints)
            for(unsigned int j=0; j<msconstraints->size(); j++)
            {
                const linearconstraint& s = *(*msconstraints)[j];

                sparselinearconstraint_ptr d(new sparselinearconstraint());
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
                           d->set_lhs(dst,v);
                        dst++;
                    }
                    d->rhs = s.rhs();
                    d->t = s.t();
                }
                if(s.lhssize() > this->fp->d.sum(stage))
                {
                    assert(stage < this->fp->T());
                    fsrcexpconstraints.push_back(s);
                    fdstexpconstraints.push_back(d.get());
                    fups.push_back(up);
                }
                fconstraints.push_back(d);
            }*/
    }

private:
    std::vector<unsigned int> foffsets;
};

/// @}


} // namespace

#endif // DE_H
