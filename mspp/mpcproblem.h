#ifndef MPCPROBLEM_H
#define MPCPROBLEM_H

#include "mspp/msproblem.h"


namespace mspp
{

/// Adding the previous "u"
template <typename R>
class mpcvprestriction: public tilde
{
private:
    virtual bool is_included(unsigned int i, unsigned int s) const
    {
        if(i==s-1)
            return true;
        else
            return R().included(i,s);
    }
    virtual bool is_included(unsigned int i, unsigned int j, unsigned int s) const
    {
        unsigned int varsadded = i==0 ? 1 : 2;  // last stage souhol not come here
        bool uincluded = i==s-1;
        if(j==0)
            return uincluded;
        else if(j==1 && varsadded == 2)
            return false;
        else
            return R().included(i,j-varsadded,s);
    }
};


template<class P>
class mmpcvarequivalent: public msproblem
               <expectation, linearfunction, typename P::G_t, typename P::V_t,
                           typename P::X_t,
                           mpcvprestriction<typename P::Rx_t>,
                           typename P::Rxi_t>
{
    static msproblemstructure
        ps(const msproblemstructure& sps)
    {
        assert(sps.size()>0);
        msproblemstructure r(sps.size());
        r[0]=sps[0]+1;
        int i=1;
        for(; i< sps.size()-1; i++)
            r[i]=sps[i]+2;
        r[i] = sps[i]+1;
        return r;
    }

public:

    mmpcvarequivalent(const P& sp) :
        msproblem<expectation, linearfunction,
             typename P::G_t, typename P::V_t,
             typename P::X_t,
             mpcvprestriction<typename P::Rx_t>,
             typename P::Rxi_t>(ps(sp.d)),
           fsp(new P(sp)),
           flambda( sp.rho().lambda),
           falpha( sp.rho().alpha)
    {
        static_assert(std::is_same<typename P::C_t,mmpcvar>::value);
        assert(fsp->T()>0);
        assert(falpha > 0);
    }

    mmpcvarequivalent(const P& sp,double lambda, double alpha) :
        msproblem<expectation, linearfunction,
             typename P::G_t, typename P::V_t,
             typename P::X_t,
             mpcvprestriction<typename P::Rx_t>,
             typename P::Rxi_t>(ps(sp.d)),
           fsp(new P(sp)),
           flambda(lambda),falpha(alpha)
    {
        static_assert(std::is_same<typename P::R_t,expectation>::value);
        assert(fsp->T()>0);
        assert(falpha > 0);
    }

    virtual std::string varname_is(unsigned int stage, unsigned int i) const
    {
        if(stage == 0)
        {
            if(i==0)
                return "u";
            else
                return fsp->varname(stage,i-1);
        }
        else if(stage==this->T())
        {
            if(i==0)
                return "theta";
            else
                return fsp->varname(stage,i-1);
        }
        else
        {
            if(i==0)
                return "u";
            else if(i==1)
                return "theta";
            else
                return fsp->varname(stage,i-2);
        }
    }

    virtual linearfunction f_is(
            unsigned int k,
            const vectors<typename P::X_t>& barxi) const
    {
        linearfunction r = this->newf(k);
        if(k==0)
        {
            linearfunction orig = fsp->f_tbd(0,barxi);

            assert(orig.xdim()+1==r.xdim());

            unsigned int i=0;
            r.setc(i++,1);
            for(;i-1<orig.xdim(); i++)
                r.setc(i,orig.c(i-1));
        }
        else if(k < this->T())
        {
            assert(this->barxdim(k)>=this->tildexdim(k)+2);
            r.setc(this->tildexdim(k),1);
            r.setc(this->tildexdim(k)+1,1);
        }
        else // k==T
        {
            assert(this->barxdim(k)>=this->tildexdim(k)+1);
            r.setc(this->tildexdim(k),1);
        }
        return r;
    }

    virtual void xset_is(
            unsigned int k,
            const vectors<typename P::X_t>& barxi,
            ranges<typename P::V_t>& r,
            msconstraints<typename P::G_t>& g
            ) const
    {
        using G_t = typename P::G_t;
        using V_t = typename P::V_t;

        ranges<typename P::V_t> srcr;
        msconstraints<G_t> srcgs;
        fsp->xset_tbd(k,barxi,srcr,srcgs);


        unsigned int m; // # of new variables in this stage

        if(k==0 || k==this->T())
            m=1;
        else
            m=2;

        unsigned int dst=0;
        for(;dst<m;)
            r[dst++]=range<V_t>(range<V_t>::realt);
        for(unsigned int src=0; src<srcr.size(); )
            r[dst++] = srcr[src++];

        for(typename msconstraints<G_t>::iterator it = srcgs.begin();
             it != srcgs.end(); it++)
        {
            const linearmsconstraint& srcg = *it;

            linearmsconstraint& dstg = this->addg(g,k);

            unsigned int src = 0;
            unsigned int dst = 0;

            if(k>1) // we copy all coefs
            {
                assert(this->barxdimupto(k-2,k)==this->fsp->barxdimupto(k-2,k));
                for(;src<this->barxdimupto(k-2,k);)
                    dstg.setlhs(dst++,srcg.lhs(src++));
            }
            if(k>0) // here we make space for previous stage u
            {
                assert(this->barxdimupto(k-1,k)==fsp->barxdimupto(k-1,k)+1);
                dstg.setlhs(dst++,0);
                for(;src<this->fsp->barxdimupto(k-1,k);)
                    dstg.setlhs(dst++,srcg.lhs(src++));
            }
            // now we make space for m vars
            for(unsigned int i=0; i<m; i++)
                dstg.setlhs(dst++,0);
            for(;dst<this->barxdim(k);)
                dstg.setlhs(dst++,srcg.lhs(src++));

            assert(src==this->fsp->barxdim(k));
            assert(dst==this->barxdim(k));

            dstg.setrhs(srcg.rhs());
            dstg.settype(srcg.t());
        }

        // next we add "mu" and "nu" cosntraints
        if(k>0)
        {
            typename P::F_t srcf = fsp->f_tbd(k,barxi);
            assert(srcf.xdim() == this->fsp->barxdim(k));

            G_t muc(this->barxdim(k));
            G_t nuc(this->barxdim(k));

            unsigned int src = 0;
            unsigned int dst = 0;

            if(k>1) // we add coefs up to stage k-2
            {
                for(;src<this->fsp->barxdimupto(k-2,k);)
                {
                    muc.setlhs(dst,-mu() * srcf.c(src));
                    nuc.setlhs(dst++,-nu() * srcf.c(src++));
                }
            }

            // next we add a coef for u
            muc.setlhs(dst,mu());
            nuc.setlhs(dst++,nu());

            for(;src<this->fsp->barxdimupto(k-1,k);)
            {
                muc.setlhs(dst,-mu() * srcf.c(src));
                nuc.setlhs(dst++,-nu() * srcf.c(src++));
            }

            if(k<this->T())
            {
                // skip u in this stage
                muc.setlhs(dst,0);
                nuc.setlhs(dst++,0);
            }

            // add coef for this stage theta
            muc.setlhs(dst,1);
            nuc.setlhs(dst++,1);

            // add remaining coefs
            for(; src< this->fsp->barxdim(k);)
            {
                muc.setlhs(dst,-mu() * srcf.c(src));
                nuc.setlhs(dst++,-nu() * srcf.c(src++));
            }

            assert(src==this->fsp->barxdim(k));
            assert(dst==this->barxdim(k));

            muc.settype(constraint::geq);
            muc.setrhs(0);
            nuc.settype(constraint::geq);
            nuc.setrhs(0);

            g.add(muc);
            g.add(nuc);
        }
    }
    double mu() const { return 1.0 - flambda; }
    double nu() const { return 1.0 - flambda + flambda / falpha; }
public:
    ptr<P> fsp;
    double flambda;
    double falpha;
};


}


#endif // MPCPROBLEM_H
