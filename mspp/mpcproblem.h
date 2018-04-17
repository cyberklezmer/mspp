#ifndef MPCPROBLEM_H
#define MPCPROBLEM_H

#include "mspp/msproblem.h"

namespace mspp
{

template<class S>
class mmpcvarproblem: public msproblem
               <typename S::V_t, linearfunction, typename S::G_t, typename S::C_t>
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
    mmpcvarproblem(const S& sp, double alpha, double lambda) :
        msproblem<typename S::V_t,
                  linearfunction,
                  typename S::G_t,
                  typename S::C_t>(ps(sp.d)),
           fsp(new S(sp)), falpha(alpha), flambda(lambda)
    {
        assert(fsp->T()>0);
        assert(falpha > 0);
    }

    virtual std::string varname_is(unsigned int stage, unsigned int i) const
    {
        if(i<fsp->d[stage])
            return fsp->varname(stage,i);
        if(stage == 0)
            return "u";
        if(stage == this->T())
            return "theta";
        else if(i==fsp->d[stage])
            return "theta";
        else
            return "u";
    }

    virtual linearfunction f_is(
            unsigned int k,
            const typename S::C_t& barxi) const
    {
        unsigned int dk = this->d[k];
        linearfunction r(std::vector<double>(dk,0));
        if(!k)
        {
            typename S::F_t orig = fsp->f(0,barxi);

            for(unsigned int i=0;i<dk-1; i++)
                r.setc(i,orig.c(i));
        }
        r.setc(dk-1,1);
        if(k > 0 && k < this->T())
            r.setc(dk-2,1);
        return r;
    }

    virtual void stageinfo_is(
            unsigned int k,
            const typename S::C_t& barxi,
            vardefs<typename S::V_t>& r,
            msconstraints<typename S::G_t>& g
            ) const
    {
        bool laststage = k == this->T();

        unsigned int newsize = this->d.sum(k);

        vardefs<typename S::V_t> srcr;

        using G_t = typename S::G_t;

        msconstraints<G_t> srcgs;
        fsp->stageinfo(k,barxi,srcr,srcgs);

        unsigned int i=0;
        for(; i<srcr.size(); i++)
            r[i] = srcr[i];

        if(k>0)
            r[i++]=realvar(realvar::R);
        if(!laststage)
            r[i]=realvar(realvar::R);


        for(typename msconstraints<G_t>::iterator srcg = srcgs.begin();
             srcg != srcgs.end(); srcg++)
        {
//                    std::cout << "sd: " << flp->sumd(stage) << std::endl;
            G_t& dstg = this->addg(g,k);
            dstg.setrhs(srcg->rhs());
            dstg.settype(srcg->t());

            unsigned int src = 0;
            unsigned int dst = 0;

            for(unsigned int i=0; i<=k; i++)
            {
                for(unsigned int j=0; j<fsp->d[i]; j++)
                    dstg.setlhs(dst++,srcg->lhs(src++));
                if(i<this->T())
                    dstg.setlhs(dst++,0); // u
                if(i > 0)
                    dstg.setlhs(dst++,0); // theta
            }
        }

        if(k)
        {
            typename S::F_t c = fsp->f(k,barxi);

            G_t& muc = this->addg(g,k);
            G_t& nuc = this->addg(g,k);
            muc.settype(constraint::geq);
            muc.setrhs(0);
            nuc.settype(constraint::geq);
            nuc.setrhs(0);

            unsigned int i = newsize-(laststage ? 1 : 2);
            muc.setlhs(i,1);
            nuc.setlhs(i--,1);

            for(int j=fsp->d[k]-1; j>=0; j--)
            {
                muc.setlhs(i,-mu() * c.c(j));
                nuc.setlhs(i--,-nu() * c.c(j));
            }

            muc.setlhs(i,mu());
            nuc.setlhs(i,nu());
        }
    }
    double mu() const { return 1.0 - flambda; }
    double nu() const { return 1.0 - flambda + flambda / falpha; }
public:
    ptr<S> fsp;
    double falpha;
    double flambda;
};


}


#endif // MPCPROBLEM_H
