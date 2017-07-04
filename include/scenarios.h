#ifndef SCENARIOS_H
#define SCENARIOS_H

#include "probability.h"


class nxtree: public indexedtree
{
public:
    nxtree(const std::vector<unsigned int>& ns)
           : fns(ns), fstageoffsets(ns.size()), fstagesizes(ns.size())
    {
        fnumnodes = 0;
        for(unsigned int i=0, s = 1; i<ns.size(); i++)
        {
            s *= ns[i];
            fstageoffsets[i] = fnumnodes;
            fstagesizes[i] = s;
            fnumnodes += s;
        }
    }
public:

    virtual unsigned int numbranches(const path& p) const
    {
        assert(p.size() < depth());
        return fns[p.size()];
    }

    virtual unsigned int depth() const
    {
        return fns.size();
    }
private:
    virtual unsigned int numnodes(unsigned int k) const
    {
        assert(k < depth());
        return fstagesizes[k];
    }

    virtual unsigned int stageoffset(unsigned int k) const
    {
        assert(k < depth());

        return fstageoffsets[k];
    }
    virtual unsigned int relnodeindex(const path& p, unsigned int k) const
    {
        assert(p.size() > 0);
        assert(p.size() <= depth());
        assert(k < depth());

        unsigned int a=0;

        unsigned int n=1;

        for(int i=k; i>=0; i--)
        {
            a += p[i] * n;
            n*=fns[i];
        }

        assert(a < fnumnodes);

        return a;
    }
protected:
    virtual unsigned int nodeindex(const path& p) const
    {
        unsigned int k=p.size()-1;
        return stageoffset(k)+relnodeindex(p,k);
    }

    virtual unsigned int totalnumnodes() const
    {
        return fnumnodes;
    }

private:
    std::vector<unsigned int> fns;
    std::vector<unsigned int> fstageoffsets;
    std::vector<unsigned int> fstagesizes;
    unsigned int fnumnodes;
};


class uniformtreeprobability: public treeprobability
{
public:
    /** Default destructor */
    uniformtreeprobability(const tree_ptr& t) : treeprobability(t) {}

    virtual prob operator()(const path& ap) const
    {
        assert(ap.size()>0);
        assert(ap.size()<=ft->depth());
        path p(ap);
        p.resize(ap.size()-1);
        return 1.0 / ft->numbranches(p);
    }

};

using uniformtreeprobability_ptr = std::shared_ptr<uniformtreeprobability>;


class homogeneoustree: public tree
{
public:

    homogeneoustree(unsigned int ad, unsigned int an )
           : fn(an), fd(ad)
    {
    }
public:

    virtual unsigned int numbranches(const path& p) const
    {
        return fn;
    }

    virtual unsigned int depth() const
    {
        return fd;
    }
private:
    unsigned int fn;
    unsigned int fd;
};

using homogeneoustree_ptr = std::shared_ptr<homogeneoustree>;

class iidtreeprobability: public treeprobability
{
public:
    iidtreeprobability(const homogeneoustree_ptr& t, const std::vector<prob>& aps)
        : treeprobability(t), fps(aps) {}

    virtual prob operator()(const path& ap) const
    {
        assert(ap.size()>0);
        assert(ap.size()<=ft->depth());
        assert(ap[ap.size()-1]<fps.size());
        return fps[ap[ap.size()-1]];
    }
private:
    std::vector<prob> fps;
};

using iidtreeprobability_ptr = std::shared_ptr<iidtreeprobability>;

template <typename Xi>
class iidtreemapping: public treemapping<Xi>
{
public:
    iidtreemapping(const homogeneoustree_ptr& t, const std::vector<Xi>& av) :
       treemapping<Xi>(t), fv(av) {}
    virtual Xi operator()(const path& p) const
    {
        assert(p.size()>0);
        assert(p.size()<=this->ft->depth());
        assert(p[p.size()-1]<fv.size());
        return fv[p[p.size()-1]];
    }
private:
    std::vector<Xi> fv;
};

template<typename Xi>
using iidtreemapping_ptr = std::shared_ptr<iidtreemapping<Xi>>;


template <typename Xi>
class idscenariotree : public scenariotree<Xi>
{
    static std::vector<unsigned int> nx(const std::vector<std::vector<Xi>>& axs )
    {
        std::vector<unsigned int> r;
        for(unsigned int i=0; i<axs.size(); i++)
            r.push_back(axs[i].size());
        return r;
    }

public:
    idscenariotree(const std::vector<std::vector<Xi>>& axs ) :
        ft(new nxtree(nx(axs))), fxs(axs), fprobs(axs.size())
    {
        for(unsigned int i=0; i<fxs.size();i++)
            fprobs[i]=1.0 / (double) fxs[i].size();
    }

    virtual unsigned int depth() const { return ft->depth(); }

    virtual Xi x(const path& p)  const
    {
        assert(p.size()>0);
        assert(p.size()<=this->ft->depth());
        assert(p[p.size()-1]<fxs[p.size()-1].size());
        return fxs[p.size()-1][p[p.size()-1]];
    }

    virtual scenario<Xi> s(const path& p) const
    {
        assert(p.size()>0);
        assert(p.size()<=this->ft->depth());
        assert(p[p.size()-1]<fxs[p.size()-1].size());

        scenario<Xi> r;
        for(unsigned int i=0; i<p.size(); i++)
        {
            r.push_back(fxs[i][p[i]]);
        }
        return r;
    }

    virtual prob p(const path& p) const
    {
        assert(p.size()>0);
        assert(p.size()<=this->ft->depth());
        assert(p[p.size()-1]<fxs[p.size()-1].size());
        return fprobs[p.size()-1];
    }
    virtual prob up(const path& p) const
    {
        assert(p.size()>0);
        assert(p.size()<=this->ft->depth());
        assert(p[p.size()-1]<fxs[p.size()-1].size());
        prob pr = 1.0;
        for(unsigned int i=0; i<p.size(); i++)
            pr *= this->fprobs[i];
        return pr;
    }

    virtual const tree_ptr& t() const { return ft; }

private:
    tree_ptr ft;
    std::vector<std::vector<Xi>> fxs;
    std::vector<prob> fprobs;
};

template<typename Xi>
using idscenariotree_ptr = std::shared_ptr<idscenariotree<Xi>>;


template <typename Xi>
class iidscenariotree : public scenariotree<Xi>
{
public:
    iidscenariotree(unsigned int depth, std::vector<Xi>& ax ) :
        ft(new homogeneoustree(depth,ax.size())), fx(ax),
         fprob(1.0 / (double) fx.size())
      {}

    virtual unsigned int depth() const { return ft->depth(); }

    virtual Xi x(const path& p)  const
    {
        assert(p.size()>0);
        assert(p.size()<=this->ft->depth());
        assert(p[p.size()-1]<fx.size());
        return fx[p[p.size()-1]];
    }

    virtual scenario<Xi> s(const path& p) const
    {
        assert(p.size()>0);
        assert(p.size()<=this->ft->depth());
        assert(p[p.size()-1]<fx.size());

        scenario<Xi> r;
        for(unsigned int i=0; i<p.size(); i++)
        {
            r.push_back(fx[p[i]]);
        }
        return r;
    }

    virtual prob p(const path& p) const
    {
        assert(p.size()>0);
        assert(p.size()<=this->ft->depth());
        assert(p[p.size()-1]<fx.size());
        return fprob;
    }
    virtual prob up(const path& p) const
    {
        assert(p.size()>0);
        assert(p.size()<=this->ft->depth());
        assert(p[p.size()-1]<fx.size());
        return pow(this->fprob,p.size());
    }

    virtual const tree_ptr& t() const { return ft; }

private:
    tree_ptr ft;
    std::vector<Xi>& fx;
    prob fprob;
};

template<typename Xi>
using iidscenariotree_ptr = std::shared_ptr<iidscenariotree<Xi>>;


class shiftedtree: public tree
{
public:
    shiftedtree(const tree_ptr& t): ft(t) {}
    virtual unsigned int depth() const { return ft->depth()+1; }

    static path orig(const path& src)
    {
        path dst(src.size()-1);
        for(unsigned int i=0; i<dst.size(); i++)
            dst[i] = src[i+1];
        return dst;
    }

    virtual unsigned int numbranches(const path& ap) const
    {
        if(ap.size()==0)
            return 1;
        else
            return ft->numbranches(orig(ap));
    }
private:
    const tree_ptr& ft;
};

using shiftedtree_ptr = std::shared_ptr<shiftedtree>;


template<typename Xi>
class shiftedscenariotree : public scenariotree<Xi>
{
public:
    shiftedscenariotree(const scenariotree_ptr<Xi>& orig, const Xi& x0 ) :
        forig(orig), fx0(x0), fst(new shiftedtree(orig->t()))
    {}

    virtual unsigned int depth() const
    {
        return fst->depth();
    }

    virtual Xi x(const path& p)  const
    {
        if(p.size()==1)
        {
            assert(p[0]==0);
            return fx0;
        }
        else
            return forig->x(shiftedtree::orig(p));
    }

    virtual scenario<Xi> s(const path& p) const
    {
        assert(p.size()>0);
        if(p.size()==1)
        {
            scenario<Xi> r;
            r.push_back(fx0);
            return r;
        }
        else
        {
            scenario<Xi> r = forig->s(shiftedtree::orig(p));
            r.insert(r.begin(),fx0);
            return r;
        }
    }

    virtual prob p(const path& p) const
    {
        if(p.size()==1)
        {
            assert(p[0]==0);
            return 1.0;
        }
        else
            return forig->p(shiftedtree::orig(p));
    }
    virtual prob up(const path& p) const
    {
        assert(p.size()>0);
        return p.size() == 1 ? 1.0 : forig->up(shiftedtree::orig(p));
    }

    virtual const tree_ptr& t() const { return fst; }

private:
    scenariotree_ptr<Xi> forig;
    tree_ptr fst;
    Xi fx0;
};

template<typename Xi>
using shiftedscenariotree_ptr = std::shared_ptr<shiftedscenariotree<Xi>>;


#endif // SCENARIOS_H

