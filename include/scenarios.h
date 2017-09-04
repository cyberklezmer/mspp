#ifndef SCENARIOS_H
#define SCENARIOS_H

#include "probability.h"

namespace mspp
{

class indexedtree : public tree
{
public:
    virtual unsigned int totalnumnodes() const = 0;
    virtual unsigned int nodeindex(const path& p) const = 0;
};

using indexedtree_ptr = std::shared_ptr<indexedtree>;

template <typename L>
class generalprocess: public process<L>
{
public:
    generalprocess(const indexedtree_ptr& t) :
       process<L>(t), fdata(t->totalnumnodes()) {}
    virtual L operator()(const path& p) const
    {
        return fdata[it().nodeindex(p)];
    }
    virtual L& operator()(const path& p)
    {
        return fdata[it().nodeindex(p)];
    }
protected:
    const indexedtree& it() const
       { return dynamic_cast<indexedtree&> (*process<L>::ft); }
private:
    std::vector<L> fdata;
};

template<typename L>
using generalprocess_ptr = std::shared_ptr<generalprocess<L>>;



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

using nxtree_ptr = std::shared_ptr<nxtree>;


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
class idprocess: public process<Xi>
{
    static tree_ptr maketree(const std::vector<std::vector<Xi>>& av)
    {
        std::vector<unsigned int> n;
        for(unsigned int i=0; i<av.size(); i++)
            n.push_back(av[i].size());
        return tree_ptr(new nxtree(n));
    }

public:
    idprocess(const std::vector<std::vector<Xi>>& av) :
       process<Xi>(maketree(av)), fv(av) {}
    virtual Xi operator()(const path& p) const
    {
        assert(p.size()>0);
        assert(p.size()<=this->ft->depth());
        assert(p[p.size()-1]<fv[p.size()-1].size());
        return fv[p.size()-1][p[p.size()-1]];
    }
private:
    std::vector<std::vector<Xi>> fv;
};

template<typename Xi>
using idprocess_ptr = std::shared_ptr<idprocess<Xi>>;


template <typename Xi>
class iidprocess: public process<Xi>
{
public:
    iidprocess(const homogeneoustree_ptr& t, const std::vector<Xi>& av) :
       process<Xi>(t), fv(av) {}
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
using iidprocess_ptr = std::shared_ptr<iidprocess<Xi>>;


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

template <typename L>
class shiftedprocess : public process<L>
{
public:
    shiftedprocess(const process_ptr<L>& orig, const L& x0) :
        process<L>(tree_ptr(new shiftedtree(orig->t()))),
        forig(orig), fx0(x0)
    {}

    virtual L operator()(const path& p) const
     { return p.size()==1 ? fx0 : (*forig)(shiftedtree::orig(p)) ; }
protected:
    process_ptr<L> forig;
    L fx0;
};


class shiftedtreeprobability : public treeprobability
{
public:
    shiftedtreeprobability(const treeprobability_ptr& orig) :
        treeprobability(tree_ptr(new shiftedtree(orig->t()))), forig(orig)
    {}

    virtual prob operator()(const path& p) const
     { return p.size()==1 ? 1.0 : (*forig)(shiftedtree::orig(p)) ; }

protected:
    treeprobability_ptr forig;
};



template<typename Xi>
class shiftedscenariotree : public modularscenariotree<Xi>
{
public:
    shiftedscenariotree(const process_ptr<Xi>& x,
                        const treeprobability_ptr& p,
                        const Xi& x0 ) :
        modularscenariotree<Xi>(
            process_ptr<Xi>(new shiftedprocess<Xi>(x,x0)),
            treeprobability_ptr(new shiftedtreeprobability(p)))
    {}
};

template<typename Xi>
using shiftedscenariotree_ptr = std::shared_ptr<shiftedscenariotree<Xi>>;

template<typename Xi>
inline modularscenariotree_ptr<Xi>
      idscenariotree(std::vector<std::vector<Xi>>& dist)
{
      std::vector<unsigned int> n;
      for(unsigned int i=0; i<dist.size(); i++)
          n.push_back(dist[i].size());
      nxtree_ptr t(new nxtree(dist));
      uniformtreeprobability_ptr p(new uniformtreeprobability(t) );
      idprocess_ptr<Xi> m(p,idprocess<Xi>(dist));
      return new modularscenariotree<Xi>(m,p);
}

}

#endif // SCENARIOS_H

