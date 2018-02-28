#ifndef SCENARIOS_H
#define SCENARIOS_H

#include "mspp/commons.h"

namespace mspp
{

class indexedtree : public tree
{
public:
    virtual unsigned int totalnumnodes() const = 0;
    virtual unsigned int nodeindex(const path& p) const = 0;
};

using indexedtree_ptr = std::shared_ptr<indexedtree>;


template <typename Xi>
class generalprocess: public process<Xi>
{
public:
    generalprocess(const indexedtree_ptr& t) :
       process<Xi>(t), fdata(t->totalnumnodes()) {}
    virtual Xi operator()(const path& p) const
    {
        return fdata[it().nodeindex(p)];
    }
    virtual Xi& operator()(const path& p)
    {
        return fdata[it().nodeindex(p)];
    }
protected:
    const indexedtree& it() const
       { return dynamic_cast<indexedtree&> (*process<Xi>::ft); }
private:
    std::vector<Xi> fdata;
};

template<typename Xi>
using generalprocess_ptr = std::shared_ptr<generalprocess<Xi>>;


class nktree: public indexedtree
{
public:
    nktree(const std::vector<unsigned int>& ns)
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

    virtual unsigned int dupto(unsigned int k) const
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
        return dupto(k)+relnodeindex(p,k);
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

using nktree_ptr = std::shared_ptr<nktree>;


class uniformprobability: public probability
{
public:
    /** Default destructor */
    uniformprobability(const tree_ptr& t) : probability(t) {}

    virtual prob operator()(const path& ap) const
    {
        assert(ap.size()>0);
        assert(ap.size()<=ft->depth());
        path p(ap);
        p.resize(ap.size()-1);
        return 1.0 / ft->numbranches(p);
    }

};

using uniformprobability_ptr = std::shared_ptr<uniformprobability>;


class ntree: public tree
{
public:

    ntree(unsigned int ad, unsigned int an )
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

using ntree_ptr = std::shared_ptr<ntree>;


template <typename Xi>
class nkidprocess: public process<Xi>
{
    static tree_ptr maketree(const std::vector<std::vector<Xi>>& av)
    {
        std::vector<unsigned int> n;
        for(unsigned int i=0; i<av.size(); i++)
            n.push_back(av[i].size());
        return tree_ptr(new nktree(n));
    }

public:
    nkidprocess(const std::vector<std::vector<Xi>>& av) :
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
using nkidprocess_ptr = std::shared_ptr<nkidprocess<Xi>>;


template <typename Xi>
class iidprocess: public process<Xi>
{
public:
    iidprocess(const ntree_ptr& t, const std::vector<Xi>& av) :
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

class iidprobability: public probability
{
public:
    iidprobability(const ntree_ptr& t, const std::vector<prob>& aps)
        : probability(t), fps(aps) {}

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

using iidprobability_ptr = std::shared_ptr<iidprobability>;

}

#endif // SCENARIOS_H
