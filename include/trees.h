#ifndef TREES_H
#define TREES_H

#include "commons.h"

using path = std::vector<unsigned int>;

class tree;

class treecallback
{
public:
    virtual void callback(const path& p) = 0;
};

class tree : public object
{
public:
    virtual unsigned int depth() const = 0;
    virtual unsigned int numbranches(const path& p) const = 0;
    void foreachnode(treecallback *callee)
    {
         if(depth()>0)
         {
             path p(0);
             doforeachnode(callee,p);
         }
     }

private:
    void doforeachnode(
             treecallback *callee,
             path& ap
            )
    {
        unsigned int k=ap.size();
        unsigned int n=numbranches(ap);

        path p(ap);
        p.push_back(0);
        bool cb = k+1 < depth();
        for(p[k]=0; p[k]<n; p[k]++)
        {
            callee->callback(p);
            if(cb)
                doforeachnode(callee, p);
        }
    }
};

using tree_ptr = std::shared_ptr<tree>;


class indexedtree : public tree
{
public:
    virtual unsigned int totalnumnodes() const = 0;
    virtual unsigned int nodeindex(const path& p) const = 0;
};

using indexedtree_ptr = std::shared_ptr<indexedtree>;


template <typename L>
using scenario = std::vector<L>;

template <typename L>
class treemapping : public object
{
public:
    treemapping(const tree_ptr& t) : ft(t) {}

    // tbd possibly make a reference
    virtual L operator()(const path& p) const = 0;
    virtual scenario<L> s(const path& p) const
    {
        scenario<L> r;
        path np;
        for(unsigned int i=0; i<p.size(); i++)
        {
            np.push_back(p[i]);
            r.push_back((*this)(np));
        }
        return r;
    }
    const tree_ptr& t() const { return ft; }
protected:
    tree_ptr ft;
};

template <typename L>
using treemapping_ptr = std::shared_ptr<treemapping<L>>;


template <typename L>
class generaltreemapping: public treemapping<L>
{
public:
    generaltreemapping(const indexedtree_ptr& t) :
       treemapping<L>(t), fdata(t->totalnumnodes()) {}
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
       { return dynamic_cast<indexedtree&> (*treemapping<L>::ft); }
private:
    std::vector<L> fdata;
};

template<typename L>
using generaltreemapping_ptr = std::shared_ptr<generaltreemapping<L>>;



#endif // TREES_H
