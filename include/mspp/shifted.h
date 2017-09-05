#ifndef SHIFTED_H
#define SHIFTED_H

#include "mspp/commons.h"

namespace mspp
{

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


class shiftedprobability : public probability
{
public:
    shiftedprobability(const probability_ptr& orig) :
        probability(tree_ptr(new shiftedtree(orig->t()))), forig(orig)
    {}

    virtual prob operator()(const path& p) const
     { return p.size()==1 ? 1.0 : (*forig)(shiftedtree::orig(p)) ; }

protected:
    probability_ptr forig;
};


template<typename Xi>
class shiftedscenariotree : public modularscenariotree<Xi>
{
public:
    shiftedscenariotree(const process_ptr<Xi>& x,
                        const probability_ptr& p,
                        const Xi& x0 ) :
        modularscenariotree<Xi>(
            process_ptr<Xi>(new shiftedprocess<Xi>(x,x0)),
            probability_ptr(new shiftedprobability(p)))
    {}
};

template<typename Xi>
using shiftedscenariotree_ptr = std::shared_ptr<shiftedscenariotree<Xi>>;

}

#endif // SHIFTED_H
