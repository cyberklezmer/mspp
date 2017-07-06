#ifndef PROBABILITY_H
#define PROBABILITY_H

#include "trees.h"

using prob = double;

class treeprobability: public treemapping<prob>
{
public:
    treeprobability(const tree_ptr& t) : treemapping(t) {}
    virtual prob up(const path& ap) const
    {
        prob pr=1.0;
        path p;
        for(unsigned int i=0; i<ap.size(); i++)
        {
            p.push_back(ap[i]);
            pr*=(*this)(p);
        }
        return pr;
    }
};

using treeprobability_ptr = std::shared_ptr<treeprobability>;


template <typename Xi>
class distribution : public object
{
};



template <typename Xi>
class mcdistribution : virtual public distribution<Xi>
{
public:
    virtual Xi draw() = 0;
};

template <typename Xi>
class discretizabledistribution : virtual public distribution<Xi>
{
public:
    virtual std::vector<Xi> d( unsigned int n) = 0;
};

template <typename Xi>
class scenariotree : public object
{
public:
    scenariotree()
      {}

    virtual const treemapping_ptr<Xi>& tm() const = 0;
    virtual const treeprobability_ptr& tp() const = 0;

    const tree_ptr& t() const { return tp()->t(); }
    virtual unsigned int depth() const { return tp()->t()->depth();}
    Xi x(const path& p) const  { return (*tm())(p); }
    scenario<Xi> s(const path& p) const { return tm()->s(p); }
    virtual prob p(const path& p) const { return (*tp())(p); }
    virtual prob up(const path& p) const { return tp()->up(p); }
};

template<class Xi>
using scenariotree_ptr=std::shared_ptr<scenariotree<Xi>>;



template<class Xi>
class modularscenariotree : public scenariotree<Xi>
{
public:
    modularscenariotree(const treemapping_ptr<Xi>& x,
                        const treeprobability_ptr& p)
         : fx(x), fp(p)
    {
    }

    virtual const treemapping_ptr<Xi>& tm() const { return fx; }
    virtual const treeprobability_ptr& tp() const { return fp; }

protected:
    treemapping_ptr<Xi> fx;
    treeprobability_ptr fp;
};

template<class Xi>
using modularscenariotree_ptr=std::shared_ptr<modularscenariotree<Xi>>;

#endif // PROBABILITY_H
