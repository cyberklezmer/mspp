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

    virtual unsigned int depth() const = 0;

    virtual Xi x(const path&) const = 0;
    virtual scenario<Xi> s(const path&) const = 0;
    virtual prob p(const path&) const = 0;
    virtual prob up(const path&) const = 0;

    virtual const tree_ptr& t() const = 0;
};

template<class Xi>
using scenariotree_ptr=std::shared_ptr<scenariotree<Xi>>;

template<class Xi>
class modularscenariotree : public scenariotree<Xi>
{
public:
    modularscenariotree(const treemapping_ptr<Xi>& x, const treeprobability_ptr& p)
         : fx(x), fp(p)
      {
          assert(fx->t().get()==fp->t().get());
      }
    virtual unsigned int depth() const
       { return fp->t()->depth(); }
    virtual Xi x(const path& p) const
       { return (*fx)(p); }
    virtual scenario<Xi> s(const path& p) const
       { return fx->s(p); }
    virtual prob p(const path& p) const
        { return (*fp)(p); }
    virtual prob up(const path& p) const
        { return fp->up(p); }
    const tree_ptr& t() const
        { return fp->t(); }
private:
    treemapping_ptr<Xi> fx;
    treeprobability_ptr fp;
};

template<class Xi>
using modularscenariotree_ptr=std::shared_ptr<modularscenariotree<Xi>>;


#endif // PROBABILITY_H
