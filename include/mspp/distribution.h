#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include "mspp/commons.h"

namespace mspp
{

/// \addtogroup Distributions
/// @{

using probability = double;

template <typename X, typename C>
class distribution: public object
{
public:
    C c(const scenario<X> s) const { return C(s); }
};

template <typename X, typename C, typename P>
class parametricdistribution: virtual public distribution<X,C>
{
public:
    parametricdistribution(const P& p): fp(p) {}
protected:
    const P& p() { return fp; }
private:
    P fp;
};

template <typename X, typename D>
class processdistribution : public object
{
public:
    processdistribution(const X& z0, D& d, unsigned int T)
      : fz0(z0), fd(T)
    {
        for(unsigned int i=0; i<T; i++)
            fd[i]=d;
    }
    processdistribution(const X& z0, const std::vector<D>& d)
      : fz0(z0), fd(d)
    {
    }

    X z0() const { return fz0; }
    const D& d(unsigned int k) const { assert(k>0 && k<fd.size()); return fd[k]; }
    unsigned int T() const { return fd.size(); }
private:
    X fz0;
    std::vector<D> fd;
};


/// \addtogroup mcdisc MC distributions
/// \ingroup Distributions
/// @{

template <typename X, typename C>
class mcdistribution: virtual public distribution<X, C>
{
public:
    X draw(const C& c) const
    {
        return do_draw(c);
    }

    X draw(const scenario<X>& s) const
    {
        return do_draw(distribution<X,C>::c(s));
    }
private:
    virtual X do_draw(const C& c) const = 0;
};



/// @}

/// \addtogroup discretedisc Discrete distributions
/// \ingroup Distributions
/// @{

template <typename X>
struct atom
{
    X x;probability p;
};

template <typename X, typename C>
class ddistribution: virtual public distribution<X,C>
{
public:
    ddistribution() {}

    void atoms(std::vector<atom<X>>& a, const C& c) const
    {
        a.resize(0);
        atoms_are(a,c);
    }
    void atoms(std::vector<atom<X>>& a, const scenario<X>& s) const
    {
        a.resize(0);
        atoms_are(a,distribution<X,C>::c(s));
    }

protected:
    virtual void atoms_are(std::vector<atom<X>>& a, const C& s) const = 0;
};



template <typename X>
struct indexedatom {X x;probability p; unsigned int i; };

template <typename X>
class indexedhistory : public std::vector<X>
{
public:
    double uncprob() const
    {
        assert(std::vector<X>::size());
        assert(*this[0].p==1);
        double r = 1.0;
        for(unsigned int i=1; i<std::vector<X>::size(); i++)
            r *= *this[i].p;
        return r;
    }

    scenario<X> s() const
    {
        scenario<X> r(std::vector<X>::size());
        for(unsigned int i=0; i<std::vector<X>::size(); i++)
            r[i]=*this[i].x;
        return r;
    }
};


template<typename X>
class sctreecallback
{
public:
    virtual void callback(const indexedhistory<X>& s) = 0;
};

template <typename X>
class scenariotree: public object
{
private:
    virtual unsigned int T_is() const = 0;
    virtual X root_is() const = 0;
    virtual void branches_are(std::vector<atom<X>>& bchs, const indexedhistory<X>& s) const = 0;
public:
    unsigned int T() { return T_is(); }
    void foreachnode(sctreecallback<X> *callee) const
    {
        indexedhistory<X> s;
        s.push_back({root_is(),1,0});
        callee->callback(s);
        if(T_is())
            doforeachnode(callee, s);
    }
private:
    void doforeachnode(sctreecallback<X> *callee, indexedhistory<X> s) const
    {
        unsigned int k=s.size();
        std::vector<atom<X>> itms;
        branches_are(itms,s);
        s.resize(k+1);
        bool cb = k < this->T();
        unsigned int i=0;
        for(typename std::vector<atom<X>>::iterator it = itms.begin();
              it != itms.end();
              it++,i++)
        {
            s[k]={it->x,it->p, i};
            callee->callback(s);
            if(cb)
                doforeachnode(callee, s);
        }
    }
};


template <typename X, typename D>
class distrscenariotree : public scenariotree<X>
{
public:
    distrscenariotree(const processdistribution<X,D>& p):fp(p)
    {}
private:
    virtual void branches_are(std::vector<atom<X>>& bchs,
                              const indexedhistory<X>& s) const
    {
        unsigned int k=s.size();
        assert(k);
        assert(fp.T());
        fp.d(k).atoms(bchs,s.s());
    }
    virtual X root_is() const
    {
        return fp.z0();
    }
    virtual unsigned int  T_is() const { return fp.T(); }
private:
    processdistribution<X,D> fp;
};

/// @}


/// @} - distributions


} // namespace


#endif // DISTRIBUTION_H
