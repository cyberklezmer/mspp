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
    using X_t = X;
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

template <typename X>
class imcdistribution: public mcdistribution<X, emptycondition<X>>
{
public:
    X draw() const
    {
        do_draw();
    }
private:
    virtual X do_draw(const emptycondition<X>&) const
    {
        return do_draw();
    }
    virtual X do_draw() const = 0;
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
        atoms(a,distribution<X,C>::c(s));
    }

protected:
    virtual void atoms_are(std::vector<atom<X>>& a, const C& s) const = 0;
};

template <typename X>
class iddistribution: virtual public ddistribution<X,emptycondition<X>>
{
public:
    iddistribution() {}

    void atoms(std::vector<atom<X>>& a) const
    {
        a.clear();
        atoms_are(a);
    }
    void atoms(std::vector<atom<X>>& a, const emptycondition<X>&) const
    {
        atoms(a);
    }

    void atoms(std::vector<atom<X>>& a, const scenario<X>&) const
    {
        atoms(a);
    }
protected:
    virtual void atoms_are(std::vector<atom<X>>& a, const emptycondition<X>&) const
      { atoms_are(a); }
    virtual void atoms_are(std::vector<atom<X>>& a) const = 0;
};

/// @} - discrete distributions


/// @} - distribution

/// \addtogroup processes Proceses
/// @{

template <typename D>
class processdistribution : public object
{
public:
    processdistribution(const typename D::X_t& z0, const D& d, unsigned int T)
      : fz0(z0)
    {
        for(unsigned int i=0; i<T; i++)
            fd.push_back(d);
    }
    processdistribution(const typename D::X_t& z0, const std::vector<D>& d)
      : fz0(z0), fd(d)
    {
    }

    typename D::X_t z0() const { return fz0; }
    const D& d(unsigned int k) const { assert(k>0 && k<=fd.size()); return fd[k-1]; }
    unsigned int T() const { return fd.size(); }
private:
    typename  D::X_t fz0;
    std::vector<D> fd;
};


template <typename X>
struct indexedatom {X x;probability p; unsigned int i; };

template <typename X>
class indexedhistory : public std::vector<indexedatom<X>>
{
public:
    double uncprob() const
    {
        assert(std::vector<indexedatom<X>>::size());
        assert((*this)[0].p==1);
        double r = 1.0;
        for(unsigned int i=1; i<std::vector<indexedatom<X>>::size(); i++)
            r *= (*this)[i].p;
        return r;
    }

    scenario<X> s() const
    {
        scenario<X> r(std::vector<indexedatom<X>>::size());
        for(unsigned int i=0; i<std::vector<indexedatom<X>>::size(); i++)
            r[i]=(*this)[i].x;
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
public:
    using X_t = X;
public:
    unsigned int T() const { return T_is(); }
    void foreachnode(sctreecallback<X> *callee) const
    {
        indexedhistory<X> s;
        s.push_back({root_is(),1,0});
        callee->callback(s);
        if(T_is())
            doforeachnode(callee, s);
    }
private:
    virtual unsigned int T_is() const = 0;
    virtual X root_is() const = 0;
    virtual void branches_are(std::vector<atom<X>>& bchs, const indexedhistory<X>& s) const = 0;

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


template <typename D>
class distrscenariotree : public scenariotree<typename D::X_t>
{
    using X=typename D::X_t;
public:
    distrscenariotree(const processdistribution<D>& p):fp(p)
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
    processdistribution<D> fp;
};

/// @} - processes


} // namespace


#endif // DISTRIBUTION_H
