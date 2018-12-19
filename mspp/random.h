#ifndef RANDOM_H
#define RANDOM_H

#include "mspp/commons.h"
#include <cmath>

namespace mspp
{

/// \addtogroup Distributions
/// @{

using probability = double;
const double probabilitytolerance = 1e-15;

template <typename I=double, typename C=novalue>
class distribution: public object
{
public:
    distribution() {}
    using I_t = I;
    using C_t = C;
};

template <typename F, typename S, typename C>
class jdistribution: virtual public distribution
        <pair<F, S>, C>
{
public:
    using F_t = F;
    using S_t = S;
};

template <typename X, typename C>
class vdistribution : virtual public distribution<vector<X>,C>
{
public:
    using X_t = X;
    vdistribution(unsigned int dim) : fdim(dim) {}
    unsigned int dim() const { return fdim; }
private:
    const unsigned int fdim;
};


/// \addtogroup mdisc MC distributions
/// \ingroup Distributions
/// @{

template <typename I, typename C>
class mdistribution: virtual public distribution<I, C>
{
public:
    I draw(const C& c) const
    {
        return do_draw(c);
    }
    I draw() const
    {
        static_assert(std::is_same<C,novalue>::value);
        return do_draw(novalue());
    }
private:
    virtual I do_draw(const C& c) const = 0;
};

/// @}

/// \addtogroup discretedisc Discrete distributions
/// \ingroup Distributions
/// @{

template <typename I>
struct atom
{
    using I_t = I;
    I x;
    probability p;
};

template <typename I, typename C>
class ddistribution: public mdistribution<I,C>
{
public:
    atom<I> operator () (unsigned int i, const C& c) const
    {
        return atom_is(i,c);
    }
    atom<I> operator () (unsigned int i) const
    {
        static_assert(std::is_same<C,novalue>::value);
        return atom_is(i,novalue());
    }
private:
    virtual atom<I> atom_is(unsigned int i, const C& c) const  = 0;
    virtual I do_draw(const C& c) const
    {
        double u = sys::uniform();
        double sum = 0;
        for(unsigned int i=0; i<numeric_limits<unsigned int>::max() ; i++)
        {
            atom<I> a = (*this)(i,c) ;
            sum += a.p;
            if(sum >= u)
                return a.x;
        }
        assert(0);
    }
};

/// /p listdefined=true indicates that atoms are primarily determined by
/// method /p atoms_are
template <typename I, typename C, bool listdef = false>
class fdistribution: virtual public ddistribution<I,C>
{
public:
    static constexpr bool flistdef = listdef;

    void atoms(vector<atom<I>>& a, const C& c) const
    {
        a.resize(0); // better safe than sorry
        atoms_are(a,c);
        probability p=0;
        #ifndef NDEBUG
            for(unsigned int i=0; i< a.size(); i++)
                p += a[i].p;
            assert(fabs(p-1.0) < probabilitytolerance);
        #endif
    }

    unsigned int natoms(const C& c) const
    {
        return natoms_is(c);
    }

private:
    virtual void atoms_are(vector<atom<I>>& a, const C& c) const
    {
        assert(!listdef);
        a.resize(0);
        unsigned int s= natoms_is(c);
        for(unsigned int i=0; i<s; i++)
            a.push_back((*this)(i,c));
    }
    virtual unsigned int natoms_is(const C& c) const
    {
        vector<atom<I>> a;
        atoms_are(a,c);
        assert(a.size());
        return a.size();
    }
    virtual atom<I> atom_is(unsigned int i, const C& c) const
    {
        assert(listdef);
        vector<atom<I>> a;
        atoms_are(a,c);
        assert(a.size());
        return a[i];
    }
};

template <typename I>
class ldistribution: public fdistribution<I,novalue, true>
{   
public:
    ldistribution(const vector<atom<I>>& atoms) :
        fatoms(atoms)
    {
    }

    ldistribution(const vector<I>& values) :
        fatoms(values.size())
    {
        assert(values.size());
        probability p=1.0/fatoms.size();
        for(unsigned int i=0; i<fatoms.size(); i++)
        {
            fatoms[i].x = values[i];
            fatoms[i].p = p;
        }
    }
private:
    virtual void atoms_are(vector<atom<I>>& a, const novalue&) const
    {
        a = fatoms;
    }

    virtual unsigned int natoms_is(const novalue&) const
    {
        return fatoms.size();
    }

    virtual atom<I> atom_is(unsigned int i, const novalue&) const
    {
        return fatoms[i];
    }
private:
    vector<atom<I>> fatoms;
};

template <typename X>
class lvdistribution: public ldistribution<vector<X>>,
        public vdistribution<X,novalue>
{
public:
    lvdistribution(const vector<atom<vector<X>>>& atoms) :
        ldistribution<vector<X>>(atoms),
        vdistribution<X,novalue>(atoms[0].x.size())
    {
        assert(atoms.size());
        for(unsigned int i=1; i<atoms.size(); i++)
            assert(atoms[i].x.size()==atoms[0].x.size());
    }

    lvdistribution(const vector<vector<X>>& values) :
        ldistribution<vector<X>>(values), vdistribution<X,novalue>(values[0].size())
    {
        assert(values.size());
        for(unsigned int i=1; i<values.size(); i++)
            assert(values[i].size()==values[0].size());
    }
};


template <typename I>
class diracdistribution: public ldistribution<I>
{
public:
    diracdistribution(const I& a)
        : ldistribution<I>({a}) {}
    I x() const
    {
        atom<I> a = (*this)(0);
        return a.x;
    }
};


template <typename X, bool listdefined>
lvdistribution<X> operator *(const fdistribution<X,novalue,listdefined>& x,
                             const fdistribution<X,novalue,listdefined>& y)
{
    assert(x.natoms(novalue()));
    assert(y.natoms(novalue()));
    vector<atom<vector<X>>> a(x.natoms(novalue())*y.natoms(novalue()));

    unsigned int dst = 0;
    for(unsigned int i=0; i<x.natoms(novalue()); i++)
        for(unsigned int j=0; j<y.natoms(novalue()); j++)
        {
            atom<X> xa = x(i,novalue());
            atom<X> ya = y(j,novalue());
            vector<X> r = { xa.x, ya.x };
            a[dst++]= { r, xa.p*ya.p};
        }
    return lvdistribution<X>(a);
}

/// @} - discrete distributions

/// \addtogroup iterativedists Iteratively defined distributions
/// \ingroup Distributions
/// @{


/// Iterative joint distribution, \p M is a mapping from the condition
/// of F into (condition of F, condition of G)
template <typename F, typename S, typename M>
class ijdistribution:
     public jdistribution<typename F::I_t, typename S::I_t, typename F::C_t>
{
public:
    using F_t = F;
    using S_t = S;
    using M_t = M;

    ijdistribution(const F& f, const S& s) : ff(f), fs(s)
    {
        static_assert(std::is_base_of<mapping<typename M::D_t,typename M::R_t>,M>::value);
        static_assert(std::is_same<typename M::D_t,
               pair<typename F::C_t,typename F::I_t>>::value);
        static_assert(std::is_convertible<typename M::R_t&,typename S::C_t&>::value);
    }
    const F& first() const { return ff; }
    const S& second() const { return fs; }
private:
    F ff;
    S fs;
};

// mapping for unconditional ijdistribution
template <typename D,typename R, typename M>
class uijmapping : public mapping<pair<nothing,D>,R>
{
public:
    virtual R operator() (const pair<nothing,D>& p) const
    {
        return M()(p.second);
    }
};

/// with no condition
template <typename F, typename S, typename M>
using uijdistribution = ijdistribution<F,S,uijmapping<typename F::I_t,typename S::C_t,M>>;

template <typename F, typename S, typename M>
class mijdistribution: public ijdistribution<F,S,M>,
                       public mdistribution<pair<typename F::I_t,typename S::I_t>,typename F::C_t>
{
public:
    using I_t = typename ijdistribution<F,S,M>::I_t;
    using C_t = typename ijdistribution<F,S,M>::C_t;
    mijdistribution(const F& d, const S& e) :
      ijdistribution<F,S,M>(d,e)
    {
    }
    virtual I_t do_draw(const C_t& c) const
    {
        typename F::I_t i = this->first().draw(c);
        pair<C_t,typename F::I_t> p = { c, i };
        typename S::I_t j = this->second().draw(M()(p));
        return {i,j};
    }
};

template <typename F, typename S, typename M>
class fijdistribution: public mijdistribution<F,S,M>,
        public fdistribution<pair<typename F::I_t,typename S::I_t>,
                                          typename F::C_t>
{
public:
    using I_t = typename ijdistribution<F,S,M>::I_t;
    using C_t = typename ijdistribution<F,S,M>::C_t;
    fijdistribution(const F& d, const S& e) :
      mijdistribution<F,S,M>(d,e)
    {
        static_assert(std::is_base_of<
            fdistribution<typename F::I_t, typename F::C_t>,F>::value);
        static_assert(std::is_base_of<
            fdistribution<typename S::I_t, typename S::C_t>,S>::value);
    }
private:
    virtual void atoms_are(vector<atom<I_t>>& a, const C_t& c) const
    {
        a.resize(0);
        vector<atom<typename F::I_t>> fa;
        this->first().atoms(fa,c);

        for(unsigned int i=0; i<fa.size(); i++)
        {
            vector<atom<typename S::I_t>> sa;
            this->second().atoms(sa,M()(pair<C_t,typename F::I_t>(c,fa[i].x)));
            for(unsigned int j=0; j<sa.size(); j++)
                a.push_back({ { fa[i].x, sa[j].x}, fa[i].p * sa[j].p});
        }
    }
};


template <typename F, typename S, typename M>
using umijdistribution = mijdistribution<F,S,uijmapping<typename F::I_t,typename S::C_t,M>>;


/// Iterativey defined vector distribution
/// Z is a mapping from a history into Z::R_t (the condition of D)
/// the first (at index 0) distrubiton may be conditioned, too
template <typename D,typename Z>
class ivdistribution:
    virtual public vdistribution<typename D::I_t,typename Z::R_t>
{
    static void constexpr check()
    {
        static_assert(std::is_base_of<zeta<typename Z::X_t, typename Z::R_t>,Z>::value);
        static_assert(std::is_same<typename Z::X_t,typename D::I_t>::value);
    }
public:
    using X_t = typename D::I_t;
    using Z_t = Z;
    ivdistribution(const vector<D>& d)
       : vdistribution<typename D::I_t,typename Z::R_t>(d.size()), fd(d)
    {
        check();
        assert(d.size());
    }
    ivdistribution(const D& d, unsigned int dim)
       : vdistribution<typename D::I_t,typename Z::R_t>(dim), fd(dim,d)
    {
        check();
        assert(dim);
    }
    void atoms(const scenario<X_t>& s, vector<atom<X_t>>& a) const
    {
        assert(s.size());
        assert(s.size() <= this->dim());
        fd[s.size()-1].atoms(a,Z()(s));
    }
    const D& d(unsigned int i) const
    {
        assert(i < fd.size());
        return fd[i];
    }
private:
    vector<D> fd;
};

/*
template <typename D, typename E, typename Z>
class dipdistribution:
        public ijdistribution<D,E,Z>,
        public ddistribution<pair<typename D::I_t,typename E::I_t>,
                                               typename D::C_t>
{
public:
    using I_t = typename iterativedistribution<D,E>::I_t;
    using C_t = typename iterativedistribution<D,E>::C_t;
    diterativedistribution(const D& d, const E& e) :
        iterativedistribution<D,E>(d,e)
    {
    }
private:
    virtual void atoms_are(vector<atom<I_t>>& a, const C_t& c) const
    {
        a.clear();
        vector<atom<typename D::I_t>> da;
        this->d().atoms(da,c);
        unsigned int nd=da.size();
        for(unsigned int i=0; i<nd; i++)
        {
            atom<typename D::I_t> x = da[i];

            vector<atom<typename E::I_t>> ea;
            this->e().atoms(ea,x.x);
            unsigned int ne=ea.size();

            for(unsigned int j=0; j<ne; j++)
            {
                atom<typename E::I_t> y = ea[j];
                a.push_back({{x.x,y.x},x.p*y.p});
            }
        }
    }

    virtual atom<I_t> atom_is(unsigned int i, const C_t& c) const
    {
        vector<atom<I_t>> a;
        atoms_are(a,c);
        assert(i<a.size());
        return a[i];
    }

    virtual unsigned int natoms_is(const C_t& c) const
    {
        vector<atom<I_t>> a;
        atoms_are(a,c);
        return a.size();
    }

    virtual I_t do_draw(const C_t& c) const
    {
        I_t r;
        r.i = this->d().draw(c);
        r.j = this->e().draw(r.i);
        return r;
    }
};

*/

/// @} - iterative distributions

/// @} - distribution


} // namespace


#endif // RANDOM_H

