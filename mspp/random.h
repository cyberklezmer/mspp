#ifndef RANDOM_H
#define RANDOM_H

#include "mspp/commons.h"
#include <cmath>

namespace mspp
{

/// \addtogroup Distributions
/// @{

const double probabilitytolerance = 1e-15;

/// \brief Base class.
///
///  \tparam I type of the values
///  \tparam C type of the condition
template <typename I=double, typename C=novalue>
class distribution: public object
{
public:
    distribution() {}
    using I_t = I;
    using C_t = C;
};

/// \brief Joint distribution
/// \tparam F type of the first coordinate
/// \tparam S type of the second coordinate
/// \tparam C type of the condition
template <typename F, typename S, typename C>
class jdistribution: virtual public distribution
        <pair<F, S>, C>
{
public:
    using F_t = F;
    using S_t = S;
};


/// \brief Vector distribution
/// \tparam X type of the vector coordinates
/// \tparam C type of the condition
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

/// \brief Monte carlo distribution
/// \tparam I type of the values
/// \tparam C type of the condition

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

/// \addtogroup realdists Real distributions
/// \ingroup Distributions
/// @{

template <typename C>
using realdistribution=distribution<double, C>;

template <typename C>
class cdfdistribution: virtual public realdistribution<C>
{
public:
    probability cdf(double x) const
    {
        return this->cdf(x,novalue());
    }
    probability cdf(double x, const C& c) const
    {
        return cdf_is(x,c);
    }
private:
    virtual probability cdf_is(double, const C& c) const = 0;
};


template <typename C>
class qdistribution: virtual public cdfdistribution<C>,
        public mdistribution<double,C>
{
public:
    double quantile(probability p) const
    {
        return this->quantile(p,novalue());
    }
    double quantile(probability p, const C& c) const
    {
        return quantile_is(p,c);
    }
private:
    virtual double quantile_is(probability, const C&) const = 0;
       // a numeric method may be implemnted here
    virtual double do_draw(const C& c) const
    {
        probability u = sys::uniform();
        return quantile(u,c);
    }
};

template <typename D>
class scaleddistribution : public qdistribution<novalue>
{
    void check()
    {
        assert(fsd > 0);
        static_assert(std::is_base_of<qdistribution<novalue>,D>::value);
    }

public:
    scaleddistribution(double m, double sd) : fm(m), fsd(sd)
    {
        check();
    }

    scaleddistribution(const D& d, double m, double sd)
        : fm(m), fsd(sd), fd(d)
    {
        check();
    }
    /// returns the distribution which whas scalled
    const D& srcd() const
    {
        return fd;
    }
    const double m() const { return fm; }
    const double sd() const { return fsd; }
private:
    virtual double quantile_is(probability p, const novalue&) const
    {
        return fm + fsd*fd.quantile(p);
    }
    virtual probability cdf_is(double x, const novalue&) const
    {
        return fd.cdf((x-fm) / fsd);
    }

    double fm;
    double fsd;
    D fd;
};



template <typename D>
class ardistribution : public qdistribution<double>
{
    void check()
    {
        static_assert(std::is_base_of<qdistribution<novalue>,D>::value);
    }

public:
    ardistribution(double a=1) : fa(a)
    {
        check();
    }

    ardistribution(const D& d, double a=1) : fa(a), fd(d)
    {
        check();
    }
    double a() const { return fa; }
    const D& srcd() const { return fd; }
private:
    virtual double quantile_is(probability p, const double& c) const
    {
        return fa * c + fd.quantile(p);
    }
    virtual probability cdf_is(double x, const double& c) const
    {
        return fd.cdf(x-fa * c);
    }

    double fa;
    D fd;
};

/// @} // real distributions


/// \addtogroup discretedisc Discrete distributions
/// \ingroup Distributions
/// @{



/// \brief Atom of #ddistribution
/// \tparam I type of the values
template <typename I>
struct atom
{
    using I_t = I;
    I x;
    probability p;
};

/// \brief Discrete distribution
/// \tparam I type of the values
/// \tparam C type of the condition
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
        return this->drawuniv(c);
    }
protected:
    /// generates an uniform variable \p u and goes through the probabilities
    /// until their sum exceeds \p u. May be very slow
    I drawuniv(const C& c) const
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

/// \brief Finite discrete distribution
/// \tparam I type of the values
/// \tparam C type of the condition
/// \tparam listdef vay of atoms definition:
///
/// \p listdef=true indicates
/// that the atoms are primarily determined by an overloading of
/// \p atoms_are while \p atom_is and \p natom_is need not be overoaded
/// as they use \p atoms_are. If
/// \p listdef=false then the atoms are
/// determided by overloading of \p natom_is and \p atom_is while \p atoms_are needs not
/// be overloaded.
template <typename I, typename C, bool listdef = false>
class fdistribution: virtual public ddistribution<I,C>
{
public:
    static bool constexpr flistdef=listdef;
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

    void atoms(vector<atom<I>>& a) const
    {
        static_assert(std::is_same<C,novalue>::value);
        return atoms(a,novalue());
    }

    unsigned int natoms(const C& c) const
    {
        return natoms_is(c);
    }

    unsigned int natoms() const
    {
        static_assert(std::is_same<C,novalue>::value);
        return natoms_is(novalue());
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
    /// returning \p false means that slow \ref drawuniv will be used in MC.
    virtual bool is_equiprobable(const C& c) const { return false; }
    virtual I do_draw(const C& c) const
    {
        if(is_equiprobable(c))
        {
            double u = sys::uniform();
            unsigned int n = natoms(c);
            unsigned int i = static_cast<unsigned int>(floor(n*u));
            return (*this)(i,c).x;
        }
        else
            return this->drawuniv(c);
    }
};


/// \brief List defined distribution
/// \tparam I type of the values
///
///
template <typename I>
class ldistribution: public fdistribution<I,novalue, true>
{   
public:
    ldistribution(const vector<atom<I>>& atoms) :
        fatoms(atoms), fequiprobable(false)
    {
        assert(fatoms.size());
        probability tp = fatoms[0];
        for(unsigned int i=1; i<fatoms.size(); i++)
        {
            if(tp != fatoms[i].p)
                return;
        }
        fequiprobable = true;
    }

    ldistribution(const vector<I>& values) :
        fatoms(values.size()), fequiprobable(true)
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
    bool fequiprobable;
    virtual bool is_equiprobable(const novalue&) const { return fequiprobable; }
};

/// \brief Alternative list defined distribution
/// \tparam I type of the values
///
/// The atoms sets are based on the condition (of type <tt>unsigned int</tt>).
template <typename I>
class altldistribution: public fdistribution<I,unsigned int, true>
{
    vector<ptr<ldistribution<I>>> makedists(const vector<vector<I>>& values)
    {
        vector<ptr<ldistribution<I>>> ret(values.size());
        for(unsigned int i=0; i< values.size(); i++)
            ret[i].reset(new ldistribution<I>(values[i]));
        return ret;
    }
public:
/*    altldistribution(const vector<vector<atom<I>>>& atoms) :
        fdists(makedists(atoms))
    {
    }*/

private:
    virtual void atoms_are(vector<atom<I>>& a, const unsigned int& c) const
    {
        assert(c<fdists.size());
        fdists[c]->atoms(a,novalue());
    }

    virtual unsigned int natoms_is(const unsigned int& c) const
    {
        assert(c<fdists.size());
        return fdists[c]->natoms(novalue());
    }

    virtual atom<I> atom_is(unsigned int i,  unsigned int& c) const
    {
        assert(c<fdists.size());
         return fdists[c]->atom(i,novalue());
    }

    virtual bool is_equiprobable(const unsigned int& c) const
    { return fdists[c]->is_equiprobable(); }

    vector<ptr<ldistribution<I>>> fdists;

};


/// \brief List defined vector distribution
/// \tparam X type of the values' components
///
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


/// \brief Dirac distribution
/// \tparam I type of the value
///
template <typename I>
class diracdistribution: public ldistribution<I>
{
    std::vector<I> point(const I& a) const
    {
       std::vector<I> ret = {a};
       return ret;
    }
public:
    diracdistribution(const I& a)
        : ldistribution<I>(point(a)) {}
    I x() const
    {
        atom<I> a = (*this)(0);
        return a.x;
    }
};

/// \brief Product distribution
template <typename X, bool listdefined>
lvdistribution<X> operator *(const fdistribution<X,novalue,listdefined>& x,
                             const fdistribution<X,novalue,listdefined>& y)
{
    assert(x.natoms());
    assert(y.natoms());
    vector<atom<vector<X>>> a(x.natoms()*y.natoms());

    unsigned int dst = 0;
    for(unsigned int i=0; i<x.natoms(); i++)
        for(unsigned int j=0; j<y.natoms(); j++)
        {
            atom<X> xa = x(i);
            atom<X> ya = y(j);
            vector<X> r = { xa.x, ya.x };
            a[dst++]= { r, xa.p*ya.p};
        }
    return lvdistribution<X>(a);
}

/// @} - discrete distributions

/// \addtogroup iterativedists Iteratively defined distributions
/// \ingroup Distributions
/// @{


/// \brief Iterative joint distribution,
/// \tparam F the first distribution
/// \tparam S the second distribution
/// \tparam M a mapping from the
/// <tt>pair<condition of F, condition of S></tt> into
/// <tt>condition of \p F</tt>
///
/// Tne condtition of the distribution is the same type as that of \p F.
/// The distribution itself is defined iteratively by <em>F</em>|(condition of \e F) and
/// <em>S</em>|<em>M</em>(value of \e F,condition of <em>F</em>)

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

/// Helper class, used as \p M in \p uijdistribution
template <typename D,typename R, typename M>
class uijmapping : public mapping<pair<novalue,D>,R>
{
public:
    virtual R operator() (const pair<nothing,D>& p) const
    {
        return M()(p.second);
    }
};

/// \brief Unconditional iterative distribution
/// \tparam F the first distribution
/// \tparam S the second distribution
/// \tparam M s mspping from the value of \p F into the condition of \f S
///
/// By default, \p M is the identity from value of \p F into condition of \p S.
/// If those types differ, then \p M should be explicitly specified.

template <typename F, typename S, typename M=idmapping<typename F::I_t>>
using uijdistribution = ijdistribution<F,S,uijmapping<typename F::I_t,typename S::C_t,M>>;

/// \brief MC iterative distribution
/// \tparam F the first distribution
/// \tparam S the second distribution
/// \tparam M a mapping from the condition of \p F into
/// <tt>(condition of F, condition of S)</tt>
///
/// See \ref ijdistribution
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

/// \brief Monte Carlo iterative distribution with finite support
/// \tparam F the first distribution
/// \tparam S the second distribution
/// \tparam M a mapping from the condition of \p F into
/// <tt>(condition of F, condition of S)</tt>
///
/// See \ref ijdistribution
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

/// \brief Unconditional MS iterative distribution
/// \tparam F the first distribution
/// \tparam S the second distribution
/// \tparam M s mspping from the value of \p F into the condition of \f S
///
/// See \ref uijdistribution

template <typename F, typename S, typename M>
using umijdistribution = mijdistribution<F,S,uijmapping<typename F::I_t,typename S::C_t,M>>;


/// \brief Iteratively defined vector conditional vector distribution
/// \tparam D conditional distribution of the vector's component
/// \tparam Z mapping from a history of the distribution into the condition of \p D
///
/// The distribution is defned iteratively, the \e i-th component being
/// conditioned by <tt>{Z}(s)</tt>$ where \p s is the scenario
/// containing previous values of the vector.
/// Template parameter \p Z has to be a descendant of \ref zeta<\p X, \p C>
/// where \p X is the type of \p D's values and \p C is the type
/// of its condition.
///

template <typename E, typename D, typename Z>
class ivdistribution:
    virtual public vdistribution<typename D::I_t,novalue>
{
    static void constexpr check()
    {
        static_assert(std::is_base_of<zeta<typename Z::X_t, typename Z::R_t>,Z>::value);
        static_assert(std::is_base_of<
                      distribution<typename D::I_t,typename D::C_t>,D>::value);
        static_assert(std::is_base_of<
                   distribution<typename E::I_t,novalue>,E>::value);
        static_assert(std::is_same<typename Z::R_t,typename D::C_t>::value);
        static_assert(std::is_same<typename Z::X_t,typename D::I_t>::value);
        static_assert(std::is_same<typename E::I_t,typename D::I_t>::value);
    }
public:
    using X_t = typename D::I_t;
    using E_t = E;
    using Z_t = Z;
//    ivdistribution(const E& e, const vector<D>& d)
//       : vdistribution<typename D::I_t,novalue>(d.size()+1), fd(d), fe(e)
//    {
//       check();
//    }
    ivdistribution(const E& e, const D& d, unsigned int dim)
       : vdistribution<typename D::I_t,novalue>(dim), fe(e), fd(dim-1,d)
    {
        check();
        assert(dim);
    }

    ivdistribution(const X_t& x0, const D& d, unsigned int dim)
       : vdistribution<typename D::I_t,novalue>(dim), fe(diracdistribution<X_t>(x0)), fd(dim-1,d)
    {
        check();
        assert(dim);
    }

    vector<X_t> draw() const
    {
        vector<X_t> ret;
        X_t x = fe.draw(novalue());
        ret.push_back(x);
        typename Z::R_t c = Z(ret);
        for(unsigned int i=0; i<this->dim()-1; i++)
        {
            x = this->fd[i].draw(c);
            ret.push_back(x);
            c = Z(ret);
        }
        return ret;
    }

    /// may be used only if \p D is a descendant of \ref fdistribution
    void atoms(const vector<X_t>& s, vector<atom<X_t>>& a) const
    {
        if(!s.size())
            fe.atoms(a);
        else
        {
            assert(s.size() < this->dim());
            fd[s.size()-1].atoms(a,Z(s));
        }
    }
    const D& d(unsigned int i) const
    {
        assert(i>0);
        assert(i-1 < fd.size());
        return fd[i-1];
    }
    const E& e() const
    {
        return fe;
    }
private:
    E fe;
    vector<D> fd;
};


template <typename E, typename D, typename Z, typename M>
class mivdistribution:
    public ivdistribution<E,D,Z>
{
    static void constexpr check()
    {
        static_assert( std::is_base_of<distribution<typename M::I_t,novalue>, M>::value);
        static_assert( std::is_same<typename M::I_t,typename D::I_t>::value);
    }
public:
    mivdistribution(const E& e, const D& d, unsigned int dim)
       : ivdistribution<E,D,Z>(e,d,dim),
         vdistribution<typename D::I_t,novalue>(dim)
    {
        check();
        assert(dim);
    }

    mivdistribution(const typename E::I_t& x0, const D& d, unsigned int dim)
        : ivdistribution<E,D,Z>(x0,d,dim),
          vdistribution<typename D::I_t,novalue>(dim)
    {
        check();
        assert(dim);
    }
    M md(unsigned int i) const
    {
        assert(i>0);
        assert(i<this->dim());
        return md_is(i);
    }
private:
    virtual M md_is(unsigned int i) const = 0;
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

