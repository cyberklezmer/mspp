#ifndef COMMONS_H
#define COMMONS_H

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>
#include <assert.h>
#include <limits>
#include <random>

static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required");

static_assert(__cplusplus >= 201703, "C++17 language standard required");

// static_assert(RAND_MAX >=2147483647, "RAND_MAX too small");


using namespace std;

namespace mspp
{

/** \addtogroup sysclasses System Classes
 *  @{
 */

/// Base class
class object
{
public:
    virtual ~object() {}
};

/// \ingroup Distributions
using probability = double;


/// Exception
class exception : public std::exception, public object
{
public:
    exception(const std::string &error = "unknown mspp error", unsigned int erno = 0)
        : fmsg(error), ferno(erno)
    {
    }

    const std::string& msg() const  { return fmsg; }
    unsigned int erno() const { return ferno; }
private:
    std::string fmsg;
    unsigned int ferno;
};

/// \brief System class, implementing output and random generation.
///
/// Has only static members, needs not to be instatiated.
class sys : public object
{
    static sys& self()
    {
       static sys s;
       return s;
    }
public:
    sys() : fout(0), flog(0), ferr(0), funiform(0.0,1.0)
    {
    }
    ~sys()
    {
    }

    /// Sets an existing stream as the text output of the library.
    static void setout(std::ostream& o)
    {
        self().fout = &o;
    }

    /// Resets the standard ouptut of the library,
    /// redirecting it back to \p std::cout
    static void resetout()
    {
        self().fout = 0;
    }

    static void setlog(std::ostream& l)
    {
        self().flog = &l;
    }
    static void resetlog()
    {
        self().flog=0;
    }

    static void seterr(std::ostream& e)
    {
        self().ferr = &e;
    }
    static void reseterr()
    {
        self().ferr = 0;
    }

    static std::ostream& out()
    {
        return self().fout ? *(self().fout) : std::cout;
    }
    static std::ostream& log()
    {
        return self().flog ? *(self().flog) : std::clog;
    }
    static std::ostream& err()
    {
        return self().ferr ? *(self().ferr) : std::cerr;
    }

    /// Resets the seed of the random generator according to computer time.
    static void seed()
    {
        self().fengine.seed(time(0));
    }

    /// Resets the seed of the random generator.
    static void seed(unsigned int aseed)
    {
        self().fengine.seed(aseed);
    }

    /// Returns a pesudo-random observaton from the uniform distribution on [0,1]
    static probability uniform()
    {
        return self().funiform(self().fengine);
    }

private:
    std::ostream* fout;
    std::ostream* flog;
    std::ostream* ferr;
    std::default_random_engine fengine;
    std::uniform_real_distribution<double> funiform;

};

/// Heap pointer to \p T
template <typename T>
using ptr=std::shared_ptr<T>;


/// @}


/// \addtogroup general General Definitions
///  @{


/// Value reperesenting infinity with respect to \p V
template <typename V>
constexpr V inf() { return std::numeric_limits<V>::infinity(); }

/// Value reperesenting minus infinity with respect to \p V
template <typename V>
constexpr V minf() { return -inf<V>(); }

/// Maximal value of \p V
template <typename V>
constexpr V max() { return std::numeric_limits<V>::max(); }

/// Minimal value of \p V
template <typename V>
constexpr V min() { return std::numeric_limits<V>::lowest(); }


//template <typename T>
//using vector=std::vector<T>;

/// Vector of vectors
template <typename T>
class vectors : public vector<vector<T>>
{
public:
    vectors() {}
    vectors(unsigned int dim) : vector<vector<T>>(dim) {}
    vectors(const vector<vector<T>>& v): vector<vector<T>>(v) {}

    void converttoonedim(vector<T>& d)
    {
        d.clear();
        for(typename vectors<T>::iterator i = this->begin();
             i != this->end();
             i++)
            for(typename vector<T>::iterator j = i->begin();
                 j != i->end();
                 j++)
                d.push_back(*j);
    }

    operator vector<T>() const
    {
        vector<T> d;
        converttoonedim(d);
        return d;
    }

    unsigned int totaldim() const
    {
        unsigned int r=0;
        for(typename vectors<T>::iterator i = this->begin();
             i != this->end();
             i++)
            r += (*this)[i].size();
        return r;
    }
};

/// \addtogroup pvar Variables
/// @{

/// Real decision variable
using realvar=double;
/// Integer decision variable
using intvar=long;
/// Binary decision variable
using binvar=bool;
/// Mixed decision variable
union mixedvar
{
    realvar x;
    intvar n;
    binvar b;
    mixedvar(const realvar& ax) : x(ax) {}
    mixedvar(const intvar& an) : n(an) {}
    mixedvar(const binvar& ab) : b(ab) {}
    operator realvar() const { return x; }
    operator intvar() const { return n; }
    operator binvar() const { return b; }
//    mixedvar& operator =(realvar& v) { x=v; return *this; }
//    mixedvar& operator =(intvar& v) { n=v; return *this; }
};


/// Vector of variables of type \p V
template <typename V>
using variables=vector<V>;


/// Range (interval of possible values) of variables \p V
template <typename V>
class range : public object
{
public:
    using V_t = V;

    enum type { realt, intt, bint };

    range()
    {
        if constexpr( std::is_same<V,realvar>::value)
           ft = realt;
        else if constexpr( std::is_same<V,intvar>::value)
           ft = intt;
        else if constexpr( std::is_same<V,binvar>::value)
           ft = bint;
        else
        {
           assert(0);
        }

        setlimits();
    }


    range(type t, V l=minf<V>(), V h=inf<V>())
    {
        if constexpr( std::is_same<V,realvar>::value)
        {
           assert(t==realt);
           ft = realt;
        }
        else if constexpr( std::is_same<V,intvar>::value)
        {
            assert(t==intt);
            ft = intt;
        }
        else if constexpr( std::is_same<V,binvar>::value)
        {
            assert(t==bint);
            ft = bint;
        }
        else if constexpr( std::is_same<V,mixedvar>::value)
        {
           assert(t==realt || t==intt || t==bint);
           ft = t;
        }
        else
            assert(0);
        setlimits(l,h);
    }

/*    void setreal(realvar l=minf<realvar>(), realvar h=inf<realvar>())
    {
        static_assert(std::is_same<V,mixedvar>::value);
        ft = realt;
        setlimits(l,h);
    }

    void setint(intvar l=minf<intvar>(), intvar h=inf<intvar>())
    {
        static_assert(std::is_same<V,mixedvar>::value);
        ft = intt;
        setlimits(l,h);
    }
*/
    void setlimits(V l=minf<V>(), V h=inf<V>())
    {
        assert(ft==realt || ft==intt);
        fl=l; fh=h;
    }

    void setpositive()
    {
       setlimits(0);
    }
    V l() const
    {
        assert(ft==realt || ft==intt);
        return fl;
    }
    V h() const
    {
        assert(ft==realt || ft==intt);
        return fh;
    }
    type t() const
    {
        assert(ft>= realt && ft <= bint);
        return ft;
    }
    bool ishinf() const
    {
        assert(ft==realt || ft==intt);
        if(ft==realt)
        {
            realvar h = fh;
            return h==inf<realvar>();
        }
        else
        {
            intvar h = fh;
            return h==inf<int>();
        }
    }
    bool islinf() const
    {
        assert(ft==realt || ft==intt);
        if(ft==realt)
        {
            realvar l = fl;
            return l==-inf<realvar>();
        }
        else
        {
            intvar l = fl;
            return l==-inf<int>();
        }
    }
private:
    V fl;
    V fh;
    type ft;
};


///@}



/// \addtogroup cns Vector restrictions
/// @{

/// Base class for restrictions of decision variables history
///
///
class tilde: public object
{
public:
    bool included(unsigned int i, unsigned int s) const
    {
        assert(i<s);
        return is_included(i,s);
    }
    bool included(unsigned int i, unsigned int j, unsigned int s) const
    {
        assert(i<s);
        return is_included(i,j,s);
    }

private:    
    virtual bool is_included(unsigned int i, unsigned int s) const = 0;
    virtual bool is_included(unsigned int i, unsigned int j, unsigned int s) const = 0;
};

/// Not restricting variable restriction
class allx: public tilde
{
private:
    virtual bool is_included(unsigned int, unsigned int ) const { return true; }
    virtual bool is_included(unsigned int, unsigned int, unsigned int) const { return true; }
};

/// Restriction of variables to those from the last period
class lastx: public tilde
{
private:
    virtual bool is_included(unsigned int i, unsigned int s) const
       { return i==s-1; }
    virtual bool is_included(unsigned int i, unsigned int j, unsigned int s) const
       { return i==s-1; }
};


/// @}


///@}


/// \addtogroup problems Problems
/// @{

/// Base class of decision problems' constraints
class constraint: public object
{
public:
    enum type {eq, geq, leq};
    constraint(unsigned int xdim, type t=eq)
        : fxdim(xdim), ft(t) {}

    type t() const { return ft; }
    void settype(type t) { ft = t; }
    unsigned int xdim() const { return fxdim; }
private:
    type ft;
    unsigned int fxdim;
};

const constraint::type eq = constraint::eq;
const constraint::type leq = constraint::leq;
const constraint::type geq = constraint::geq;

/// Base class for decision criterions
class criterion: public object
{
};

///@}


/// \addtogroup fns Functions
/// @{

/// A class bearing no information (alternative to \p void)
struct nothing {};

constexpr nothing na;

/// \brief Base class for mappings
/// \tparam D domain
/// \tparam R range
template <typename D,typename R>
class mapping: public object
{
public:
    using D_t=D;
    using R_t=R;
    virtual R operator() (const D& ) const = 0;
};

/// Function (mapping from \p D to real numbers)
template <typename D>
using function = mapping<D,double>;

/// Identity mapping
template <typename I>
class idmapping: public mapping<I,I>
{
public:
    virtual I operator() (const I& i) const { return i; }
};

/// \brief Mapping having no value
///
/// Used as a condition of unconditional distributions, for instance.
template <typename I>
class nomapping: public mapping<I,nothing>
{
public:
    virtual nothing operator() (const I& i) const { return na; }
};

/// Function from Eucleidean space
class efunction: public function<vector<double>>
{
public:
    efunction(unsigned int xdim) : fxdim(xdim) {}
    unsigned int xdim() const { return fxdim; }
    double operator() (const vector<double>& x) const
    {
        assert(x.size() == fxdim);
        return value_is(x);
    }
private:
    unsigned int fxdim;
    virtual double value_is(const vector<double>& x) const = 0;
};

/// Convex function from an Eucleidean space
class convexefunction: public efunction
{
public:
    convexefunction(unsigned int xdim) : efunction(xdim) {}
};

/// Function with defined subdifferential
class subdifefunction: public convexefunction
{
public:
    subdifefunction(unsigned int dim) : convexefunction(dim) {}

    vector<double> sg(const vector<double>& x) const
    {
        return sg_is(x);
    }
private:
   virtual vector<double> sg_is(const vector<double>& x) const = 0;
};

/// Linear function
class linearfunction: public subdifefunction
{
public:
    linearfunction(unsigned int dim):
        subdifefunction(dim), fc(dim) {}

    linearfunction(vector<double> c):
        subdifefunction(c.size()), fc(c) {}

    void setc(const vector<double>& c)
    {
        assert(c.size()==xdim());
        fc = c;
    }
    void setc(unsigned int i, double c)
    {
        assert(i<fc.size());
        fc[i] = c;
    }
    double c(unsigned int i) const
    {
        assert(i < fc.size());
        return fc[i];
    }
    const vector<double>& c() const
    {
        return fc;
    }

private:
    virtual vector<double> sg_is(const vector<double>& x) const
    {
        assert(x.size()==xdim());
        return fc;
    }
    virtual double value_is(const vector<double>& x) const
    {
        assert(x.size()==xdim());
        double s=0;
        for(unsigned int i=0; i<fc.size(); i++)
            s += fc[i] * x[i];
        return s;
    }
private:
    vector<double> fc;
};

///@}



/// \ingroup Processes

/// Scenario (history) of a stochastic process with values in \p X
template <typename X=double>
using scenario = vector<X>;

/// \ingroup Processes
/// \brief Base class for scenario restrictions
/// \tparam X type of the scenario components
/// \tparam R type into which the scenario is transformed
///
/// Scenario restrictions are mappings transforming scenarios
/// into some other space, which server either as conditions
/// of stochastic processes or as random parameters of decision problems.
///
template <typename X,typename R>
class zeta: public mapping<scenario<X>,R>
{
public:
    using R_t=R;
    using X_t=X;
};

/// \ingroup Processes

/// \brief Not restricting scenario restriction
/// \tparam X scenario component type
///
/// Transforms the scenario into itself.
template <typename X>
class allxi: public zeta<X, scenario<X>>
{
public:
    virtual scenario<X> operator() (const scenario<X>& s) const
       { return s; }
};

/// \ingroup Processes
/// \brief Restriction of a scenario into a function of its last value
/// \tparam X scenario component type
/// \tparam M mapping the last value is transformed with
/// (identity by default)
///
/// The range of \p lastxi is the same that the one of \p M.
template <typename X, typename M=idmapping<X>>
class lastxi: public zeta<X,typename M::R_t>
{
public:
    using M_t = M;
    lastxi() {}
    lastxi(const M& m) : fm(m) {}
    const M& m() const { return fm; }
    virtual typename M::R_t operator() (const scenario<X>& s) const
       { assert(s.size()); return fm(s[s.size()-1]);}
private:
   M fm;
};

/// \ingroup Processes

/// \brief Muting (complete restriction) of a scenario.
/// \tparam X scenario component type
template <typename X>
class noxi: public zeta<X,nothing>
{
public:
    virtual nothing operator() (const scenario<X>&) const
       { return na; }
};

} // namespace


/// \defgroup sms Solution Methods


#endif // COMMONS_H
