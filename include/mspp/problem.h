#ifndef PROBLEM_H
#define PROBLEM_H

#include "mspp/commons.h"

namespace mspp
{

/// \addtogroup problems Problems
/// @{

class problemstructure: public std::vector<unsigned int>
{
public:
    problemstructure() {}
    problemstructure(const std::vector<unsigned int>& v ) :
       std::vector<unsigned int>(v) {}

    unsigned int sum(unsigned int k) const
    {
        assert(k < size());
        unsigned int s=0;
        for(int i=0; i<=k; i++)
            s+=operator[](i);
        return s;
    }

    unsigned int sum() const { return size() ? sum(size()-1) : 0; }

    unsigned int upto(unsigned int k) const
    {
        assert(k < size());
        return k==0 ? 0 : sum(k-1);
    }
};

using variable=double;

class vardef : public object
{
};

template <typename V>
using vardefs = std::vector<V>;



class cobase : public object
{
protected:
    cobase() { assert(0); } // should never be called, but has to be here to circumvent https://gcc.gnu.org/bugzilla/show_bug.cgi?id=58751
    cobase(const problemstructure& ps, unsigned int k) : fps(ps), fk(k) {}
    const problemstructure& ps() const { return fps; }
    const unsigned int k() const { return fk; }
private:private:
    void operator =(const cobase&); // to prevent assignment
    problemstructure fps;
    unsigned int fk;
};

using objective = cobase;

class constraint: public cobase
{
public:
    enum type {eq, geq, leq};
    constraint(const problemstructure& ps, unsigned int k):
        cobase(ps,k), ft(eq) {}

    type t() const { return ft; }
    void settype(type t) { ft = t; }
private:
    type ft;
};


template <typename G>
using ghconstraints = std::vector<ptr<G>>;

/// Base abstract class describing multistage problems
template<typename V, typename F, typename G, typename C>
class problem : public object
{
    public:
        problem(const problemstructure& ps)
            : d(ps)
        {
        }

        ///@{
        /// @name Accessors
        unsigned int T() const { return d.size()-1; }

        virtual void f(
                unsigned int k,
                const C& barxi,
                ptr<F> f) const
        {
            f= newf_is(k);
            f_is(k,barxi,*f);
        }


        /**
         * @brief
         * @param k stage
         * @param barxi history of the random vector up to #k
         * @param xs return value - variable ranges
         * @param g return value - constraints
         */
        void constraints(
                unsigned int k,
                const C& barxi,
                ptr<vardefs<V>>& r,
                ptr<ghconstraints<G>>& gh) const
        {
            r.reset(new vardefs<V>(d[k]));
            for(unsigned int i=0; i < d[k]; i++)
                (*r)[i] = newv_is(k,i);
            gh.reset(new ghconstraints<G>);
            this->constraints_are(k,barxi,*r,*gh);
        }

/*        std::string varname(unsigned int k) const
        {
            assert(k<d.sum());
            unsigned int s;
            for(s=0; s<T();s++)
            {
                if(d.upto(s+1) > k)
                    break;
            }
            std::ostringstream st;
            st << varname(s,k-d.upto(s)) << "_" << s;
            return st.str();
        }
*/
        virtual std::string varname(unsigned int stage, unsigned int i) const
        {
            return varname_is(stage,i);
        }

        const problemstructure d;
        ///@}

protected:
        G& addgh(ghconstraints<G>& g, unsigned int k) const
        {
            g.push_back(newg_is(k));
            return **g.rbegin();
        }

private:

        ///@{
        /// @name Interface towards descendants
        virtual ptr<G> newg_is(unsigned int k) const
        {
            return ptr<G>(new G(d,k));
        }

        virtual ptr<F> newf_is(unsigned int k) const
        {
            return ptr<F>(new F(d,k));
        }

        virtual ptr<V> newv_is(unsigned int k, unsigned int i) const
        {
            return ptr<V>(new V);
        }

        virtual void constraints_are(
                unsigned int k,
                const C& xi,
                vardefs<V>& xs,
                ghconstraints<G>& g
                ) const = 0;

        virtual void f_is(
                unsigned int k,
                const C& xi,
                F& f) const = 0;

        virtual std::string varname_is(unsigned int stage, unsigned int i) const
        {
            std::ostringstream s;
            s << "x" << stage << "_" << i;
            return s.str();
        }
};


/// \addtogroup pvar Variables
/// \ingroup problems
/// @{


class realvar : public vardef
{
public:
    enum type { R, Rplus, Rminus };

    void set(double l=minf, double h=inf)
    { fl=l; fh=h; }

    void set(type at)
    {
        switch(at)
        {
        case R: fl=minf; fh=inf; break;
        case Rplus: fl=0; fh=inf; break;
        case Rminus: fl=minf; fh=0; break;
        }
    }
    double l() const { return fl; }
    double h() const { return fh; }
private:
    void operator =(const realvar&); // to prevent assignment
    double fl;
    double fh;
};

///@}


/// \addtogroup Objectives
/// \ingroup problems
/// @{

class evobjective : virtual public objective
{
public:
    double operator()(const std::vector<variable> x) const
    {
        return objvalue_is(x);
    }
private:
    virtual double objvalue_is(const std::vector<variable>& x) const = 0;
};

class convexobjective: virtual public objective
{
};

class sgcobjective: public convexobjective, public evobjective
{
public:
    std::vector<double> sg(const std::vector<variable> x) const
    {
        return sg_is(x);
    }

private:
   virtual std::vector<double> sg_is(const std::vector<variable>& x) const = 0;
};


class linearobjective: public sgcobjective
{
public:
    linearobjective(const problemstructure& ps, unsigned int k):
        objective(ps,k), fc(ps[k-1]) {}

    void set(const std::vector<double>& c)
    {
        assert(c.size()==fc.size());
        fc = c;
    }
    double& operator[](unsigned int i)
    {
        assert(i < fc.size());
        return fc[i];
    }

private:
    virtual std::vector<double> sg_is(const std::vector<variable>& x) const
    {
        assert(x.size()==fc.size());
        return fc;
    }
    virtual double objvalue_is(const std::vector<variable>& x) const
    {
        assert(x.size()==fc.size());
        double s=0;
        for(unsigned int i=0; i<fc.size(); i++)
            s += fc[i] * x[i];
        return s;
    }
private:
    std::vector<double> fc;
};


///@}


/// \addtogroup Constraints
/// \ingroup problems
/// @{

class linearconstraint : virtual public constraint
{
public:

    linearconstraint( const problemstructure& ps, unsigned int k )
        : flhssize(ps.sum(k)), frhs(0) {}

    virtual double lhs(unsigned int k) const = 0;
    unsigned int lhssize() const { return flhssize; }

    double rhs() const { return frhs; }
    void setrhs(double rhs ) { frhs = rhs; }
private:
    unsigned int flhssize;
    double frhs;
};

class fulllinearconstraint : public linearconstraint
{
public:
    fulllinearconstraint& operator=(fulllinearconstraint& c) { return c; }
       // added because gcc deletes its default version


    fulllinearconstraint( const problemstructure& ps, unsigned int k )
        : constraint(ps,k), linearconstraint(ps,k),
           flhs(lhssize(),0) {}

    double lhs(unsigned int k) const { assert(k<flhs.size()); return flhs[k]; }
    void setlhs(const std::vector<double>& lhs)
    {
        assert(flhs.size() == lhs.size());
        flhs = lhs;
    }
    void set(const std::vector<double>& lhs,type t,double rhs)
    {
        setlhs(lhs);
        settype(t);
        setrhs(rhs);
    }

private:
    std::vector<double> flhs;
};


class interstagelinearconstraint : public linearconstraint
{
public:
    interstagelinearconstraint& operator=(interstagelinearconstraint& c) { return c; }
       // added because gcc deletes its default version

    interstagelinearconstraint( const problemstructure& ps, unsigned int k )
        : constraint(ps,k), linearconstraint(ps,k),
          flast(k ? ps[k-1] : 0,0), fcurrent(ps[k],0) {}

    const std::vector<double>& last() const  { return flast; }
    double last(unsigned int k) const { assert(k<flast.size()); return flast[k]; }
    unsigned int lastsize() const { return flast.size(); }

    const std::vector<double>& current() const  { return fcurrent; }
    double current(unsigned int k) const { assert(k<fcurrent.size()); return fcurrent[k]; }
    unsigned int currentsize() const { return fcurrent.size(); }

    double lhs(unsigned int k) const
    {
        assert(k<lhssize());
        unsigned int cindex = lhssize()-fcurrent.size();
        unsigned int lindex = cindex-fcurrent.size();
        if(k<lindex)
            return 0;
        else if(k>=cindex)
            return fcurrent[k-cindex];
        else
            return flast[k-lindex];
    }

    void setlast(const std::vector<double>& last)
    {
        assert(flast.size() == last.size());
        flast = last;
    }
    void setcurrent(const std::vector<double>& current)
    {
        assert(fcurrent.size() == current.size());
        fcurrent = current;
    }
    void set(const std::vector<double>& last,
             const std::vector<double>& current, type t,double rhs)
    {
        setlast(last);
        setcurrent(current);
        settype(t);
        setrhs(rhs);
    }

private:
    std::vector<double> fcurrent;
    std::vector<double> flast;
};


///@}


///@}

} // namespace

#endif // PROBLEM_H
