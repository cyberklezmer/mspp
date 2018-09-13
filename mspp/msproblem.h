#ifndef MSPROBLEM_H
#define MSPROBLEM_H

#include "mspp/linear.h"

namespace mspp
{


/// \addtogroup msproblems Multistage Problems
/// \ingroup problems
/// @{


class msproblemstructure: public vector<unsigned int>
{
public:
    msproblemstructure() {}
    msproblemstructure(unsigned int dim):vector<unsigned int>(dim) {}
    msproblemstructure(const vector<unsigned int>& v ) :
       vector<unsigned int>(v) {}

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


template <typename G>
class msconstraints
{
public:
    msconstraints(unsigned int n) : fn(n) {}
    G& add(const G& g)
    {
        assert(g.xdim()==fn);
        fgs.push_back(g);
        return *fgs.rbegin();
    }
    G& add()
    {
        fgs.push_back(G(fn));
        return *fgs.rbegin();
    }
    G& operator [] (unsigned int i)
    {
        assert(i<fgs.size());
        return fgs[i];
    }

    operator vector<G> () const
    {
        return fgs;
    }


private:
    unsigned int fn;
    vector<G> fgs;
};


template <typename V>
class ranges
{
public:
    ranges(unsigned int n) : frs(n)
    {
    }
    range<V>& operator [] (unsigned int i)
    {
        assert(i<frs.size());
        return frs[i];
    }
    operator vector<range<V>> () const
    {
        return frs;
    }

private:
    vector<range<V>> frs;
};


/// Base abstract class describing multistage problems
template<typename C, typename F, typename G,
          typename V, typename I=double,
          typename R=allx, typename Z=allxi<I>>
class msproblem : public object
{
    public:

        using C_t = C;
        using F_t = F;
        using G_t = G;
        using V_t = V;
        using I_t = I;
        using Z_t = Z;
        using R_t = R;
        using O_t = typename Z::C_t;


        msproblem(const msproblemstructure& ps, C arho=C()):
             d(ps), rho(arho)
        {
        }

        msproblem(const vector<unsigned int>& ps, C arho=C()):
             d(ps), rho(arho)
        {
        }

        ///@{
        /// @name Accessors
        unsigned int T() const { return d.size()-1; }

        virtual F f(const scenario<I>& s) const
        {
            assert(s.size());
            return f(s.size()-1,zeta(s));
        }

        virtual F f(
                unsigned int k,
                const O_t& zeta) const
        {
            F f(this->xdim(k));
            f_is(k,zeta,f);
            assert(f.xdim()==xdim(k));
            return f;
        }

        void x( const scenario<I>& s,
                vector<range<V_t>>& r,
                vector<G>& gh) const
        {
            assert(s.size());
            x(s.size()-1,zeta(s),r,gh);
        }

        /**
         * @brief
         * @param k stage
         * @param barxi history of the random vector up to #k
         * @param r return value - variable ranges
         * @param gh return value - msconstraints
         */

        void x( unsigned int k,
                const O_t& zeta,
                vector<range<V_t>>& v,
                vector<G>& gh) const
        {
            ranges<V_t> r(d[k]);
            msconstraints<G_t> cs(this->barxdim(k));
            this->x_is(k,zeta,r, cs);
            for(unsigned int i=0; i<gh.size(); i++)
                assert(!(gh[i].constantinlast(d[k])));
            v = r;
            gh=cs;
        }

        virtual std::string varname(unsigned int stage, unsigned int i) const
        {
            return varname_is(stage,i);
        }

        const C rho;

        const msproblemstructure d;
        ///@}

        ///@{
        /// @name Interface towards descendants
protected:

    O_t zeta(const scenario<I>& s) const
    {
        assert(s.size());
        return Z(s);
    }
public:
    bool includedinbarx(unsigned int i, unsigned int j, unsigned int k) const
    {
        assert(i<=k);
        assert(j<=d[i]);
        R r;
        return i==k ? true : r.included(i,j,k);
    }
    unsigned int barxdimupto(unsigned int i, unsigned int k) const
    {
        assert(k<d.size());
        unsigned int s = 0;
        for(unsigned int ii=0; ii<=i; ii++ )
            for(unsigned int jj=0; jj<d[ii]; jj++)
                if(includedinbarx(ii,jj,k))
                    s++;
        return s;
    }

    unsigned int tildexdim(unsigned int k) const
    {
        assert(k>0);
        return barxdimupto(k-1, k);
    }

    unsigned int barxdim(unsigned int k) const
    {
        return barxdimupto(k, k);
    }
    unsigned int xdim(unsigned int k) const
    {
        assert(k <= this->T());
        return d[k];
    }



private:

    virtual void f_is(
            unsigned int k,
            const O_t& zeta,
            F& f
            ) const = 0;

    virtual void x_is(
            unsigned int k,
            const O_t& zeta,
            ranges<V_t>& r,
            msconstraints<G>& g
            ) const = 0;

    virtual std::string varname_is(unsigned int stage, unsigned int i) const
    {
        std::ostringstream s;
        s << "x" << stage << "_" << i;
        return s.str();
    }

    ///@}
};

class expectation : public criterion
{
};


class mpmcvar: public criterion
{
public:
    mpmcvar(double alambda, double aalpha=0.05) :
        lambda(alambda), alpha(aalpha)
    {}
    const double lambda;
    const double alpha;
};

class nestedmcvar: public criterion
{
public:
    nestedmcvar(double alambda, double aalpha=0.05) :
        lambda(alambda), alpha(aalpha)
    {}
    const double lambda;
    const double alpha;
};

class msconstraint : virtual public constraint
{
public:
    msconstraint(type t=eq) {}

    virtual bool constantinlast(unsigned int i) const = 0;
};


class linearmsconstraint : public msconstraint, public linearconstraint
{
public:

    linearmsconstraint(unsigned int lhssize) :
       linearconstraint(lhssize), constraint(lhssize)
    {}

    linearmsconstraint(vector<double> lhs, constraint::type t=constraint::eq, double rhs = 0) :
       linearconstraint(lhs,rhs), constraint(lhs.size(),t)
    {
    }

    virtual bool constantinlast(unsigned int d) const
    {
        bool iszero = true;
        for(unsigned int i=0; i < d; i++)
            if(lhs(xdim()-i-1))
            {
                iszero = false;
                break;
            }
        return iszero;
    }
    void set(const vector<double>& lhs, type t, double rhs)
    {
        setlhs(lhs);
        settype(t);
        setrhs(rhs);
    }
};


///@}



} // namespace

#endif // MSPROBLEM_H
