#ifndef MSPROBLEM_H
#define MSPROBLEM_H

#include "mspp/linear.h"

namespace mspp
{


/// \addtogroup msproblems Multistage Problems
/// \ingroup problems
/// @{

/// Determines dimensions of decision variables.
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
        for(unsigned int i=0; i<=k; i++)
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

/// \brief Vector-like container of constraints
/// \tparam G constraint type
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

/// Container of decision variables ranges
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
    unsigned int size() const { return frs.size(); }
private:
    vector<range<V>> frs;
};


/// \brief Base abstract class - multistage problem
/// \tparam O risk criterion
/// \tparam F objective
/// \tparam G constraint type
/// \tparam Y random parameter type
/// \tparam V variable type
/// \tparam R decision variables restriction
template<typename O, typename F, typename G,
          typename Y, typename V, typename R=allx>
class msproblem : public object
{
    public:

        using O_t = O;
        using F_t = F;
        using G_t = G;
        using Y_t = Y;
        using V_t = V;
        using R_t = R;


        msproblem(const msproblemstructure& ps, O arho=O()):
             fd(ps), rho(arho)
        {
        }

        msproblem(const vector<unsigned int>& ps, O arho=O()):
             fd(ps), rho(arho)
        {
        }

        ///@{
        /// @name Accessors
        unsigned int T() const { return fd.size()-1; }

        virtual F f(
                unsigned int k,
                const Y_t& xi) const
        {
            F f(this->xdim(k));
            f_is(k,xi,f);
            assert(f.xdim()==xdim(k));
            return f;
        }


        /**
         * @brief
         * @param k stage
         * @param barxi history of the random vector up to #k
         * @param r return value - variable ranges
         * @param gh return value - msconstraints
         */

        void x( unsigned int k,
                const Y_t& xi,
                vector<range<V_t>>& v,
                vector<G>& gh) const
        {
            ranges<V_t> r(fd[k]);
            msconstraints<G_t> cs(this->barxdim(k));
            this->x_is(k,xi,r, cs);
            for(unsigned int i=0; i<gh.size(); i++)
                assert(!(cs[i].constantinlast(fd[k])));
            v = r;
            gh=cs;
        }

        double maxf() const
        {
            return maxf_is();
        }

        double minf() const
        {
            return minf_is();
        }

        std::string varname(unsigned int stage, unsigned int i) const
        {
            std::ostringstream s;
            s << varname_is(stage,i) << "_" << stage;
            return s.str();
        }

        void varnames(vector<vector<std::string>>& n) const
        {
           n.resize(T()+1);
           for(unsigned int i=0; i<=T(); i++)
           {
               unsigned int m=d()[i];
               n[i].resize(m);
               for(unsigned int j=0; j<m; j++)
                   n[i][j] = varname(i,j);
           }
        }

        const O rho;

        const msproblemstructure& d() const { return fd;}
        ///@}

        ///@{
        /// @name Interface towards descendants
public:
    bool includedinbarx(unsigned int i, unsigned int j, unsigned int k) const
    {
        assert(i<=k);
        assert(j<=fd[i]);
        R r;
        return i==k ? true : r.included(i,j,k);
    }
    unsigned int barxdimupto(unsigned int i, unsigned int k) const
    {
        assert(k<fd.size());
        unsigned int s = 0;
        for(unsigned int ii=0; ii<=i; ii++ )
            for(unsigned int jj=0; jj<fd[ii]; jj++)
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
        return fd[k];
    }

private:

    virtual void f_is(
            unsigned int k,
            const Y_t& xi,
            F& f
            ) const = 0;

    virtual void x_is(
            unsigned int k,
            const Y_t& xi,
            ranges<V_t>& r,
            msconstraints<G>& g
            ) const = 0;

    virtual double maxf_is() const
    {
        return max<double>();
    }

    virtual double minf_is() const
    {
        return min<double>();
    }

    virtual std::string varname_is(unsigned int stage, unsigned int i) const
    {
        std::ostringstream s;
        s << "X" << "_" << i;
        return s.str();
    }

    ///@}

    msproblemstructure fd;
};

/// Expectation decision criterion
class expectation : public criterion
{
};

/// Multi-period mean-CVaR
class mpmcvar: public criterion
{
public:
    mpmcvar(double alambda, double aalpha=0.05) :
        lambda(alambda), alpha(aalpha)
    {}
    const double lambda;
    const double alpha;
};

/// Nested mean-CVaR
class nestedmcvar: public criterion
{
public:
    nestedmcvar(double alambda, double aalpha=0.05) :
        lambda(alambda), alpha(aalpha)
    {}
    const double lambda;
    const double alpha;
};

/// \brief Multistage constraint
///
/// Unlike \ref constraint, \p msconstraint needs not to be constant
/// in the ``last-stage'' decision variables
class msconstraint : virtual public constraint
{
public:
    msconstraint(type t=eq) {}

    virtual bool constantinlast(unsigned int i) const = 0;
};

/// Linear multi-stage constraint
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
