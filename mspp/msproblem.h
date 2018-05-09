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
class msconstraints: public vector<G>
{
public:
    G& add(const G& g)
    {
        vector<G>::push_back(g);
        return *vector<G>::rbegin();
    }
};


/// Base abstract class describing multistage problems
template<typename C, typename F, typename G,
          typename V, typename X, typename Rx=everything, typename Rxi=Rx>
class msproblem : public object
{
    public:

        using C_t = C;
        using F_t = F;
        using G_t = G;
        using V_t = V;
        using X_t = X;
        using Rxi_t = Rxi;
        using Rx_t = Rx;

        msproblem(const msproblemstructure& ps, C rho=C()):
             d(ps), frho(rho)
        {
        }

        msproblem(const vector<unsigned int>& ps, C rho=C()):
             d(ps), frho(rho)
        {
        }

        ///@{
        /// @name Accessors
        unsigned int T() const { return d.size()-1; }

        virtual F f(const scenario<X>& s) const
        {
            assert(s.size());
            return f(s.size()-1,barxi(s));
        }

        virtual F f(
                unsigned int k,
                const subvectors<X>& barxi) const
        {
            F f = f_is(k,barxi);
            assert(f.xdim()==barxdim(k));
            return f;
        }



        void x( const scenario<X>& s,
                ranges<V_t>& r,
                msconstraints<G>& gh) const
        {
            assert(s.size());
            x(s.size()-1,barxi(s),r,gh);
        }

        /**
         * @brief
         * @param k stage
         * @param barxi history of the random vector up to #k
         * @param r return value - variable ranges
         * @param gh return value - msconstraints
         */

        void x( unsigned int k,
                const subvectors<X>& barxi,
                ranges<V>& r,
                msconstraints<G>& gh) const
        {
            r.clear();
            r.resize(d[k]);
            gh.clear();
            this->x_is(k,barxi,r,gh);
            for(unsigned int i=0; i<gh.size(); i++)
                assert(!(gh[i].constantinlast(d[k])));
        }

        virtual std::string varname(unsigned int stage, unsigned int i) const
        {
            return varname_is(stage,i);
        }

        C rho() const { return frho; }

        const msproblemstructure d;
        ///@}

        ///@{
        /// @name Interface towards descendants
protected:

    subvectors<X> barxi(const scenario<X>& s) const
    {
        assert(s.size());
        Rxi r;
        return r(s, s.size()-1);
    }
public:
    bool includedinbarx(unsigned int i, unsigned int j, unsigned int k) const
    {
        assert(i<=k);
        assert(j<=d[i]);
        Rx r;
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
protected:
    /// caution - reference stops to be valid after any change of gs.
    G& addg(msconstraints<G>& gs, unsigned int k) const
    {
        gs.add(G(barxdim(k)));
        return *gs.rbegin();
    }

    F newf(unsigned int k) const
    {
        return F(this->barxdim(k));
    }

private:

    virtual F f_is(
            unsigned int k,
            const subvectors<X>& barxi) const = 0;

    virtual void x_is(
            unsigned int k,
            const subvectors<X>& barxi,
            ranges<V>& r,
            msconstraints<G>& g
            ) const = 0;

    virtual std::string varname_is(unsigned int stage, unsigned int i) const
    {
        std::ostringstream s;
        s << "x" << stage << "_" << i;
        return s.str();
    }

    ///@}

    const C frho;
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
       linearconstraint(lhssize)
    {}

    linearmsconstraint(vector<double> lhs, constraint::type t=constraint::eq, double rhs = 0) :
       linearconstraint(lhs,rhs), constraint(t)
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
