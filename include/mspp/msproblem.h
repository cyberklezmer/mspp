#ifndef MSPROBLEM_H
#define MSPROBLEM_H

#include "mspp/commons.h"

namespace mspp
{



/// \addtogroup msproblems Multistage Problems
/// \ingroup problems
/// @{


class msproblemstructure: public std::vector<unsigned int>
{
public:
    msproblemstructure() {}
    msproblemstructure(unsigned int dim):std::vector<unsigned int>(dim) {}
    msproblemstructure(const std::vector<unsigned int>& v ) :
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

template <typename G>
class msconstraints: public std::vector<G>
{
public:
    G& add(const G& g)
    {
        std::vector<G>::push_back(g);
        return *std::vector<G>::rbegin();
    }
};


/// Base abstract class describing multistage problems
template<typename V, typename F, typename G, typename C>
class msproblem : public object
{
    public:

        using V_t = V;
        using G_t = G;
        using F_t = F;
        using C_t = C;

        msproblem(const msproblemstructure& ps):d(ps)
        {
        }

        msproblem(const std::vector<unsigned int>& psv):d(psv)
        {
        }


        ///@{
        /// @name Accessors
        unsigned int T() const { return d.size()-1; }

        virtual F f(
                unsigned int k,
                const scenario<typename C::X_t>& s) const
        {
            return f(k,C(s));
        }

        virtual F f(
                unsigned int k,
                const C& barxi) const
        {
            return f_is(k,barxi);
        }


        /**
         * @brief
         * @param k stage
         * @param barxi history of the random vector up to #k
         * @param r return value - variable ranges
         * @param gh return value - msconstraints
         */
        void stageinfo(
                unsigned int k,
                const scenario<typename C::X_t>& s,
                vardefs<V>& r,
                msconstraints<G>& gh) const
        {
            stageinfo(k,C(s),r,gh);
        }

        void stageinfo(
                unsigned int k,
                const C& barxi,
                vardefs<V>& r,
                msconstraints<G>& gh) const
        {
            r.clear();
            r.resize(d[k]);
            gh.clear();
            this->stageinfo_is(k,barxi,r,gh);
        }

        virtual std::string varname(unsigned int stage, unsigned int i) const
        {
            return varname_is(stage,i);
        }

        const msproblemstructure d;
        ///@}

        ///@{
        /// @name Interface towards descendants
protected:

   G& addg(msconstraints<G>& gs, unsigned int k) const
   {
       gs.add(G(d,k));
       return *gs.rbegin();
   }

private:
        virtual F f_is(
                unsigned int k,
                const C& barxi) const = 0;

        virtual void stageinfo_is(
                unsigned int k,
                const C& barxi,
                vardefs<V>& r,
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


class linearmsconstraint : public constraint
{
public:

    linearmsconstraint(const msproblemstructure& ps,
                                             unsigned int k)
        : flhssize(ps.sum(k)), frhs(0) {}
    linearmsconstraint(unsigned int lhssize, constraint::type t, double rhs )
        : constraint(t), flhssize(lhssize), frhs(rhs) {}

    virtual double lhs(unsigned int k) const = 0;
    unsigned int lhssize() const { return flhssize; }

    double rhs() const { return frhs; }
    void setrhs(double rhs ) { frhs = rhs; }
    virtual void setlhs(unsigned int i, double v) = 0;
private:
    unsigned int flhssize;
    double frhs;
};

class fulllinearconstraint : public linearmsconstraint
{
public:

    fulllinearconstraint( const msproblemstructure& ps, unsigned int k )
        : linearmsconstraint(ps,k), flhs(lhssize(),0) {}

    fulllinearconstraint( const std::vector<double> lhs, constraint::type t,
                                      double rhs)
        : linearmsconstraint(lhs.size(),t,rhs), flhs(lhs) {}

    double lhs(unsigned int k) const { assert(k<flhs.size()); return flhs[k]; }
    void setlhs(const std::vector<double>& lhs)
    {
        assert(flhs.size() == lhs.size());
        flhs = lhs;
    }
    void setlhs(unsigned int i, double v)
    {
        assert(i<flhs.size());
        flhs[i] = v;
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


class interstagelinearconstraint : public linearmsconstraint
{
public:

    interstagelinearconstraint( const msproblemstructure& ps, unsigned int k )
        : linearmsconstraint(ps,k), flast(k ? ps[k-1] : 0,0), fcurrent(ps[k],0) {}

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



} // namespace

#endif // MSPROBLEM_H