#ifndef PROCESS_H
#define PROCESS_H

#include "mspp/random.h"

namespace mspp
{

/// \addtogroup processes Proceses
/// @{

/// \brief Distribution of a stochastic process with
/// a deterministic initial value
/// \tparam D distribution of individual stochastic components
/// \tparam Z mapping transforming scenarios into \p D's condition
///
/// The distribution as a whole is not conditioned
template <typename D, typename Z>
using processdistribution
 =ivdistribution<diracdistribution<typename D::I_t>,D,Z>;


template <typename D, typename Z, typename M>
using mprocessdistribution
 =mivdistribution<diracdistribution<typename D::I_t>,D,Z,M>;

/// \brief I.i.d. process distribution
/// \tparam D distribution
template <typename D>
using iidprocessdistribution
  = processdistribution<D, noxi<typename D::I_t>>;

/// \brief Markov process distribution
/// \tparam D distribution
template <typename D>
using markovprocessdistribution
  = processdistribution<D, lastxi<typename D::I_t>>;

/// \brief Markov process distribution with known marginals
/// \tparam D distribution
template <typename D, typename M>
using mmarkovprocessdistribution
  = mprocessdistribution<D, lastxi<typename D::I_t>, M>;


template<typename X, typename S=void>
class tdcallback : public object
{
public:
    virtual void callback(const scenario<X>& path,
                          probability up, S* state=0) const = 0;
};

/// \brief Discrete distribution with a tree structure.
/// \tparam X type of individual time components
/// \tparam E type of nodes
/// \tparam A type of state variable, used during recursions
///
/// Its content may be gone through only by recursion.
///
template <typename X, typename K=atom<X>, typename A=nothing>
class treedistribution: virtual public vdistribution<X,nothing>
{
public:
    using X_t = X;
    using K_t = K;
    using A_t = A;

    treedistribution(unsigned int T) : fdim(T+1)
    {
        assert(T>0);
    }

    template <typename S=void>
    void foreachnode(const tdcallback<X, S> *callee, S* state=0) const
    {
        vector<K> e;
        scenario<X> s;
        A mystate;
        beforefe(mystate);
        doforeachnode(callee, e, s, state,mystate);
        afterfe(mystate);
    }

    void branches(const vector<K>& e, vector<K>& es) const
    {
        branches_are(e, es);
    }
    unsigned int T() const { return fdim-1; }
private:

    template <typename S=void>
    void doforeachnode(const tdcallback<X,S> *callee,
                       const vector<K>& ae,
                       const scenario<X> as,
                       S* state,
                       A& mystate) const
    {
        unsigned int k=ae.size();
        vector<K> e(ae);
        vector<K> es;

        branches_are(e,es);

        scenario<X> s(as);
        s.resize(k+1);
        e.resize(k+1);
        bool rec = k < this->dim()-1;

        for(unsigned int i=0;i<es.size();i++)
        {
            e[k]=es[i];
            probability up;
            index2sinfo(e,s[k],up,mystate);
            callee->callback(s,up,state);
            if(rec)
                doforeachnode(callee, e, s, state, mystate);
        }
    }

    virtual void branches_are(const vector<K>& e, vector<K>& es) const = 0;
    virtual void beforefe(A&) const {}
    virtual void afterfe(A&) const {}
protected:
    virtual void index2sinfo(const vector<K>& e,
                                X& s,
                                probability& up,
                                A& ) const
    {
       if constexpr(std::is_same<K,atom<X_t>>::value)
       {
           assert(e.size());
           up = 1;
           for(unsigned int i=0; i<e.size(); i++)
               up *= e[i].p;
           s = e[e.size()-1].x;
       }
       else
            assert(0);
    }
private:
    unsigned int fdim;
    virtual unsigned int dim_is() const { return fdim; }
};

/// \brief \ref treedistribution defined by \ref processdistribution
/// \tparam D the components distribution (descendant of \ref fdistribution)
/// \tparam Z the mapping transforming scenarios into the \p D's condition
template <typename D, typename Z>
class fdprocessdistribution :
        public processdistribution<D,Z>,
        public treedistribution<typename D::I_t>
{
public:
    using D_t = D;
    using X_t=typename D::I_t;
    fdprocessdistribution(const X_t& xi0, const D& d, unsigned int T) :
       treedistribution<X_t>(T),
       processdistribution<D, Z> (xi0, d, T+1)
    {
        static_assert(
           std::is_base_of<fdistribution<X_t,typename D::C_t,
                    D::flistdef>,D>::value
                    );
        assert(T>0);
    }
    fdprocessdistribution(const X_t& xi0, const vector<D>& d) :
       treedistribution<X_t>(d.size()),
       processdistribution<D, Z> (xi0, d)
    {
        static_assert(
           std::is_base_of<fdistribution<X_t,typename D::C_t>,D>::value
                    );
    }
private:
    virtual void branches_are(const vector<atom<X_t>>& e,
         vector<atom<X_t>>& es) const
    {
        unsigned int k=e.size();
        assert(this->T());
        scenario<X_t> s;
        for(unsigned int i=0; i<k; i++)
            s.push_back(e[i].x);
        this->atoms(s,es);
    }
    virtual unsigned int dim_is() const
       { return this->T()+1; }
};


/// @} - processes

} // namespace

#endif // PROCESS_H
