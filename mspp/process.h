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

/*
class processdistribution :
    public ivdistribution<diracdistribution<typename D::I_t>,D,Z,M>
{
public:
    using X_t = typename D::I_t;
    using D_t = D;
    using Z_t = Z;
    processdistribution(const X_t& xi0, const D& d, unsigned int T) :
      ivdistribution<diracdistribution<X_t>,D,Z,M>
         (diracdistribution<X_t>(xi0),d,T),
      vdistribution<X_t,novalue>(T+1)
    {
    }
    processdistribution(const X_t& xi0, const vector<D>& d) :
        ivdistribution<diracdistribution<X_t>,D,Z,M>
          (diracdistribution<X_t>(xi0),d)
    {
    }
    unsigned int T() const { return this->dim()-1; }

    const X_t x0() const
    {
        return this->e().x();
    }
};
*/

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
template <typename X, typename E=atom<X>, typename A=novalue>
class treedistribution: virtual public vdistribution<X,novalue>
{
public:
    using X_t = X;
    using E_t = E;
    using A_t = A;

    treedistribution(unsigned int T) : vdistribution<X,novalue>(T+1) {}

    template <typename S=void>
    void foreachnode(const tdcallback<X, S> *callee, S* state=0) const
    {
        vector<E> e;
        scenario<X> s;
        A mystate;
        beforefe(mystate);
        doforeachnode(callee, e, s, state,mystate);
        afterfe(mystate);
    }

    void branches(const vector<E>& e, vector<E>& es) const
    {
        branches_are(e, es);
    }
private:

    template <typename S=void>
    void doforeachnode(const tdcallback<X,S> *callee,
                       const vector<E>& ae,
                       const scenario<X> as,
                       S* state,
                       A& mystate) const
    {
        unsigned int k=ae.size();
        vector<E> e(ae);
        vector<E> es;

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

    virtual void branches_are(const vector<E>& e, vector<E>& es) const = 0;
    virtual void beforefe(A&) const {}
    virtual void afterfe(A&) const {}
protected:
    virtual void index2sinfo(const vector<E>& e,
                                X& s,
                                probability& up,
                                A& ) const
    {
       if constexpr(std::is_same<E,atom<X_t>>::value)
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
    using X_t=typename D::I_t;
    fdprocessdistribution(const X_t& xi0, const D& d, unsigned int T) :
       treedistribution<X_t>(T),
       processdistribution<D, Z> (xi0, d, T),
       vdistribution<X_t,novalue>(T+1)
    {
        static_assert(
           std::is_base_of<fdistribution<X_t,typename D::C_t,
                    D::flistdef>,D>::value
                    );
        assert(T);
    }
    fdprocessdistribution(const X_t& xi0, const vector<D>& d) :
       treedistribution<X_t>(d.size()+1),
       processdistribution<D, Z> (xi0, d),
       vdistribution<X_t,novalue>(d.size()+1)
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
};


/// @} - processes

} // namespace

#endif // PROCESS_H
