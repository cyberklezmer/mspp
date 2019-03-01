#ifndef MSPPTEST_H
#define MSPPTEST_H

#include <cmath>
#include <fstream>
#include "mspp/lpsolver.h"
#include "mspp/random.h"
#include "mspp/process.h"
#include "mspp/cdists.h"
#include "mspp/de.h"

using namespace mspp;

template <typename V>
class csvlpsolver : public lpsolver<V>
{
public:

    virtual void solve(
            const lpproblem<V>& lp,
            vector<double>& sol,
            double& objvalue) const
    {
        using namespace std;
        ofstream f("problem.csv");
        f << "variables" << endl;
        for(int i=0;i<lp.vars.size(); i++)
            f << lp.varnames[i] << ",";
        f << endl;
        f << endl;

        f << "objective function" << endl;
        for(int i=0;i<lp.f.size(); i++)
           f << lp.f[i] << ",";
        f << endl;

        f << "lower consttraints" << endl;
        for(int i=0;i<lp.vars.size(); i++)
            if(lp.vars[i].islinf())
                f << ",";
            else
                f << lp.vars[i].l() << ",";
        f << endl;

        f << "upper consttraints" << endl;
        for(int i=0;i<lp.vars.size(); i++)
            if(lp.vars[i].ishinf())
                f << ",";
            else
                f << lp.vars[i].h() << ",";

        f << endl;

        for(int y=0; y<3; y++)
        {
            constraint::type t = (constraint::type) y;
            f << "consttraints ";
            switch(t)
            {
                case constraint::eq:
                    f << "=";
                break;
                case constraint::geq:
                    f << ">=";
                break;
                case constraint::leq:
                    f << "<=";
                break;
            }
            f << endl;
            for(int i=0; i<lp.constraints.size(); i++)
            {
                const sparselinearconstraint& c = lp.constraints[i];
                if(c.t() == t)
                {
                    ptr<vector<double>> clhs(c.lhs());

                    int j=0;
                    for(; j<clhs->size(); j++)
                        f << (*clhs)[j] << ",";
                    for(; j<lp.vars.size(); j++)
                        f << "0,";
                    f << c.rhs << endl;
                }
             }
        }
        std::cout << "Test solver only produces problem.csv file.";
    }
};


template <typename T>
inline ostream& operator<<(ostream& o, const std::vector<T>& v)
{
    o << "{";
    if(v.size())
        for(unsigned int j=0; ; j++)
        {
            o << v[j];
            if(j==v.size()-1)
                break;
            o << ",";
        }
    o << "}";

    return o;
}

template <typename T>
inline ostream& operator<<(ostream& o, const vectors<T>& v)
{
    o << "{";
    if(v.size())
    for(unsigned int i=0;;)
    {
        o << v[i];
        if(++i == v.size())
            break;
        o << ",";
    }
    o << "}" << std::endl;
}

template <typename F, typename S>
inline ostream& operator<<(ostream& o, const pair<F,S>& p)
{
    o << "{" << p.first << "," << p.second << "}";
    return o;
}



template <typename X>
class printtreedistcb : public tdcallback<X,ostream>
{
    virtual void callback(const scenario<X>& xi,
                          probability up,
                          ostream* os) const
    {
        assert(xi.size());
        unsigned int k=xi.size()-1;

        for(unsigned int i; i<k; i++)
            *os << '\t';
//        double x = xi[k-1];
        *os << "{" << k << "," << up << "," << xi[k] << "}"  << endl;
    }
};



template <typename P>
inline void printtreed(ostream& o, const P& v)
{
    static_assert(is_base_of<
       treedistribution<typename P::X_t,typename P::K_t,typename P::A_t>,
                   P>::value);
    printtreedistcb<typename P::X_t> cb;
    v.foreachnode(&cb,&o);
}

/*
 * template <typename P>
inline ostream& operator<<(ostream& o, const P& v)
{
//    static_assert(is_base_of<
//       treedistribution<typename P::X_t,typename P::K_t,typename P::A_t>,
//                   P>::value);
//    printtreedistcb<typename P::X_t> cb;
//    v.foreachnode(&cb,&o);
    return o;
}
*/

#endif // MSPPTEST_H
