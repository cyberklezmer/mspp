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
void print(std::ostream& o, vector<vector<T>> v)
{
    o << "[";
    for(unsigned int i=0; i<v.size(); i++)
    {
        o << "[";
        if(v[i].size())
            for(unsigned int j=0; ; j++)
            {
                o << v[i][j];
                if(j==v[i].size()-1)
                    break;
                o << ",";
            }
        o << "]";
    }
    o << "]" << std::endl;
}





#endif // MSPPTEST_H
