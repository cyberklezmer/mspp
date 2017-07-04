#ifndef TESTS_H
#define TESTS_H

#include "linear.h"
#include "scenarios.h"
#include "mmpcvar.h"
#include <fstream>

class csvlpsolver : public lpsolver
{
public:

    virtual void solve(const varinfo_list& vars,
            const linearfunction& objective,
            const constraint_list<linearconstraint>& constraints,
            std::vector<double>& x,
            double& ) const
    {
        using namespace std;
        ofstream f("problem.csv");
        f << "variables" << endl;
        for(int i=0;i<vars.size(); i++)
            f << vars[i].n << ",";
        f << endl;
        f << endl;

        f << "objective function" << endl;
        for(int i=0;i<objective.coefs.size(); i++)
           f << objective.coefs[i] << ",";
        f << endl;

        f << "lower consttraints" << endl;
        for(int i=0;i<vars.size(); i++)
            if(vars[i].l == minf)
                f << ",";
            else
                f << vars[i].l << ",";
        f << endl;

        f << "upper consttraints" << endl;
        for(int i=0;i<vars.size(); i++)
            if(vars[i].h == inf)
                f << ",";
            else
                f << vars[i].h << ",";

        f << endl;

        for(int y=0; y<3; y++)
        {
            linearconstraint::type t = (linearconstraint::type) y;
            f << "consttraints ";
            switch(t)
            {
                case linearconstraint::eq:
                    f << "=";
                break;
                case linearconstraint::geq:
                    f << ">=";
                break;
                case linearconstraint::leq:
                    f << "<=";
                break;
            }
            f << endl;
            for(int i=0; i<constraints.size(); i++)
            {
                const linearconstraint& c = constraints[i];
                if(c.t == t)
                {
                    int j=0;
                    for(; j<c.lhs.size(); j++)
                        f << c.lhs[j] << ",";
                    for(; j<vars.size(); j++)
                        f << "0,";
                    f << c.rhs << endl;
                }
             }
        }
        throw "Test solver only produces a csv file";
    }
};



#endif // TESTS_H
