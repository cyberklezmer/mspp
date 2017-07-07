#include "lpsolvers.h"
#include <ilcplex/ilocplex.h>

using namespace std;

void cplexlpsolver::solve(const varinfo_list& vars,
        const linearfunction& objective,
        const constraint_list<sparselinearconstraint_ptr>& constraints,
        const std::vector<std::string>& varnames,
        std::vector<double>& sol,
        double& objvalue) const
{
    IloEnv env;

    try
    {
        IloModel model(env);

        IloNumVarArray x(env);

        for(int i=0; i<vars.size(); i++)
        {
            const varinfo& v = vars[i];

            const IloNum l = v.l == -inf ? -IloInfinity : v.l;
            const IloNum h = v.h == inf ? IloInfinity : v.h;

            x.add(IloNumVar(env,l,h,IloNumVar::Float,varnames[i].c_str()));
        }

        IloObjective obj = IloMinimize(env);

        assert(objective.coefs.size()==vars.size());

        for(int i=0; i<objective.coefs.size(); i++)
            obj.setLinearCoef(x[i], objective.coefs[i]);

        model.add(obj);

        for(int k=0; k< constraints.size();
            k++)
        {
            IloExpr v(env);

            const sparselinearconstraint& c = *constraints[k];

//            assert(c.lhs.size()==vars.size());

            for(int i=0; i<c.numnz(); i++)
            {
                const sparselinearconstraint::iitem& it = c.nz(i);
                if(it.value) // tbd ensure in lc
                    v+= it.value*x[it.index];
            }
            switch(c.t)
            {
            case linearconstraint::eq:
                model.add(v==c.rhs);
                break;
            case linearconstraint::leq:
                model.add(v<=c.rhs);
                break;
            case linearconstraint::geq:
                model.add(v>=c.rhs);
                break;
            default:
                assert(1);
            }
            v.end();
        }

        IloCplex cplex(model);
        cplex.setOut(sys::log());
        cplex.setWarning(sys::log());
        cplex.exportModel("model.lp"); //lp, mps, sav


       // Optimize the problem and obtain solution.
       int serr;
       if ( !(serr = cplex.solve()) )
       {
          env.error() << "Error " << serr << " when optimizing LP" << endl;
          throw(-1);
       }

       for(int i=0; i<vars.size(); i++)
           sol[i] = cplex.getValue(x[i]);
       objvalue = cplex.getObjValue();
    }
    catch (IloException& e)
    {
       cerr << "Concert exception caught: " << e << endl;
       throw;
    }
    catch (...)
    {
       cerr << "Unknown exception caught" << endl;
       throw;
    }

    env.end();
}
