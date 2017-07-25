#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <cmath>
#include "probability.h"


class meanvardistribution: virtual public distribution<std::vector<double>>
{
public:
    meanvardistribution(const std::vector<double>& m,
                        const std::vector<std::vector<double>>& sqV,
                        const std::vector<std::vector<double>>& stddists = {})
        : fm(m), fsqV(sqV), fstddists(stddists)
    {
        assert(fm.size() == fsqV.size());
        for(unsigned int i=0; i<fsqV.size(); i++)
            assert(fsqV[i].size()==m.size());
        if(fstddists.size() == 0)
        {
            for(unsigned int i=0; i<m.size(); i++)
                fstddists.push_back({-1,1});
        }
        else
        {
            for(unsigned int i=0; i<m.size(); i++)
            {
                const double tol = 0.0001;
                unsigned int n = fstddists[i].size();
                if(!n)
                    assert(1);
                else if(n>1)
                {
                    double s=0;
                    double s2=0;
                    for(unsigned int j=0; j< n; j++)
                    {
                        double x = fstddists[i][j];
                        s+=x;
                        s2+=x*x;
                    }
                    assert(fabs(s) < tol);
                    assert(fabs(s2-n) < tol);
                }
                else
                    assert(fabs(fstddists[i][0]) < tol);
            }
        }
        assert(fstddists.size() == fm.size());
    }
    virtual std::vector<std::vector<double>> d()
    {
        std::vector<std::vector<double>> r;
        std::vector<double> x;
        grid(x,r);
        return r;
    }
private:
    void grid(std::vector<double> x, std::vector<std::vector<double>>& r)
    {
        unsigned int k=x.size();
        if(k==fm.size())
        {
            std::vector<double> y(k,0);
            for(unsigned int i=0; i<k; i++)
                for(unsigned int j=0; j<k; j++)
                    y[i] += fsqV[i][j]*x[j];
            for(unsigned int i=0; i<k; i++)
                y[i] += fm[i];
            r.push_back(y);
        }
        else
        {
            for(unsigned int i=0; i<fstddists[k].size(); i++)
            {
                std::vector<double> y(x);
                y.push_back(fstddists[k][i]);
                grid(y,r);
            }
        }
    }
    const std::vector<double> fm;
    std::vector<std::vector<double>> fsqV;
    std::vector<std::vector<double>> fstddists;
};

/*
class gaussiandistribution:
         public mcdistribution<std::vector<double>>,
         public discretizabledistribution<std::vector<double>>
{
public:
    gaussiandistribution(const std::vector<double>& m,
                         const std::vector<std::vector<double>>& sqV)
        : fm(m), fsqV(sqV)

    {
        assert(fm.size() == fsqV.size());
        for(unsigned int i=0; i< fm.size(); i++)
            assert(fsqV[i].size()==fm.size());
    }

    virtual std::vector<std::vector<double>> d()
    {
        meanvardistribution mvd(fm,fsqV);
        return mvd.d(n);
    }

    virtual std::vector<double> draw()
    {
        assert(0); // not yet implemented

        return std::vector<double>(fm.size());
    }

private:
    const std::vector<double> fm;
    std::vector<std::vector<double>> fsqV;
};
*/
#endif // DISTRIBUTIONS_H
