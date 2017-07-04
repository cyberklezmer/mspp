#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <cmath>
#include "probability.h"


class meanvardistribution: public discretizabledistribution<std::vector<double>>
{
public:
    meanvardistribution(const std::vector<double>& m,
                        const std::vector<std::vector<double>>& sqV)
        : fm(m), fsqV(sqV)
    {}
    virtual std::vector<std::vector<double>> d(unsigned int n)
    {
        assert(n==1 || n==pow(2,fm.size()));

        std::vector<std::vector<double>> r;
        if(n==1)
            r.push_back(fm);
        else
        {
            std::vector<int> x;
            grid(x,r);
        }
        return r;
    }
private:
    void grid(std::vector<int> x, std::vector<std::vector<double>>& r)
    {
        if(x.size()==fm.size())
        {
            std::vector<double> y(fm.size(),0);
            for(unsigned int i=0; i<fm.size(); i++)
                for(unsigned int j=0; j<fm.size(); j++)
                    y[i] += fsqV[i][j]*x[j];
            for(unsigned int i=0; i<fm.size(); i++)
                y[i] += fm[i];
            r.push_back(y);
        }
        else
        {
            std::vector<int> hi(x), lo(x);
            hi.push_back(1);
            lo.push_back(-1);
            grid(hi,r);
            grid(lo,r);
        }
    }
    const std::vector<double> fm;
    std::vector<std::vector<double>> fsqV;
};

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

    virtual std::vector<std::vector<double>> d(unsigned int n)
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

#endif // DISTRIBUTIONS_H
