#ifndef HMCTEST_H
#define HMCTEST_H

#include "almproblem.h"
#include "mspp/sddp.h"
#include "mspp/de.h"
#include "mspp/cdists.h"
#include "mspp/hmcapprox.h"

namespace mspp
{

class exppmapping: public mapping<pair<double,double>,vector<double>>
{
public:

  virtual vector<double> operator() (const pair<double,double>& a) const
  {
    return { exp(a.first),exp(a.second) };
  }
};

class expmapping: public mapping<double,vector<double>>
{
public:

  virtual vector<double> operator() (const double& a) const
  {
    return { exp(a) };
  }
};


template <typename O, bool cont>
void hmctestetadet(unsigned int T)
{
    assert(T>=1);
    assert(T<=2);

    double xim = -0.3;
    double xisd = 0.5;

    arnormalprocessdistribution xipd(0,xim,xisd,1.0,T);

    onedcovering c1({0+xim},{-0.67+xim, 0.67+xim});
    onedcovering c2({-0.2+2*xim,0.2+2*xim},{-0.33+2*xim,2*xim, 0.33+ 2*xim});

    vector<onedcovering> c({c1});
    if(T==2)
        c.push_back(c2);

    using dhat
       = dhmcapproximation<arnormalprocessdistribution,onedcovering,20>;

    dhat dha(xipd,c);

    using zetat = hmczeta<double,expmapping>;
static_assert(std::is_same<zetat::R_t,vector<double>>::value);
    using apt = almproblem<true, /* nestedmcvar */ mpmcvar>;
    apt ap(0.5,0.05,T);


    desolution<apt,dhat,zetat,O> dx(ap,dha);

    std::cout <<  "obj=." << dx.obj() << std::endl;

    vector<unsigned int> zd(T,1);
    if(!cont)
    {
        msddpsolution<apt,dhat, zetat,O> sx(ap,dha,zd);

        std::cout << "hmctest disc: lb=" <<
             sx.obj().lb() << " < " << "opt=" << dx.obj()
                  << " < ubm=" << sx.obj().ubm() << ", ubb="
                  << sx.obj().ubb() << std::endl;

        if(fabs(sx.obj().lb()-dx.obj()) > 0.05)
        {
            std::cerr << "hmctest disc: opt="
                 << dx.obj() << " expected, " <<
                 sx.obj().lb() << "<" << sx.obj().ubm() <<
                 " achieved." << std::endl;
            throw;
        }
    }
    else
    {
        using chat
           = chmcapproximation<arnormalprocessdistribution,onedcovering>;

        chat cha(xipd,c);

        msddpsolution<apt,chat, zetat,O> cx(ap,cha,zd);

        std::cout << "hmctest cont: lb=" <<
             cx.obj().lb() << " < " << "opt=" << dx.obj()
                  << " < ubm=" << cx.obj().ubm() << ", ubb="
                  << cx.obj().ubb() << std::endl;

        if(fabs(cx.obj().lb()-dx.obj()) > 0.1)
        {
            std::cerr << "malmtest: opt="
                 << dx.obj() << " expected, " <<
                 cx.obj().lb() << "<" << cx.obj().ubm() <<
                 " achieved." << std::endl;
            throw;
        }
    }
    std::cout <<  "passed." << std::endl;
}

template <typename O, bool diraceta = false>
void hmctestetastoch(unsigned int T)
{
    assert(T>=1);
    assert(T<=2);

    double xim = -0.3;
    double xisd = 0.5;

    arnormalprocessdistribution xipd(0,xim,xisd,1.0,T);

    onedcovering c1({0+xim},{-0.67+xim, 0.67+xim});
    onedcovering c2({-0.2+2*xim,0.2+2*xim},{-0.33+2*xim,2*xim, 0.33+ 2*xim});

    vector<onedcovering> c({c1});
    if(T==2)
        c.push_back(c2);

    using apt = almproblem<false, /* nestedmcvar */ mpmcvar>;
    apt ap(0.5,0.05,T);

    double etam = 0;
    double etasd = 0.2;

    using apn = discretization<normaldistribution,nothing>;
    using etapd = iidprocessdistribution<apn>;

    const unsigned int etaaprn = diraceta ? 1 : 20;
    etapd etap( -std::numeric_limits<double>::infinity(),
               apn(normaldistribution(etam,etasd),etaaprn),T+1);

    using dhat
       = dhmcapproximation<arnormalprocessdistribution,onedcovering,20>;

    dhat dha(xipd,c);

    using xietat = xietaprocessdist<dhat,etapd>;
    xietat xieta(dha,etap);

    using zetat = hmczeta<pair<double, double>,exppmapping>;

    vector<unsigned int> zd(T,2);

    msddpsolution<apt,xietat, zetat,O> sx(ap,xieta,zd);

    std::cout << "hmctest stoch disc: lb=" <<
         sx.obj().lb() << " < " << "opt=" << "??"
              << " < ubm=" << sx.obj().ubm() << ", ubb="
              << sx.obj().ubb() << std::endl;

    if(diraceta)
    {
        using zetat = hmczeta<double,expmapping>;

        using apt = almproblem<true, /* nestedmcvar */ mpmcvar>;
        apt ap(0.5,0.05,T);

        desolution<apt,dhat,zetat,O> dx(ap,dha);

        std::cout <<  "obj=." << dx.obj() << std::endl;

        if(fabs(sx.obj().lb()-dx.obj()) > 0.1)
        {
            std::cerr << "hmctest (dirac eta): opt="
                 << dx.obj() << " expected, " <<
                 sx.obj().lb() << "<" << sx.obj().ubm() <<
                 " achieved." << std::endl;
            throw;
        }
    }
    else
    {
        using nt = normaldistribution;
        using etapd = iidprocessdistribution<nt>;

        etapd etap( -std::numeric_limits<double>::infinity(),
                   normaldistribution(etam,etasd),T+1);

        using chat
           = chmcapproximation<arnormalprocessdistribution,onedcovering>;

        chat cha(xipd,c);

        using xietat = xietaprocessdist<chat,etapd>;
        xietat xieta(cha,etap);


        msddpsolution<apt,xietat, zetat,O> cx(ap,xieta,zd);

        std::cout << "hmctest stoch disc: lb=" <<
             cx.obj().lb() << " < " << "opt=" << "??"
                  << " < ubm=" << cx.obj().ubm() << ", ubb="
                  << cx.obj().ubb() << std::endl;

        if(fabs(cx.obj().lb()-sx.obj().lb())>0.05)
        {
            std::cerr << "hmc test: different cont approx values: opt="
                 << cx.obj().lb() << "<" << cx.obj().ubm() << " expected, " <<
                 sx.obj().lb() << "<" << sx.obj().ubm() <<
                 " achieved." << std::endl;
            throw;
        }
    }
    std::cout <<  "passed." << std::endl;

}


} // namespace


#endif // HMCTEST_H
