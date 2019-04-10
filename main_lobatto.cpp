#include<iostream>
#include<fstream>
#include<cstdlib>
#include<limits>

//#include<logout.hpp>
#include"runge-kutta.hpp"

using namespace std;

//using vec = array<long double,2>;
using nt = long double;
using Integrator = Lobatto3a<nt,2,4>;
using vec = Integrator::state_vector_t;

//������ �������� �����
ostream& operator<<(ostream& s,const vec& v)
{
    s<<"[";
    for(unsigned int i=0;i<v.size()-1;i++) s<<v[i]<<",";
    s<<v[v.size()-1]<<"]";
    return s;
}

/*ostream& operator<<(ostream& s,const vec& v)
{
    s<<"[";
    for(unsigned int i=0;i<v.size();i++) s<<"|"<<v[i].first<<","<<v[i].second<<"|";
    s<<"]";
    return s;
}*/

class LVP //Lotka-Volterra problem
{
    long double c1_,c2_;
public:
    explicit LVP(nt c1=2.0L,nt c2=1.0L):c1_(c1),c2_(c2){}
    vec operator() (const vec& x,const nt&)
    {
        return vec{x[0]*(x[1]-c1_),x[1]*(c2_-x[0])};
    }
};

/*class LVP2 //Lotka-Volterra problem
{
    long double c1_,c2_;
public:
    explicit LVP2(nt c1=2.0L,nt c2=1.0L):c1_(c1),c2_(c2){}
    vec operator() (const vec& x,const nt&)
    {
        return vec{qp_pair<nt>{x[0].first*(x[0].second-c1_),x[0].second*(c2_-x[0].first)}};
    }
};*/

/*void lvp_q(const vec& q,const vec& p,vec& y,const long double)
{
    y[0]=q[0]*(p[0]-2.0L);
}

void lvp_p(const vec& q,const vec& p,vec& y,const long double)
{
    y[0]=p[0]*(1.0L-q[0]);
}*/

int main()
{
    LVP lvp{};
    vec x={2.0L,2.0L};
    nt epsilon=10*std::numeric_limits<nt>::epsilon();
    nt dt=0.12L;
    Integrator im{};
    ofstream f("out.dat");
    f<<"0 "<<x[0]<<" "<<x[1]<<endl;
    for(unsigned int i=0;i*dt<100;i++)
    {
        x=im.do_step(x,lvp,i*dt,dt,epsilon);
        f<<i*dt<<" "<<x[0]<<" "<<x[1]<<endl;
        cout<<i*dt<<" ";
    }
    cout<<endl;
    f.close();
    cout<<"end point:"<<x<<endl;
    system("asy pict.asy");
    return 0;

    /*LVP2 lvp{};
    vec x={qp_pair<nt>{2.0L,2.0L}};
    nt epsilon=10*std::numeric_limits<nt>::epsilon();
    nt dt=0.12L;
    Integrator im{};
    ofstream f("out.dat");
    f<<"0 "<<x[0].first<<" "<<x[0].second<<endl;
    for(unsigned int i=0;i*dt<10;i++)
    {
        x=im.do_step(x,lvp,i*dt,dt,epsilon);
        f<<i*dt<<" "<<x[0].first<<" "<<x[0].second<<endl;
        cout<<i*dt<<" ";
    }
    cout<<endl;
    f.close();
    cout<<"end point:"<<x<<endl;
    system("asy pict.asy");
    return 0;*/
}
