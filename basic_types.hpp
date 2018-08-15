#define TYPES

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <functional>
#include<vector>
#include<iterator>

#include<iostream>

#ifndef FUNCTIONS
#include"functions.hpp"
#endif

#ifndef BODY
#include"Body_t.hpp"
#endif

using namespace boost::numeric::ublas;


typedef double data_type;
typedef vector<data_type> state_vector;
typedef matrix<data_type> state_matrix;

template<typename type, typename Type>
struct RungeKutta5_Fehlberg{
    Type k1,k2,k3,k4,k5,k6;
    RungeKutta5_Fehlberg(){}
    template<typename sysf>
    void do_step(sysf sysF, const Type &in, Type &out, type time, type step){
        sysF(in,k1,time);
        sysF(in+step*k1/4.,k2,time+step/4.);
        sysF(in+step*(3*k1/32.+9*k2/32.),k3,time+3*step/8);
        sysF(in+step*(1932*k1/2197-7200*k2/2197+7296*k3/2197),k4,time+12*step/13);
        sysF(in+step*(439*k1/216-8*k2+3680*k3/513-845*k4/4104),k5,time+step);
        sysF(in+step*(-8*k1/27+ 2*k2-3544*k3/2565+1859*k4/4104-11*k5/40),k6,time+step/2);
        
        out=std::move(in + step*(16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55));
    }
};