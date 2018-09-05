#define FUNCTIONS

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <functional>

#include<iostream>

#ifndef TYPES
#include"basic_types.hpp"
#endif

using namespace boost::numeric::ublas;


typedef double data_type;
typedef vector<data_type> state_vector;

template<typename type, typename Integrator, typename sysf, typename observer,typename Type, 
typename controller = std::function<int(const Type& X, const Type& Y)>, 
typename exit=std::function<int(const Type&X, const type &t)> >
inline int integrate(Integrator rk, sysf sysF, Type& X0, type t, type h, int n, observer Obs,
                exit Bad = [](const Type& X, const type& t)->int{return 1;},
                controller Err = [](const Type& X, const Type& Y)->int{return -1;}){
    int state, number_steps = 0;
    Type Y;
    type h_begin=h;
    type T=t+h*n;
    Obs(X0,t);
    while(t<T && Bad(X0,t)){
        state=-1;
        std::cout << t << '\r';
        rk.do_step(sysF,X0,Y,t,h);
        if(Err(X0,Y)==-1){
            goto obs_point;
        } 
        while(Err(X0,Y)){
            h=h/2.;
            rk.do_step(sysF, X0,Y,t,h);
            //std::cout << "h_decrease= " << h << std::endl;
            //std::cout << "Y :" << std::endl << Y << std::endl;
            state=1;
        }
        if(state==1){
            //std::cout << "STOP1" << std::endl;
            h=h*2;
            //std::cout << "h_increase= " << h << std::endl;
            rk.do_step(sysF,X0,Y,t,h);
            goto obs_point;
        }
        while(!Err(X0,Y)){
            if(h > h_begin){
                state=0;
                break;
            }
            //std::cout << "STOP" << std::endl;
            h=h*2;
            //std::cout << "h_increase= " << h << std::endl;
            rk.do_step(sysF, X0,Y,t,h);
            state=0;
        }
        if(state==0){
            h=h/2;
            //std::cout << "h_decrease= " << h << std::endl; 
            rk.do_step(sysF,X0,Y,t,h);
        }
        obs_point:
        X0=Y;
        t+=h;
        ++number_steps;
        Obs(X0,t);
    }
    std::cout << std::endl;
    return number_steps;
}

state_vector cross_product(const state_vector &, const state_vector &);

template<typename point_t>
inline void reflection(const std::vector<point_t> &points, unsigned axis){
    std::vector<point_t> points_inter_x(points);
    std::for_each(points_inter_x.begin(), points_inter_x.end(), 
    [axis](point_t point){point.coord(axis) *= -1;});
    points.insert(points.end(), std::make_move_iterator(points_inter_x.begin()), 
    std::make_move_iterator(points_inter_x.end()));
}

template<typename type>
inline std::istream& operator>>(std::istream& is, vector<type>& in){
    for(auto& it : in)
        is >> it;
    return is;
}