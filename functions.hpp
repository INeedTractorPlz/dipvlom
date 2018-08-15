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

state_vector cross_product(const state_vector &, const state_vector &);

template<typename point_t>
void reflection(const std::vector<point_t> &points, unsigned axis);

/*template<typename type, typename Integrator, typename sysf, typename observer,typename Type, 
typename controller=std::function<int(const Type& X, const Type& Y)>, 
typename exit=std::function<bool(const Type&X, const type &t)> >
int integrate(Integrator rk, sysf sysF, Type& X0, type t, type h, int n, observer Obs,
                exit Bad = [](const Type& X, const type& t)->int{return 1;},
                controller Err = [](const Type& X, const Type& Y)->bool{return -1;});
*/