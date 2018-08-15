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
typedef matrix<data_type> state_matrix;

template<typename type>
std::istream& operator>>(std::istream& is, const vector<type,type> in){
    for(auto it : in)
        is >> *it;
    return is;
}

state_vector cross_product(const state_vector &u, const state_vector &v){
    state_vector result(u.size());
    result(0) = u(1)*v(2) - u(2)*v(1);
    result(1) = u(2)*v(0) - u(0)*v(2);
    result(2) = u(0)*v(1) - u(1)*v(0);
    return result;
}

template<typename point_t>
void reflection(const std::vector<point_t> &points, unsigned axis){
    std::vector<point_t> points_inter_x(points);
    std::for_each(points_inter_x.begin(), points_inter_x.end(), 
    [axis](point_t point){point.coord(axis) *= -1;});
    points.insert(points.end(), std::make_move_iterator(points_inter_x.begin()), 
    std::make_move_iterator(points_inter_x.end()));
}
