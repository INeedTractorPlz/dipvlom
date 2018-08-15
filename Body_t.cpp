#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <functional>
#include<vector>
#include<iterator>

#include<iostream>

#ifndef FUNCTIONS
#include"functions.hpp"
#endif

#ifndef TYPES
#include"basic_types.hpp"
#endif

#ifndef BODY
#include"Body_t.hpp"
#endif


using namespace boost::numeric::ublas;


typedef double data_type;
typedef vector<data_type> state_vector;
typedef matrix<data_type> state_matrix;

grid_t::~grid_t() {}
Surface_t::~Surface_t() {}

void Body_t::calc_mass(){
        Quadrature_t  integrator_mas(*this);
        Mass = integrator_mas([](const masspt_t& masspt)->data_type{ return masspt.mass;});
}
    
void Cubic_grid_t::grid_fill(Body_t& Body, const Surface_t& Surface) const{
        data_type border = Surface.max_width(),
        el_volume = Body.grid_width*Body.grid_width*Body.grid_width;
        masspt_t masspt = {state_vector(3,0), el_volume*Body.density(state_vector(3,0))};
        for(data_type x = -border; x <= border; x += Body.grid_width){
            masspt.coord(0) = x;
            for(data_type y = -border; y <= border; y += Body.grid_width){
                masspt.coord(1) = y;
                for(data_type z = -border; z <= border; z += Body.grid_width){
                    masspt.coord(2) = z;
                    if(Surface.is_inside(masspt)){
                        masspt.mass = el_volume*Body.density(masspt.coord);
                        Body.points.push_back(masspt);
                    }
                }
            }
        }
}    
/*
Euler_angles_t operator+(const Euler_angles_t& Angles_left, const Euler_angles_t& Angles_right){
    return Euler_angles_t (Angles_left.angles + Angles_right.angles);
}


Euler_angles_t operator+=(Euler_angles_t& Angles_left, const Euler_angles_t& Angles_right){
    return (Angles_left = Angles_left + Angles_right);
}

template<typename scalar_t>
Euler_angles_t operator*(const Euler_angles_t& Angles_left, const scalar_t& scalar){
    return Euler_angles_t (Angles_left.angles*scalar);
}

template<typename scalar_t>
Euler_angles_t operator*(const scalar_t& scalar, const Euler_angles_t& Angles_left){
    return Euler_angles_t (Angles_left.angles*scalar);
}

template<typename scalar_t>
Euler_angles_t operator/(const Euler_angles_t& Angles_left, const scalar_t& scalar){
    return Euler_angles_t (Angles_left.angles/scalar);
}

Body_position_t operator+(const Body_position_t& Body_left, const Body_position_t& Body_right){
    return Body_position_t (Body_left.center + Body_right.center, Body_left.center_velocity + 
    Body_right.center_velocity, Body_left.Euler_angles.angles + Body_right.Euler_angles.angles, 
    Body_left.angular_velocity + Body_right.angular_velocity);
}

template<typename scalar_t>
Body_position_t operator*(const Body_position_t& Body_left, const scalar_t& scalar){
    return Body_position_t (Body_left.center*scalar, Body_left.center_velocity*scalar, 
    Body_left.Euler_angles.angles*scalar, Body_left.angular_velocity*scalar);
}


template<typename scalar_t>
Body_position_t operator*(const scalar_t& scalar, const Body_position_t& Body_left){
    return Body_position_t (Body_left.center*scalar, Body_left.center_velocity*scalar, 
    Body_left.Euler_angles.angles*scalar, Body_left.angular_velocity*scalar);
}


template<typename scalar_t>
Body_position_t operator/(const Body_position_t& Body_left, const scalar_t& scalar){
    return Body_position_t (Body_left.center/scalar, Body_left.center_velocity/scalar, 
    Body_left.Euler_angles.angles/scalar, Body_left.angular_velocity/scalar);
}
*/