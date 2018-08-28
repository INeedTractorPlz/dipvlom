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

void Body_t::calc_mass_and_inertia(){
        Quadrature_t  calc(*this);
        Mass = calc([](const masspt_t& masspt)->data_type{ return masspt.mass;});
        rotational_inertia = state_matrix(3,3);
        rotational_inertia(0,0) = calc([](const masspt_t& masspt)->data_type{
            return masspt.mass*(masspt.coord(1)*masspt.coord(1) + masspt.coord(2)*masspt.coord(2));});
        
        rotational_inertia(1,1) = calc([](const masspt_t& masspt)->data_type{
            return masspt.mass*(masspt.coord(0)*masspt.coord(0) + masspt.coord(2)*masspt.coord(2));});
        
        rotational_inertia(2,2) = calc([](const masspt_t& masspt)->data_type{
            return masspt.mass*(masspt.coord(0)*masspt.coord(0) + masspt.coord(1)*masspt.coord(1));});

        rotational_inertia(0,1) = calc([](const masspt_t& masspt)->data_type{
            return -masspt.mass*masspt.coord(0)*masspt.coord(1);});
        rotational_inertia(1,0) = rotational_inertia(0,1);

        rotational_inertia(0,2) = calc([](const masspt_t& masspt)->data_type{
            return -masspt.mass*masspt.coord(0)*masspt.coord(2);});
        rotational_inertia(2,0) = rotational_inertia(0,2);

        rotational_inertia(1,2) = calc([](const masspt_t& masspt)->data_type{
            return -masspt.mass*masspt.coord(1)*masspt.coord(2);});
        rotational_inertia(2,1) = rotational_inertia(0,2);
}
    
void Cubic_grid_t::grid_fill(Body_t& Body, const Surface_t& Surface) const{
        data_type i = Body.number_granulations,
        el_volume = Body.grid_width*Body.grid_width*Body.grid_width;
        masspt_t masspt = {state_vector(3,0), el_volume*Body.density(state_vector(3,0))};
        
        for(data_type x = -i/2; x <= i/2; ++x){
            masspt.coord(0) = x*Body.grid_width;
            for(data_type y = -i/2; y <= i/2; ++y){
                masspt.coord(1) = y*Body.grid_width;
                for(data_type z = -i/2; z <= i/2; ++z){
                    masspt.coord(2) = z*Body.grid_width;
                    if(Surface.is_inside(masspt)){
                        masspt.mass = el_volume*Body.density(masspt.coord);
                        Body.points.push_back(masspt);
                    }
                }
            }
        }
}
void Body_t::reduction_to_center(data_type presicion){
    Quadrature_t calc_center(*this);
    state_vector center = calc_center([](const masspt_t& masspt)->state_vector{return masspt.coord*masspt.mass;})/this->Mass;
    for(auto it : points){
        it.coord += -center;
    }
    std::cout << "center = " << center << std::endl;
    
    Jacoby_t<data_type> Jacoby(rotational_inertia, presicion);
    for(auto it : points)
        it.coord = prod(Jacoby.rotation_matrix,it.coord);
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