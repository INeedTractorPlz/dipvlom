#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <functional>
#include<vector>
#include<iterator>

#include<iostream>

#ifndef BODY
#include"Body_t.hpp"
#endif

#ifndef TYPES
#include"basic_types.hpp"
#endif

#ifndef FUNCTIONS
#include"functions.hpp"
#endif

using namespace boost::numeric::ublas;


typedef double data_type;
typedef vector<data_type> state_vector;

void Force_t::operator()(const state_vector& R, state_vector& A) const{
       A = state_vector(R.size());
       for(unsigned i=0; i < planet_ephemeris[current].size(); ++i){
           A += planet_mass[i]*(R - planet_ephemeris[current][i])/(norm_2(R - planet_ephemeris[current][i])
           *norm_2(R - planet_ephemeris[current][i]));
       } 
       A*=G;
}

state_vector Force_t::operator()(const masspt_t& masspt) const{
        state_vector F;
        //std::cout << masspt.coord << std::endl;
        //std::cout << Body.body_position.center << std::endl;
        
        (*this)(masspt.coord + Body.body_position.center, F);
        //std::cout << Body.body_position.Euler_angles.Rot << std::endl;
        F = masspt.mass*prod(Body.body_position.Euler_angles.Rot, F);
        //std::cout << F << std::endl;
        return masspt.mass*cross_product(masspt.coord,F);
}

state_vector Force_t::Full_Torque(){
        Quadrature_t integrator_Tor(Body);
        return integrator_Tor(*this)/Body.Mass;
}

void Force_t::operator()(const Body_position_t& R, Body_position_t& Derivative_Body_Position){
        state_vector L;

        data_type A=1,B=1,C=1;
        
        L = Full_Torque();
        L = prod(trans(R.Euler_angles.Rot),L);
        
        const data_type& p = R.angular_velocity(0);
        const data_type& q = R.angular_velocity(1);
        const data_type& r = R.angular_velocity(2);

        const data_type& sphi = R.Euler_angles.sphi;
        const data_type& spsi = R.Euler_angles.spsi;
        const data_type& stetta = R.Euler_angles.stetta;

        const data_type& cphi = R.Euler_angles.cphi;
        const data_type& cpsi = R.Euler_angles.cpsi;
        const data_type& ctetta = R.Euler_angles.ctetta;

        Derivative_Body_Position.center = Body.body_position.center_velocity; //x,y,z
        (*this)(R.center, Derivative_Body_Position.center_velocity); //Vx,Vy,Vz

        Derivative_Body_Position.Euler_angles.angles(0) =   //psi - прецессия
        (p*sphi + q*cphi)/stetta; 
        
        Derivative_Body_Position.Euler_angles.angles(1) =   //phi - угол собственного вращения 
        r - (p*sphi + q*cphi)*ctetta/stetta; 
        
        Derivative_Body_Position.Euler_angles.angles(2) =   //tetta - нутация 
        p*cphi - q*sphi; 

        
        Derivative_Body_Position.angular_velocity(0) = (L(0) + (B-C)*q*r)/A; //p
        Derivative_Body_Position.angular_velocity(1) = (L(1) + (C-A)*p*r)/B; //q
        Derivative_Body_Position.angular_velocity(2) = (L(2) + (A-B)*q*p)/C; //r
}