#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <functional>
#include<vector>
#include<iterator>

#include<iostream>

#ifndef BODY
#include"Body_t.hpp"
#endif

using namespace boost::numeric::ublas;


typedef double data_type;
typedef vector<data_type> state_vector;
state_vector Force_t::external_potential(unsigned number_body, data_type time) const{
        if(number_body == 1)
                return state_vector(3,0.);
        state_vector Earth(3,0.);
        data_type phi = time*M_PI*2;
        Earth(0) = cos(phi);
        Earth(1) = sin(phi);
        //std::cout << "time = " << time << std::endl;
        //std::cout << "Earth = " << Earth << std::endl;
        return Earth;
}
void Force_t::operator()(const state_vector& R, state_vector& A) const{
       A = state_vector(R.size(),0);
       for(unsigned i=0; i < planet_ephemeris[current].size(); ++i){
           A += planet_mass[i]*(R - planet_ephemeris[current][i])/(norm_2(R - planet_ephemeris[current][i])
           *norm_2(R - planet_ephemeris[current][i])*norm_2(R - planet_ephemeris[current][i]));
       } 
       A*=G;
}


void Force_t::operator()(const state_vector& R, state_vector& A, data_type time){
       A = state_vector(R.size(),0);
       state_vector r;
       for(unsigned i = 0; i < planet_mass.size(); ++i){
           r = external_potential(i, time);
           A += planet_mass[i]*(R - r)/(norm_2(R - r)*norm_2(R - r)*norm_2(R - r));
           //std::cout << "r = " << r << std::endl;
           //std::cout << "R - r = " << R - r << std::endl;
           //std::cout << "norm_2(R - r) = " << norm_2(R - r) << std::endl;
       } 
       //std::cout << "G = " << G << std::endl;
       A*=G;
       //std::cout << "A = " << A << std::endl; 
}

state_vector Force_t::operator()(const masspt_t& masspt) const{
        state_vector F;
        //std::cout << masspt.coord << std::endl;
        //std::cout << Body.body_position.center << std::endl;
        
        (*this)(masspt.coord + Body.body_position.center, F);
        //std::cout << Body.body_position.Euler_angles.Rot << std::endl;
        F = masspt.mass*prod(trans(Body.body_position.Euler_angles.Rot), F);
        //std::cout << F << std::endl;
        return masspt.mass*cross_product(masspt.coord,F);
}

state_vector Force_t::Full_Torque(){
        Quadrature_t integrator_Tor(Body);
        return integrator_Tor(*this)/Body.Mass;
}

void Force_t::operator()(const Body_position_t& R, Body_position_t& Derivative_Body_Position,
data_type time){
        state_vector L;

        const data_type& A = Body.rotational_inertia(0,0);
        const data_type& B = Body.rotational_inertia(1,1);
        const data_type& C = Body.rotational_inertia(2,2);
        
        L = Full_Torque();
        //L = prod(trans(R.Euler_angles.Rot),L);
        
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
        (*this)(R.center, Derivative_Body_Position.center_velocity, time); //Vx,Vy,Vz

        //std::cout << "Derivative_Body_Position.center_velocity = " << Derivative_Body_Position.center_velocity << std::endl;
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