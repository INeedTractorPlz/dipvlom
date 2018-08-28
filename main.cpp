#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <functional>
#include<vector>
#include<iterator>

#include<fstream>
#include<iostream>
#include<ios>

#ifndef BODY
#include"Body_t.hpp"
#endif

#ifndef TYPES
#include"Types.hpp"
#endif

#ifndef FUNCTIONS
#include"functions.hpp"
#endif


using namespace boost::numeric::ublas;

typedef double data_type;
typedef vector<data_type> state_vector;

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

template<typename type, typename Type>
struct RungeKutta4{
    Type k1,k2,k3,k4;
    RungeKutta4(){}
    template<typename sysf>
    void do_step(sysf sysF, const Type &in, Type &out, type time, type step){
        sysF(in,k1,time);
        sysF(in+step*k1/2.,k2,time+step/2.);
        sysF(in+step*k2/2.,k3,time+step/2.);
        sysF(in+step*k3,k4,time+step);
        out=std::move(in + step*(k1+2.*k2+2.*k3+k4)/6.);
    }
};

template<typename type, typename Integrator, typename sysf, typename observer,typename Type, 
typename controller = std::function<int(const Type& X, const Type& Y)>, 
typename exit=std::function<int(const Type&X, const type &t)> >
int integrate(Integrator rk, sysf sysF, Type& X0, type t, type h, int n, observer Obs,
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

int main(int argc, char *const argv[]){
    std::ios_base::sync_with_stdio(false);
    std::ifstream planet_ephemeris_f, planet_mass_f, initial_data_f, size_f, size_integrate_f;
    std::ofstream body_position_f;

    data_type pi = atan(data_type(1.))*4, presicion = 1e-06;
    data_type Time, step, G = 4*pi*pi, const_density_1 = 1.0, const_density_2 = 100., Radius;
    unsigned current = 0, number_steps, number_bodies, number_granulations;
    std::vector<std::vector<state_vector> > planet_ephemeris;
    std::vector<data_type> planet_mass;
    state_vector angles(3);

    RungeKutta4<data_type,Body_position_t> rk4;
    Body_position_t Body_position;

    initial_data_f.open("Body_initial.dat",std::ios_base::in);
        for(unsigned i=0; i<3; ++i)
            initial_data_f >> Body_position.center(i);
        for(unsigned i=0; i<3; ++i)
            initial_data_f >> Body_position.center_velocity(i);
        for(unsigned i=0; i<3; ++i)
            initial_data_f >> angles(i);
        for(unsigned i=0; i<3; ++i)
            initial_data_f >> Body_position.angular_velocity(i);
    initial_data_f.close();
    std::cout << angles << std::endl;
    

    Body_position.Euler_angles = Euler_angles_t(angles);
    std::cout << Body_position.Euler_angles.Rot << std::endl;
    
    size_f.open("size_grid.dat",std::ios_base::in);
        size_f >> number_granulations >> Radius;  
    size_f.close();

    Radius /= 1.5e+11;
    size_integrate_f.open("size_integrate.dat",std::ios_base::in);
        size_integrate_f >> Time >> number_steps;  
    size_integrate_f.close();


    planet_ephemeris_f.open("Planet_ephemeris.dat",std::ios_base::in);
    planet_ephemeris_f >> number_bodies;
    
    data_type void_;
    unsigned t = 0;
    while(!planet_ephemeris_f.eof()){
        planet_ephemeris.push_back(std::vector<state_vector>(number_bodies,state_vector(3,0)));
        for(unsigned i = 0; i < number_bodies; ++i){
            for(unsigned j=0; j < 3; ++j)
                planet_ephemeris_f >> planet_ephemeris[t][i](j);
            for(unsigned j=0; j < 3; ++j)
                planet_ephemeris_f >> void_;
        }
        ++t;
    }
    number_steps = t - 1;
    planet_ephemeris_f.close();
    
    planet_mass_f.open("Planet_mass.dat",std::ios_base::in);
    planet_mass.resize(number_bodies);
    for(unsigned i = 0; i < number_bodies; ++i)
        planet_mass_f >> planet_mass[i];
    planet_mass_f.close();
    
    /*planet_ephemeris.resize(1);
    planet_ephemeris[0].resize(number_bodies);
    planet_ephemeris[0][0] = state_vector(3,0);
    planet_ephemeris[0][1] = state_vector(3,0);
    planet_ephemeris[0][1](0) = 1.0;
    */



    std::function<data_type(const state_vector&)> density = 
    [const_density_1, const_density_2](const state_vector& v)->data_type{
        if(v(1) > 0)
            return const_density_1;
        else
            return const_density_2;
        };
    
    Sphere_t Sphere(Radius);
    Cubic_grid_t grid;
    Body_t Body(density, number_granulations, Body_position, Sphere, grid);
    std::vector<Body_position_t> Rigit_body_orbit;
    Force_t Force(G, current, Body, planet_ephemeris, planet_mass);
    
    Body.grid.grid_fill(Body,Body.Surface);
    Body.calc_mass_and_inertia();
    Body.reduction_to_center(presicion);

    /*for(auto it : Body.points)
        std::cout << it.coord << std::endl;
    */
    std::cout << Force.Full_Torque() << std::endl;
    std::cout << "Body.Mass = " << Body.Mass << std::endl;
    std::cout << "Rotational inertia:" << std::endl << Body.rotational_inertia << std::endl;

    /*Integrator_t Integrator(Rigit_body_orbit, current, Body, Force);
    
    step = Time/number_steps;
    std::cout << "Time = " << Time << std::endl;
    std::cout << "step = " << step << std::endl;
    std::cout << "number_steps = " << number_steps << std::endl;
    std::cout << "planet_ephemeris = " << planet_ephemeris.size() << std::endl;
    
    number_steps = integrate(rk4, Integrator, Body_position, data_type(0.), step, number_steps, Integrator);

    body_position_f.open("Body_position.dat",std::ios_base::trunc);
    body_position_f << number_steps << std::endl;
    for(unsigned k = 0; k < number_steps; ++k){
        body_position_f << Rigit_body_orbit[k].center << std::endl;
        body_position_f << Rigit_body_orbit[k].center_velocity << std::endl;
        body_position_f << Rigit_body_orbit[k].Euler_angles.angles << std::endl;
        body_position_f << Rigit_body_orbit[k].angular_velocity << std::endl << std::endl;
    }
    body_position_f.close();
    */
    return 0;
}