#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include<boost/lexical_cast.hpp>

#include <functional>
#include<vector>
#include<iterator>
#include<limits>

#include<fstream>
#include<iostream>
#include<ios>

#include"Body_t.hpp"
#include"basic_types.hpp"

using namespace boost::numeric::ublas;

typedef double data_type;
typedef vector<data_type> state_vector;
typedef matrix<data_type> state_matrix;

int main(int argc, char *const argv[]){
    std::ios_base::sync_with_stdio(false);
    std::ifstream planet_ephemeris_f, planet_mass_f, initial_data_f, size_f, size_integrate_f;
    std::ofstream body_position_f, angular_velocity_f, planets_f;

    data_type pi = atan(data_type(1.))*4, presicion = std::numeric_limits<data_type>::epsilon();
    data_type Time, step, G = 4*pi*pi, const_density_1 = 1.0e-09, const_density_2 = 1.0e-09;
    unsigned current = 0, number_steps, number_bodies, number_granulations;
    std::vector<std::vector<state_vector> > planet_ephemeris;
    std::vector<data_type> planet_mass;
    
    //RungeKutta4<data_type,Body_position_t> rk4;
    Lobatto<data_type,Body_position_t> lob3a;
    Body_position_t Body_position;

    Body_position.initial("Body_initial.dat");
    simple_cout(Body_position.Euler_angles.Rot);

    data_type T, h, delta, v = norm_2(Body_position.center_velocity), kappa2 = G, r = norm_2(Body_position.center);
    h = v*v/2 - kappa2/r;
    T = 2*pi*kappa2/sqrt(-2*h)/(-2*h);
    delta = T - (0.75 - floor(0.75/T)*T);
    simple_cout("h = ", h, "T = ", T, "delta = ", delta);
    
    std::vector<struct option> longOpts({
        {"number_steps",required_argument, NULL,0}
    });
    
    std::string short_opts = "G:d:T:";
    Parser_t<3,1> Parser(short_opts, &longOpts[0]);

    std::function<void(data_type&, char*)> fun_data_type = fun_for_parser<data_type>;
    std::function<void(unsigned&, char*)> fun_unsigned = fun_for_parser<unsigned>;
    
    size_f.open("size_grid.dat",std::ios_base::in);
        size_f >> number_granulations;  
    size_f.close();

    size_integrate_f.open("size_integrate.dat",std::ios_base::in);
        size_integrate_f >> Time >> number_steps;  
    size_integrate_f.close();

    planet_ephemeris_f.open("Planet_ephemeris.dat",std::ios_base::in);
    planet_ephemeris_f >> number_bodies;
    
    /*data_type void_;
    unsigned t = 0;
    while(!planet_ephemeris_f.eof()){
        planet_ephemeris.emplace_back(std::vector<state_vector>(number_bodies,state_vector(3,0)));
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
    */
    planet_mass_f.open("Planet_mass.dat",std::ios_base::in);
    planet_mass.resize(number_bodies);
    for(unsigned i = 0; i < number_bodies; ++i)
        planet_mass_f >> planet_mass[i];
    planet_mass_f.close();
    
    Parser.Parser(argc, argv, G, const_density_1, Time, number_steps, fun_data_type, 
    fun_data_type, fun_data_type, fun_unsigned);


    std::function<data_type(const state_vector&)> density = 
    [const_density_1, const_density_2](const state_vector& v)->data_type{
        if(v(1) > 0)
            return const_density_1;
        else
            return const_density_2;
    };
    
    Polygon_t Polygon;
    Polygon.initial("aster.dat");
    Cubic_grid_t grid;
    Body_t Body(density, number_granulations, Body_position, Polygon, grid, presicion);
    std::vector<Body_position_t> Rigit_body_orbit;
    Force_t Force(G, data_type(0.), current, Body, planet_ephemeris, planet_mass);
    
    
    simple_cout("Full_torque = ", Force.Full_Torque());
    simple_cout("Body.Mass = ", Body.Mass);
    

    Integrator_t Integrator(Rigit_body_orbit, current, Body, Force);
    
    step = Time/number_steps;
    simple_cout("Time = ", Time);
    simple_cout("step = ", step);
    simple_cout("number_steps = ", number_steps);
    //simple_cout("planet_ephemeris = ", planet_ephemeris.size());
    //number_steps = 1;

    simple_cout(Body_position);
    number_steps = integrate(lob3a, Integrator, Body_position, data_type(0.), step, number_steps, 
    Integrator);
    simple_cout(Body_position);
    
    body_position_f.open("Body_position.dat",std::ios_base::trunc);
    body_position_f << number_steps << std::endl;
    for(unsigned k = 0; k < number_steps; ++k){
        body_position_f << "time = " << step*k << std::endl;
        body_position_f << Rigit_body_orbit[k] << std::endl;
    }
    body_position_f.close();
    
    angular_velocity_f.open("Angular_velocity.dat", std::ios_base::trunc);
    for(unsigned k = 0; k < number_steps; ++k){
        for(auto x : Rigit_body_orbit[k].angular_velocity)
            angular_velocity_f << x << " ";
        angular_velocity_f << std::endl;
    }

    angular_velocity_f.close();
    
    data_type time = 0;
    planets_f.open("planets.dat", std::ios_base::trunc);
    for(unsigned k = 0; k < number_steps; ++k){
        state_vector Earth(3,0.);
        data_type phi = time*M_PI*2;
        Earth(0) = cos(phi);
        Earth(1) = sin(phi);
        time += step;
        for(auto x : Earth)
            planets_f << x << " ";
        for(auto x : Rigit_body_orbit[k].center)
            planets_f << x << " ";
        planets_f << std::endl;
    }

    planets_f.close();
    
    return 0;
}