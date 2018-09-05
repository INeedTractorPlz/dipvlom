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

#ifndef BODY
#include"Body_t.hpp"
#endif

using namespace boost::numeric::ublas;

typedef double data_type;
typedef vector<data_type> state_vector;
typedef matrix<data_type> state_matrix;

int main(int argc, char *const argv[]){
    std::ios_base::sync_with_stdio(false);
    std::ifstream planet_ephemeris_f, planet_mass_f, initial_data_f, size_f, size_integrate_f;
    std::ofstream body_position_f;

    data_type pi = atan(data_type(1.))*4, presicion = std::numeric_limits<data_type>::epsilon();
    data_type Time, step, G = 4*pi*pi, const_density_1 = 1.0e-04, const_density_2 = 1.0e-04, Radius;
    unsigned current = 0, number_steps, number_bodies, number_granulations;
    std::vector<std::vector<state_vector> > planet_ephemeris;
    std::vector<data_type> planet_mass;
    
    RungeKutta4<data_type,Body_position_t> rk4;
    Body_position_t Body_position;

    std::cout << "STOP_-1" << std::endl;
    Body_position.initial("Body_initial.dat");
    std::cout << Body_position.Euler_angles.Rot << std::endl;

    std::vector<struct option> longOpts({
        {"number_steps",required_argument, NULL,0}
    });
    
    std::string short_opts = "G:d:";
    Parser_t<2,1> Parser(short_opts, &longOpts[0]);

    std::function<void(data_type&, char*)> fun_data_type = fun_for_parser<data_type>;
    std::function<void(unsigned&, char*)> fun_unsigned = fun_for_parser<unsigned>;
    
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

    Parser.Parser(argc, argv, G, const_density_1, number_steps, fun_data_type, 
    fun_data_type, fun_unsigned);


    std::function<data_type(const state_vector&)> density = 
    [const_density_1, const_density_2](const state_vector& v)->data_type{
        if(v(1) > 0)
            return const_density_1;
        else
            return const_density_2;
        };
    
    //Sphere_t Sphere(Radius);
    Polygon_t Polygon;
    Polygon.initial("aster.dat");
    Cubic_grid_t grid;
    std::cout << "STOP_0" << std::endl;
    Body_t Body(density, number_granulations, Body_position, Polygon, grid, presicion);
    std::vector<Body_position_t> Rigit_body_orbit;
    Force_t Force(G, current, Body, planet_ephemeris, planet_mass);
    
    /*for(auto it : Body.points)
        std::cout << it.coord << std::endl;
    */
    std::cout << "STOP_4" << std::endl;
    
    std::cout << Force.Full_Torque() << std::endl;
    std::cout << "Body.Mass = " << Body.Mass << std::endl;
    

    Integrator_t Integrator(Rigit_body_orbit, current, Body, Force);
    
    step = Time/number_steps;
    std::cout << "Time = " << Time << std::endl;
    std::cout << "step = " << step << std::endl;
    std::cout << "number_steps = " << number_steps << std::endl;
    std::cout << "planet_ephemeris = " << planet_ephemeris.size() << std::endl;
    
    number_steps = integrate(rk4, Integrator, Body_position, data_type(0.), step, number_steps, 
    Integrator);

    body_position_f.open("Body_position.dat",std::ios_base::trunc);
    body_position_f << number_steps << std::endl;
    for(unsigned k = 0; k < number_steps; ++k){
        body_position_f << Rigit_body_orbit[k].center << std::endl;
        body_position_f << Rigit_body_orbit[k].center_velocity << std::endl;
        body_position_f << Rigit_body_orbit[k].Euler_angles.angles << std::endl;
        body_position_f << Rigit_body_orbit[k].angular_velocity << std::endl << std::endl;
    }
    body_position_f.close();
    
    std::cout << "G= " << G << std::endl; 
    std::cout << "const_density_1= " << const_density_1 << std::endl; 
    std::cout << "number_steps= " << number_steps << std::endl; 
    return 0;
}