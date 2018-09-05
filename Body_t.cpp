#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <functional>
#include<vector>
#include<iterator>
#include<string>
#include<fstream>

#include<iostream>

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
        rotational_inertia(2,1) = rotational_inertia(1,2);
}
    
void Cubic_grid_t::grid_fill(Body_t& Body, const Surface_t& Surface) const{
        data_type i = Body.number_granulations,
        el_volume = Body.grid_width*Body.grid_width*Body.grid_width;
        masspt_t masspt = {state_vector(3,0), el_volume*Body.density(state_vector(3,0))};
        
        std::cout << "STOP_1" << std::endl;
        for(data_type x = -i/2; x <= i/2; ++x){
            masspt.coord(0) = x*Body.grid_width;
            for(data_type y = -i/2; y <= i/2; ++y){
                masspt.coord(1) = y*Body.grid_width;
                for(data_type z = -i/2; z <= i/2; ++z){
                    masspt.coord(2) = z*Body.grid_width;
                    if(Surface.is_inside(masspt)){
                        masspt.mass = el_volume*Body.density(masspt.coord);
                        Body.points.emplace_back(masspt);
                    }
                }
            }
        }
        std::cout << "STOP_2" << std::endl;
}
void Body_t::reduction_to_center(data_type presicion){
    Quadrature_t calc_center(*this);
    state_vector center = calc_center([](const masspt_t& masspt)->state_vector{return masspt.coord*masspt.mass;})/this->Mass;
    for(auto it : points){
        it.coord += -center;
    }
    std::cout << "center = " << center << std::endl;
    
    Jacoby_t<data_type> Jacoby(rotational_inertia, presicion);
    std::cout << "Jacoby.rotation_matrix:" << std::endl;
    std::cout << Jacoby.rotation_matrix << std::endl;
    for(auto it : points)
        it.coord = prod(Jacoby.rotation_matrix,it.coord);
}    

void Polygon_t::initial(const std::string& file_name){
        std::ifstream aster_data(file_name);
        std::vector<state_vector> vertices{};
        char l;
        state_vector v(3);
        aster_data >> l;
        while (l == 'v') {
            aster_data >> v;
            std::cout << v << std::endl;
            vertices.push_back(v);
            aster_data >> l;
        }
        std::cout << "number of vertices " << vertices.size() << std::endl;
        while (!aster_data.eof()) {
            unsigned i, j, k;
            aster_data >> i >> j >> k;
            poligons.emplace_back(triangle_t(vertices[i-1], vertices[j-1], vertices[k-1]));
            aster_data >> l;
            std::cout << i << " " << j << " " << k << std::endl;
        }
        std::cout << "number of poligons " << poligons.size() << std::endl;
        aster_data.close();
}

void Body_position_t::initial(const std::string& file_name){
        std::ifstream initial_data_f(file_name);
        state_vector angles(3);
        for(unsigned i=0; i<3; ++i)
            initial_data_f >> center(i);
        for(unsigned i=0; i<3; ++i)
            initial_data_f >> center_velocity(i);
        for(unsigned i=0; i<3; ++i)
            initial_data_f >> angles(i);
        for(unsigned i=0; i<3; ++i)
            initial_data_f >> angular_velocity(i);
        Euler_angles = Euler_angles_t(angles);
        initial_data_f.close();
}
