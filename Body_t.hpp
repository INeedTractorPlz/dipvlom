#define BODY

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <functional>
#include<vector>
#include<iterator>
#include <algorithm>
#include<string>

#include<iostream>
#include"functions.hpp"

using namespace boost::numeric::ublas;


typedef double data_type;
typedef vector<data_type> state_vector;
typedef matrix<data_type> state_matrix;



struct masspt_t{
    state_vector coord;
    data_type mass;
};

struct triangle_t{
    state_vector vertex_1, vertex_2, vertex_3;
    state_vector normal;
    triangle_t(state_vector vertex_1, state_vector vertex_2, state_vector vertex_3) :
    vertex_1(vertex_1), vertex_2(vertex_2), vertex_3(vertex_3), 
    normal(cross_product(vertex_2 - vertex_1,vertex_3 - vertex_2)){}
    bool right_side(const masspt_t& masspt) const{
        return 0 >= inner_prod(normal, masspt.coord - vertex_1);
    }
};

struct Euler_angles_t{
    state_vector angles;
    data_type cpsi, spsi, cphi, sphi, ctetta, stetta;
    state_matrix Rot;
    Euler_angles_t(const state_vector &angles) : angles(angles){
        cpsi = cos(angles(0)); spsi = sin(angles(0));
        cphi = cos(angles(1)); sphi = sin(angles(1));
        ctetta = cos(angles(2)); stetta = sin(angles(2));
        Rot = state_matrix(3,3);
    
        Rot(0,0)=cpsi*cphi-spsi*ctetta*sphi;
        Rot(0,1)=-cpsi*sphi-spsi*ctetta*cphi;
        Rot(0,2)=sphi*stetta;

        Rot(1,0)=spsi*cphi+cpsi*ctetta*sphi;
        Rot(1,1)=-spsi*sphi+cpsi*ctetta*cphi;
        Rot(1,2)=-cpsi*stetta;

        Rot(2,0)=stetta*sphi;
        Rot(2,1)=stetta*cphi;
        Rot(2,2)=ctetta;
    }
    Euler_angles_t() {
        angles.resize(3);
    }
};

struct Body_position_t{
    state_vector center;
    state_vector center_velocity;
    Euler_angles_t Euler_angles;
    state_vector angular_velocity;
    
    Body_position_t(const state_vector& center, const state_vector& center_velocity,
    const state_vector& angles, const state_vector& angular_velocity) : center(center), 
    center_velocity(center_velocity), Euler_angles(angles), angular_velocity(angular_velocity){ }

    Body_position_t(){
        center.resize(3);
        center_velocity.resize(3);
        angular_velocity.resize(3);
    }

    void Change_position(const state_vector& center, const state_vector& center_velocity,
    const state_vector& angles, const state_vector& angular_velocity){
        this->center = center;
        this->center_velocity = center_velocity; 
        this->Euler_angles = Euler_angles_t(angles);
        this->angular_velocity = angular_velocity;
    }
    void initial(const std::string& file_name);
};

struct Body_t;

struct Surface_t{
    virtual data_type max_width() const = 0;
    virtual bool is_inside(const masspt_t& masspt) const = 0;
    virtual ~Surface_t() = 0;
};

struct grid_t{
    virtual void grid_fill(Body_t& Body, const Surface_t& Surface) const = 0;
    virtual ~grid_t() = 0;
};

struct Cubic_grid_t : grid_t{
    void grid_fill(Body_t& Body, const Surface_t& Surface) const;
};

struct Polygon_t : Surface_t{
    std::vector<triangle_t> poligons;
    bool is_inside(const masspt_t& masspt) const{
        for (auto& p:poligons)
            if (!p.right_side(masspt))
                return false;
        return true;
    }
    data_type max_width() const{
        data_type max_width;
        std::for_each(poligons.begin(), poligons.end(),[&max_width](const triangle_t& triangle){
            data_type max_ = std::max({norm_inf(triangle.vertex_1),norm_inf(triangle.vertex_2),
            norm_inf(triangle.vertex_3)});
            if(max_ > max_width) max_width = max_;
        });
        return max_width;
    }
    Polygon_t(const std::vector<triangle_t>& poligons) : poligons(poligons) {}
    Polygon_t(){}
    void initial(const std::string& file_name);
};

struct Sphere_t : Surface_t{
    data_type Radius;
    Sphere_t(data_type Radius) : Radius(Radius){}
    inline bool is_inside(const masspt_t& masspt) const{
        return (norm_2(masspt.coord) <= Radius);
    }
    
    inline data_type max_width() const{
        return Radius;
    }

};

struct Body_t{
    std::vector<masspt_t> points;
    Body_position_t body_position;    
    state_matrix rotational_inertia;
    std::function<data_type(state_vector)> density;
    data_type Mass;
    data_type grid_width, presicion;
    unsigned number_granulations;

    const Surface_t& Surface;
    const grid_t& grid;

    Body_t(const std::function<data_type(state_vector)>& density, unsigned number_granulations,
    const Body_position_t& body_position, const Surface_t& Surface, 
    const grid_t& grid, data_type presicion) : density(density), 
    number_granulations(number_granulations), body_position(body_position), 
    Surface(Surface), grid(grid), presicion(presicion){
        grid_width = 2*Surface.max_width()/number_granulations;
        grid.grid_fill(*this, Surface);
        calc_mass_and_inertia();
        std::cout << "Rotational inertia:" << std::endl << rotational_inertia << std::endl;
        reduction_to_center(presicion);
    }

    void calc_mass_and_inertia();
    void reduction_to_center(data_type presicion);
};

struct Force_t{
    
    data_type G;
    unsigned& current;
    
    Body_t& Body;
    std::vector<std::vector<state_vector> > planet_ephemeris;
    std::vector<data_type> planet_mass;

    Force_t(data_type G, unsigned& current, Body_t& Body, 
    const std::vector<std::vector<state_vector> >& planet_ephemeris, 
    const std::vector<data_type>& planet_mass) : G(G), current(current), Body(Body),
    planet_ephemeris(planet_ephemeris), planet_mass(planet_mass) {}

    void operator()(const state_vector& R, state_vector& A, data_type time);
    void operator()(const state_vector& R, state_vector& A) const;    
    state_vector operator()(const masspt_t& masspt) const;
    state_vector Full_Torque(); 
    state_vector external_potential(unsigned number_body, data_type time) const;
    void operator()(const Body_position_t& R, Body_position_t& Derivative_Body_Position, data_type time);
};

struct Quadrature_t{
    Body_t& body;
    Quadrature_t(Body_t& body) : body(body){ }
    template<typename Force_type>
    auto operator()(const Force_type& f){
        auto result = f(body.points[0]);
        //std::cout << body.points.size() << std::endl;
        //std::cout << "result[0] = " << result << std::endl;
        for(unsigned i=1; i < body.points.size(); ++i){
            result += f(body.points[i]); 
            //if(i <20)
            //std::cout << "result[" << i << "] = " << f(body.points[i]) << std::endl;
        }
        return result;
    }
};

struct Integrator_t{
    std::vector<Body_position_t>& Rigit_body_orbit;
    unsigned& current;
    Body_t& Body;
    Force_t& Force;
    Integrator_t(std::vector<Body_position_t> &Rigit_body_orbit, unsigned& current, Body_t &Body, 
    Force_t& Force) : Rigit_body_orbit(Rigit_body_orbit), current(current), 
    Body(Body), Force(Force){  }
    void operator()(const Body_position_t& R, Body_position_t& A, data_type time){
        Body.body_position = R;
        Force(R, A, time);
        //std::cout << "R.center_velocity: " << R.center_velocity << std::endl;
        //std::cout << "A.center: " << A.center << std::endl;
        //std::cout << "A.center_velocity: " << A.center_velocity << std::endl;
    }
    void operator()(const Body_position_t& R, data_type time){
        std::cout << "current = " << current << std::endl;
        Rigit_body_orbit.push_back(R);        
        ++current;
        std::cout << "R.center: " << R.center << std::endl;
        std::cout << "R.center_velocity: " << R.center_velocity << std::endl;
    }
};

inline Euler_angles_t operator+(const Euler_angles_t& Angles_left, const Euler_angles_t& Angles_right){
    return Euler_angles_t (Angles_left.angles + Angles_right.angles);
}


inline Euler_angles_t operator+=(Euler_angles_t& Angles_left, const Euler_angles_t& Angles_right){
    return (Angles_left = Angles_left + Angles_right);
}

template<typename scalar_t>
inline Euler_angles_t operator*(const Euler_angles_t& Angles_left, const scalar_t& scalar){
    return Euler_angles_t (Angles_left.angles*scalar);
}

template<typename scalar_t>
inline Euler_angles_t operator*(const scalar_t& scalar, const Euler_angles_t& Angles_left){
    return Euler_angles_t (Angles_left.angles*scalar);
}

template<typename scalar_t>
inline Euler_angles_t operator/(const Euler_angles_t& Angles_left, const scalar_t& scalar){
    return Euler_angles_t (Angles_left.angles/scalar);
}

inline Body_position_t operator+(const Body_position_t& Body_left, const Body_position_t& Body_right){
    return Body_position_t (Body_left.center + Body_right.center, Body_left.center_velocity + 
    Body_right.center_velocity, Body_left.Euler_angles.angles + Body_right.Euler_angles.angles, 
    Body_left.angular_velocity + Body_right.angular_velocity);
}

template<typename scalar_t>
inline Body_position_t operator*(const Body_position_t& Body_left, const scalar_t& scalar){
    return Body_position_t (Body_left.center*scalar, Body_left.center_velocity*scalar, 
    Body_left.Euler_angles.angles*scalar, Body_left.angular_velocity*scalar);
}


template<typename scalar_t>
inline Body_position_t operator*(const scalar_t& scalar, const Body_position_t& Body_left){
    return Body_position_t (Body_left.center*scalar, Body_left.center_velocity*scalar, 
    Body_left.Euler_angles.angles*scalar, Body_left.angular_velocity*scalar);
}


template<typename scalar_t>
inline Body_position_t operator/(const Body_position_t& Body_left, const scalar_t& scalar){
    return Body_position_t (Body_left.center/scalar, Body_left.center_velocity/scalar, 
    Body_left.Euler_angles.angles/scalar, Body_left.angular_velocity/scalar);
}