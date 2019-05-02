#ifndef BODY
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
#include"runge-kutta.hpp"
#include"force.hpp"
#include"basic_types.hpp"

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

    Body_position_t(const rg::state_vector<data_type, 12> &in_vec){
        *this = Body_position_t();
        rg::state_vector<data_type, 12> vec = in_vec;
        std::copy(vec.begin(), vec.begin()+3, center.begin());
        //simple_cout(center);
        std::copy(vec.begin()+3, vec.begin()+6, center_velocity.begin());
        //simple_cout(center_velocity);
        state_vector angles = state_vector(3);
        std::copy(vec.begin()+6, vec.begin()+9, angles.begin());
        Euler_angles = Euler_angles_t(angles);
        //simple_cout(angles);
        std::copy(vec.begin()+9, vec.end(), angular_velocity.begin());
        //simple_cout(angular_velocity);
    }

    Body_position_t(){
        center.resize(3);
        center_velocity.resize(3);
        angular_velocity.resize(3);
        Euler_angles.angles.resize(3);
    }

    //void fill_vec(){
    //}

    auto to_vec() const{
        rg::state_vector<data_type, 12> vec;
        std::copy(center.begin(), center.end(), vec.begin());
        //simple_cout(vec);
        std::copy(center_velocity.begin(), center_velocity.end(), vec.begin()+3);
        //simple_cout(vec);
        std::copy(Euler_angles.angles.begin(), Euler_angles.angles.end(), vec.begin()+6);
        //simple_cout(vec);
        std::copy(angular_velocity.begin(), angular_velocity.end(), vec.begin()+9);
        //simple_cout(vec);
        return vec;
    }

    void Change_position(const state_vector& center, const state_vector& center_velocity,
    const state_vector& angles, const state_vector& angular_velocity){
        this->center = center;
        this->center_velocity = center_velocity; 
        this->Euler_angles = Euler_angles_t(angles);
        this->angular_velocity = angular_velocity;
    }

    void initial(const std::string& file_name);

    //private:
    //rg::state_vector<data_type, 12> vec;

};

std::ostream& operator<<(std::ostream& is, const Body_position_t& in);

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
        std::cout << "Rotational inertia after reduction to center:" << std::endl << rotational_inertia << std::endl;
        
    }

    void calc_mass_and_inertia();
    void reduction_to_center(data_type presicion);
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

template<typename type>
struct Record_Functions_t{
    type &time;
    Force_t &Force;
    unsigned &current;
    std::function<void(std::ostream *, const Body_position_t&)> standart_and_time,
    angular_velocity, centers_planets, timeline, momentum;

    Record_Functions_t(type &time, unsigned &current, Force_t& Force) : 
    time(time), current(current), Force(Force) {

        this->standart_and_time = [this](std::ostream *file, const Body_position_t& R){
            *file << "time = " << this->time << std::endl;
            *file << R << std::endl;
        };

        this->angular_velocity = [](std::ostream *file, const Body_position_t& R){
            for(auto x : R.angular_velocity)
                *file << x << " ";
            *file << std::endl;
        };

        this->centers_planets = [this](std::ostream *file, const Body_position_t& R){
            for(auto x : R.center)
                    *file << x << " ";
                for(unsigned i=0; i < this->Force.number_bodies; i++){
                    for(auto x : this->Force.planet_ephemeris[this->current][i])
                        *file << x << " ";
                }
                *file << std::endl;
        };

        this->timeline = [this](std::ostream *file, const Body_position_t& R){
            *file << this->time << std::endl;
        };

        this->momentum = [this](std::ostream *file, const Body_position_t& R){
            for(auto x : this->Force.momentum)
                *file << x << " ";    
            *file << std::endl;
        };
    }
};

template<typename tuple_t>
struct Integrator_t{
    std::vector<Body_position_t> *Rigid_body_orbit = NULL;
    unsigned &current, frequency;
    Body_t& Body;
    Force_t& Force;
    tuple_t * fandf = NULL;
    
    //std::ostream *planets_f = NULL;
    
    Integrator_t(unsigned& current, Body_t &Body, Force_t& Force,
    //std::ostream *planets_f = NULL, 
    unsigned frequency = 1,
    std::vector<Body_position_t> *Rigit_body_orbit = NULL, tuple_t *fandf = NULL) :
    current(current), Body(Body), Force(Force){ 
        if(Rigid_body_orbit != NULL)
            this->Rigid_body_orbit = Rigid_body_orbit;

    //    if(planets_f != NULL)
    //        this->planets_f = planets_f;
        
        if(fandf != NULL)
            this->fandf = fandf;
        
        this->frequency = frequency;
    }
    
    void operator()(const Body_position_t& R, Body_position_t& A, data_type time){
        Force.time = time;
        Body.body_position = R;
        Force.fill_planets();
        Force(R, A);
        //std::cout << "R.center_velocity: " << R.center_velocity << std::endl;
        //std::cout << "A.center: " << A.center << std::endl;
        //std::cout << "A.center_velocity: " << A.center_velocity << std::endl;
    }

    void operator()(const Body_position_t& R, data_type time){
        Force.time = time;
        Body.body_position = R;
        //std::cout << "current = " << current << std::endl;
        if(Rigid_body_orbit != NULL)
            Rigid_body_orbit->push_back(R);
        
        Force.momentum(0) = R.angular_velocity(0)*Body.rotational_inertia(0,0);
        Force.momentum(1) = R.angular_velocity(1)*Body.rotational_inertia(1,1);
        Force.momentum(2) = R.angular_velocity(2)*Body.rotational_inertia(2,2);
        
        /*if(planets_f != NULL & current%frequency == 0){
            for(auto x : R.center)
                *planets_f << x << " ";
            for(unsigned i=0; i < Force.number_bodies; i++){
                for(auto x : Force.planet_ephemeris[current][i])
                    *planets_f << x << " ";
            }
            *planets_f << std::endl;
        }*/

        if(fandf != NULL & current%frequency == 0){
            files_and_functions(*fandf, R);
        }

        ++current;
        Force.fill_ephemeris();

        //std::cout << "R.center: " << R.center << std::endl;
        //std::cout << "R.center_velocity: " << R.center_velocity << std::endl;
    }

    template<typename Tail, typename... Args>
    void files_and_functions(std::tuple<Args ...> &tuple_fandf, Tail optarg){
        constexpr int N = std::tuple_size<tuple_t>::value/2 - 1;
        Functions_and_Args_t<N, Tail, Args...>::forward_fanda(tuple_fandf, optarg);
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

inline Body_position_t operator+=(Body_position_t& Body_left, const Body_position_t& Body_right){
    return Body_left = Body_position_t (Body_left.center + Body_right.center, Body_left.center_velocity + 
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


#endif