#ifndef TYPES
#define TYPES

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <getopt.h>

#include <functional>
#include<vector>
#include<iterator>
#include<cmath>
#include<tuple>
#include <random>
#include <chrono>

#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )

#include<iostream>

#include"functions.hpp"
#include"runge-kutta.hpp"

using namespace boost::numeric::ublas;


typedef double data_type;
typedef vector<data_type> state_vector;
typedef matrix<data_type> state_matrix;

template<int N, typename Tail, typename... Args>
struct Functions_and_Args_t
{
    static void forward_fanda(std::tuple<Args ...> fanda, Tail optarg){
        
        constexpr int size = sizeof...(Args);
        constexpr int index = size - 2*N - 2;
        std::get<index+1>(fanda)(std::get<index>(fanda), optarg);
        Functions_and_Args_t<N-1, Tail, Args...>::forward_fanda(fanda, optarg);
    }
};

template<typename Tail, typename... Args>
struct Functions_and_Args_t<0, Tail, Args...>
{
    static void forward_fanda(std::tuple<Args ...> fanda, Tail optarg){
        
        constexpr int size = sizeof...(Args);
        constexpr int index = size - 2;
        std::get<index+1>(fanda)(std::get<index>(fanda), optarg);
    }
};

template<int number_shorts, int number_longs, typename Head, typename... Args>
struct parser
{
    static void parse(std::tuple<Head, Args ...> parameters, int opt, char *optarg, 
    const std::vector<char>& short_keys, struct option* longOpts, int longIndex)
    {
        constexpr int k = sizeof...(Args) + 1 - 2*number_longs - 2*number_shorts;
        constexpr int N = k + number_shorts + number_longs;
        auto beheading = [](Head head, Args... args){
            return std::make_tuple(std::ref(args)...);
        };
        std::tuple<Args...> parameters_without_Head = std::apply(beheading, parameters);
        if(opt == short_keys[k])
            //if constexpr (sizeof...(Args) > number_longs + number_shorts)
                std::get<N>(parameters)(std::get<0>(parameters), optarg);
            //else
            //    fun_for_parser<tuple_element< 
        parser<number_shorts-1, number_longs, Args...>::parse(parameters_without_Head, 
        opt, optarg, short_keys,longOpts, longIndex);
    }
};

template<int number_longs, typename Head, typename... Args>
struct parser<0, number_longs, Head, Args...>
{
    static void parse(std::tuple<Head, Args ...> parameters, int opt, char *optarg, 
    const std::vector<char>& short_keys,  struct option* longOpts, int longIndex)
    {
        constexpr int N = sizeof...(Args) + 1 - number_longs;
        constexpr int l = N - number_longs;
        auto beheading = [](Head head, Args... args){
            return std::make_tuple(std::ref(args)...);
        };
        std::tuple<Args...> parameters_without_Head = std::apply(beheading, parameters);
        if(opt == 0)
            if(longOpts[l - short_keys.size()].name == longOpts[longIndex].name)
                std::get<N>(parameters)(std::get<0>(parameters), optarg);
        parser<0,number_longs - 1, Args...>::parse(parameters_without_Head, 
        opt, optarg, short_keys,longOpts, longIndex);
    }
};

template<typename Head, typename... Args>
struct parser<0, 0, Head, Args...>
{
    static void parse(std::tuple<Head, Args ...> parameters, int opt, char *optarg, 
    const std::vector<char>& short_keys,  struct option* longOpts, int longIndex)
    {    }
};

template<int number_shorts, int number_longs>
struct Parser_t{
    const char* short_opts;
    std::vector<char> short_keys;
    struct option* longOpts;
    Parser_t(const std::string& short_opts, struct option* longOpts) :
    short_opts(short_opts.c_str()), longOpts(longOpts){
        for(unsigned i = 0; i < short_opts.size(); i += 2){
            short_keys.push_back(short_opts[i]);
        }
    }
    template<typename ... Types>
    void Parser(int argc, char *const argv[], Types&&... Args){
        std::tuple<Types&& ...> parameters = std::make_tuple(std::ref(Args)...);
        int longIndex; 
        auto opt =  getopt_long(argc,argv,short_opts,longOpts,&longIndex);
        while(opt != -1){
            parser<number_shorts, number_longs, Types&&...>::parse(parameters,
            opt, optarg, short_keys,longOpts,longIndex);  
            opt = getopt_long(argc,argv,short_opts,longOpts,&longIndex);  
        }
    }
};

template<typename type>
struct Jacoby_t{
    //const state_matrix& input_matrix;
    matrix<type>& current_matrix;
    unsigned number_rot, matrix_size;
    matrix<type> rotation_matrix;
    vector<type> eigenvalues;
    type precision;

    type NormOffDiagonal(){
        type Norm = 0;
        for(unsigned i = 0; i < matrix_size; ++i)
            for(unsigned j = 0; (j < matrix_size && i != j); ++j)
                Norm += current_matrix(i,j)*current_matrix(i,j);
        return Norm;
    }
    type NormDiagonal(){
        type Norm = 0;
        for(unsigned i = 0; i < matrix_size; ++i)
            Norm += current_matrix(i,i)*current_matrix(i,i);
        return Norm;
    }
    Jacoby_t(matrix<type>& input_matrix, type precision) : current_matrix(input_matrix), 
    precision(precision), matrix_size(input_matrix.size1()){
        type tetta, c, s;
        rotation_matrix = matrix<type>(matrix_size, matrix_size, 0.);
        matrix<type> rot;
        for(unsigned i = 0; i < matrix_size; ++i)
            rotation_matrix(i,i) = 1.0;
        number_rot = 0;
        eigenvalues.resize(matrix_size);
        while(NormOffDiagonal()/NormDiagonal() > precision){
            for(unsigned p = 0; p < matrix_size; ++p){
                for(unsigned q = p + 1; q < matrix_size; ++q){ 
                    tetta = atan2(2*current_matrix(p,q),current_matrix(q,q)-current_matrix(p,p))/2;
                    c = cos(tetta); s = sin(tetta);
                    rot = matrix<type>(matrix_size, matrix_size, 0.);
                    for(unsigned i = 0; i < matrix_size; ++i)
                        rot(i,i) = 1.0;
                    rot(p,p) = c; rot(q,q) = c;
                    rot(p,q) = s; rot(q,p) = -s;
                    current_matrix = prod(trans(rot),current_matrix);
                    current_matrix = prod(current_matrix, rot);
                    rotation_matrix = prod(rotation_matrix, rot);
                }
            }
            std::cout << NormOffDiagonal()/NormDiagonal() << std::endl;
            ++number_rot;
        }
        
        for(unsigned i=0; i < matrix_size; ++i)
            eigenvalues(i) = current_matrix(i,i);
    }

};

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

template<typename type, typename Type>
struct RungeKutta5_Fehlberg{
    Type k1,k2,k3,k4,k5,k6;
    RungeKutta5_Fehlberg(){}
    template<typename sysf>
    void do_step(sysf sysF, const Type &in, Type &out, type time, type step){
        sysF(in,k1,time);
        sysF(in+step*k1/4.,k2,time+step/4.);
        sysF(in+step*(3*k1/32.+9*k2/32.),k3,time+3*step/8);
        sysF(in+step*(1932*k1/2197-7200*k2/2197+7296*k3/2197),k4,time+12*step/13);
        sysF(in+step*(439*k1/216-8*k2+3680*k3/513-845*k4/4104),k5,time+step);
        sysF(in+step*(-8*k1/27+ 2*k2-3544*k3/2565+1859*k4/4104-11*k5/40),k6,time+step/2);
        
        out=std::move(in + step*(16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55));
    }
};

template<typename type, typename Type, typename Regularization_t, std::size_t system_order = 12,std::size_t method_order = 4>
struct Lobatto{
    rg::Lobatto3a<type, system_order, method_order> lobatto3a;
    using state_vector_t = rg::state_vector<type, system_order>;
    type epsilon = 10*std::numeric_limits<type>::epsilon();
    Regularization_t &Regularization;
    unsigned max_iteration_number;

    Lobatto(Regularization_t &Regularization, unsigned max_iteration_number, type step_0) 
    : Regularization(Regularization), max_iteration_number(max_iteration_number){ 
        Regularization.step_0 = step_0;
        Regularization.step_wo_rand = step_0;
        Regularization.step_wo_reg = step_0;
    }
    
    /*template<class RHS>
    auto do_step(const typename type::state_vector_t& y,RHS rhs,const Number_t& t,const Number_t& dt,
    const Number_t& epsilon)
    {
        return this->do_step_impl(a,b,sigma,y,rhs,t,dt,epsilon);
    }*/
    template<typename sysf>
    void do_step(sysf sysF, const Type &in, Type &out, type time, type &step){
        std::function<state_vector_t(const state_vector_t&, const type&)> rhs =
        [&sysF](const state_vector_t &y, const type &t){
            Type out, in = Type(y);
            sysF(in, out, t);
            //simple_cout(y);
            auto out_vec = out.to_vec(); 
            if(out_vec != out_vec){
                simple_cout("System state vector contains nan at time = ", t);
                simple_cout("out_vec:", out_vec);
                exit(0);
            }
            return out_vec;
        };

        if(Regularization.step_check(step, in))
            Regularization.step_change(step);
        
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);

        std::poisson_distribution<int> distribution(Regularization.step_wo_rand/step);
        double random_number = distribution(generator);
        
        if(random_number >= 5 && step < Regularization.step_wo_rand){
            step *= 2;
            Regularization.step_wo_reg *= 2;
            //simple_cout("random_number = ", random_number, 
            //" Regularization.step_wo_rand/step = ", Regularization.step_wo_rand/step);
            //simple_cout("step = ", step);
        }
        
        rg::state_vector<data_type, 12> out_vec;
        //simple_cout("Before do_step");
        reg_point:
        try{
            out_vec = lobatto3a.do_step(in.to_vec(), rhs, time, step, epsilon, max_iteration_number);
        }catch(...){
            step /= 2;
            Regularization.step_wo_reg /= 2;
            goto reg_point;
        }

        //simple_cout("After do_step");
        out = Type(out_vec);
    }
};

template< typename type, typename Type>
struct Regularization_distance_t{
    std::vector<std::vector<state_vector> > &planet_ephemeris;
    unsigned& current;
    type epsilon, min_distance, step_0, step_wo_reg, step_wo_rand, &time, min_min_distance;
    unsigned number_changes = 0;
    std::ostream *distances_f;
    bool file_record = true;

    Regularization_distance_t(std::vector<std::vector<state_vector> > &planet_ephemeris, 
    unsigned& current, type epsilon, type &time, std::ostream *distances_f = NULL): 
    planet_ephemeris(planet_ephemeris),  current(current), epsilon(epsilon), time(time) {
        if(distances_f == NULL)
            file_record = false;
        else
            this->distances_f = distances_f;
        min_min_distance = -1;
    }

    bool step_check(type &step, const Type &in){
        std::vector<type> distances;
       
        for(unsigned i = 0; i < planet_ephemeris[current].size(); i++){
            distances.emplace_back(norm_2(in.center - planet_ephemeris[current][i]));
        }
        min_distance = *std::min_element(distances.begin(), distances.end());
        
        if(file_record){
            *distances_f << time;
            for(auto x : distances)
                *distances_f << " " << x;
            *distances_f << std::endl;
        }

        if(min_distance <= epsilon)
            return true;
        else
            return false;
    }
    void step_change(type &step){
        if(number_changes == 0){
            min_min_distance = min_distance;
        }
        number_changes++;
        
        //simple_cout("min_min_distance=", min_min_distance, "min_distance =",  min_distance);

        if(min_min_distance > min_distance || min_min_distance == -1)
            min_min_distance = min_distance;
        
        step_wo_rand = step_0*min_distance/epsilon;
        step = step_wo_reg*min_distance/epsilon;
    }
};
#endif