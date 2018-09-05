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
#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )

#include<iostream>

#ifndef FUNCTIONS
#include"functions.hpp"
#endif

using namespace boost::numeric::ublas;


typedef double data_type;
typedef vector<data_type> state_vector;
typedef matrix<data_type> state_matrix;

template<int number_shorts, int number_longs, typename Head, typename... Args>
struct parser
{
    static void parse(std::tuple<Head, Args ...> parameters, int opt, char *optarg, 
    const std::vector<char>& short_keys, struct option* longOpts, int longIndex)
    {
        constexpr int k = sizeof...(Args) + 1 - 2*number_longs - 2*number_shorts;
        constexpr int N = k + number_shorts + number_longs;
        auto lambda_get = [](Head head, Args... args){
            return std::make_tuple(std::ref(args)...);
        };
        std::tuple<Args...> parameters_without_Head = std::apply(lambda_get, parameters);
        if(opt == short_keys[k])
            std::get<N>(parameters)(std::get<0>(parameters), optarg); 
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
        auto lambda_get = [](Head head, Args... args){
            return std::make_tuple(std::ref(args)...);
        };
        std::tuple<Args...> parameters_without_Head = std::apply(lambda_get, parameters);
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
void fun_for_parser(type& parameter, char* optarg){
        std::string ss=boost::lexical_cast<std::string>(optarg);
        parameter = boost::lexical_cast<type>(ss);
}

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