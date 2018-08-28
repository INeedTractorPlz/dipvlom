#define TYPES

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <functional>
#include<vector>
#include<iterator>
#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )

#include<iostream>

#ifndef FUNCTIONS
#include"functions.hpp"
#endif

#ifndef BODY
#include"Body_t.hpp"
#endif

using namespace boost::numeric::ublas;


typedef double data_type;
typedef vector<data_type> state_vector;
typedef matrix<data_type> state_matrix;

template<typename type>
struct Jacoby_t{
    //const state_matrix& input_matrix;
    matrix<type> current_matrix;
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
    void HalfRotation(type c, type s, unsigned p, unsigned q){
    	vector<type> v(matrix_size);
    	v = column(rotation_matrix,p);
		column(rotation_matrix,p) = v*c - column(rotation_matrix,q)*s;
		column(rotation_matrix,q) = s*v + column(rotation_matrix,q)*c;
    }
    void Rotation(type c, type s, unsigned p, unsigned q){
        vector<type> v1(matrix_size),v2(matrix_size);
        v1 = row(current_matrix,p); v2 = row(current_matrix,q);
		row(current_matrix,p) = v1*c - v2*s;
		row(current_matrix,q) = v2*c + s*v1;

		current_matrix(p,p) = c*c*v1(p) + s*s*v2(q) - 2*c*s*v1(q);
		current_matrix(p,q) = (c*c - s*s)*v1(q) + c*s*(v1(p) - v2(q));		
		current_matrix(q,q) = s*s*v1(p) + c*c*v2(q) + 2*c*s*v1(q);
        current_matrix(q,p) = current_matrix(p,q);
		
        column(current_matrix,p) = row(current_matrix,p);
        column(current_matrix,q) = row(current_matrix,q);
	}
    Jacoby_t(const matrix<type>& input_matrix, type precision) : current_matrix(input_matrix), 
    precision(precision), matrix_size(input_matrix.size1()){
        type tetta, t, c, s;
        rotation_matrix = matrix<type>(matrix_size, matrix_size, 1.0);
        number_rot = 0;
        eigenvalues.resize(matrix_size);

        while(NormOffDiagonal()/NormDiagonal() > precision){
            for(unsigned p = 0; p < matrix_size; ++p){
                for(unsigned q = p; q < matrix_size; ++q){
                    tetta = (current_matrix(q,q)-current_matrix(p,p))/(2*(current_matrix(p,q)));
                    t = sign(tetta)/(fabs(tetta)+sqrt(tetta*tetta+1));
                    c = 1/sqrt(t*t+1); s = t*c;
                    Rotation(c, s, p, q);
                    HalfRotation(c, s, p, q);
                }
            }
            ++number_rot;
        }

        for(unsigned i=0; i < matrix_size; ++i)
            eigenvalues(i) = current_matrix(i,i);
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