#include<iostream>
#include<cmath>
#include <string>
#include<fstream>
#include<sstream>
#include <boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include<boost/shared_ptr.hpp>
#include <boost/numeric/odeint.hpp>
#include <time.h>
#include<ios>
#include"Array.hpp"


using namespace boost::numeric::ublas;
using namespace boost::numeric::odeint;
using namespace boost;

#define G 4*(M_PI)*(M_PI)
typedef Array<vector<double> > state_type;


template<typename T>
shared_ptr<vector<T> > cross_prod(const vector<T>& v1,const vector<T>& v2){
    if(v1.size()!=3 || v2.size()!=3){
        std::cout << "Input error, vector size must be equal 3";
        exit(0);
    }
    shared_ptr<vector<T> > result(new vector<T>(3));
    (*result)(0)=v1(1)*v2(2)-v1(2)*v2(1);
    (*result)(1)=v1(2)*v2(0)-v1(0)*v2(2);
    (*result)(2)=v1(0)*v2(1)-v1(1)*v2(0);
    return result;
}

struct NBody{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times, &energy_result;
    const vector<double>& m;
    std::vector<vector<vector<double> > >& norm_R;
    double &A,&B,&C,&Yleft,&Yright,m1,m2;
    vector<double> Rleft, Rright;
    unsigned &current;
    NBody(std::vector< state_type >& m_states, std::vector< double >& m_times, std::vector<double>& energy_result, 
    const vector<double>& m, double &A, double &B, double &C, double &Yleft, double &Yright,
    unsigned &current,std::vector<vector<vector<double> > >& norm_R)
    : m_states(m_states), m_times(m_times), energy_result(energy_result), m(m), A(A), B(B), C(C), 
    Yleft(Yleft),Yright(Yright), current(current),norm_R(norm_R){
        m1=m(0)*Yright/(Yleft+Yright);
        m2=m(0)*Yleft/(Yleft+Yright);
        Rleft.resize(3); Rright.resize(3);
        Rleft(0)=0; Rleft(1)=-Yleft;Rleft(2)=0;
        Rright(0)=0; Rright(1)=Yright;Rright(2)=0;
     }
    void operator()(const state_type& R, state_type& A,  double t){
        A[0]=vector<double>(12);
        f(A[0],R);
        vector<vector<double> > norm_r(std::move(vector<vector<double> >(R.size(),vector<double>(R.size()))));
        norm(R,norm_r);
        for(unsigned j=1;j<R.size();++j){
            A[j]=vector<double>(6);
            Fvel(A[j],R[j],norm_R[current](j));
            Faccel(A[j],R,norm_R[current](j),j);
        }
    }    
    void operator()( const state_type &R , double t )
    {
        if(t!=0) ++current;
        m_states.push_back(R.copy()); 
        m_times.push_back(t);
        norm_R.push_back(std::move(vector<vector<double> >(R.size(),vector<double>(R.size())))); 
        norm(R,norm_R[current]);
        energy_result.push_back(EnergyIntegral(R,norm_R[current]));
    }
    shared_ptr<matrix<double> > Rotate(double cpsi, double spsi, double cphi,
    double sphi,double ctetta,double stetta);
    void norm(const state_type& R, vector<vector<double> >& norm_r);
    void Fvel(vector<double>& F, const vector<double>& R,const vector<double>& norm_r);
    
    void Faccel(vector<double>& F, const state_type& R, const vector<double>& norm_r,
    unsigned j);
    
    void f(vector<double>& Derivative, const state_type& R);
    double EnergyIntegral(const state_type& R, const vector<vector<double> >& norm_r);
};


shared_ptr<matrix<double> > NBody::Rotate(double cpsi, double spsi, double cphi,
double sphi,double ctetta,double stetta){
    shared_ptr<matrix<double> > Rot(new matrix<double>(3,3));
    
    (*Rot)(0,0)=cpsi*cphi-spsi*ctetta*sphi;
    (*Rot)(0,1)=-cpsi*sphi-spsi*ctetta*cphi;
    (*Rot)(0,2)=sphi*stetta;

    (*Rot)(1,0)=spsi*cphi+cpsi*ctetta*sphi;
    (*Rot)(1,1)=-spsi*sphi+cpsi*ctetta*cphi;
    (*Rot)(1,2)=-cpsi*stetta;

    (*Rot)(2,0)=stetta*sphi;
    (*Rot)(2,1)=stetta*cphi;
    (*Rot)(2,2)=ctetta;
    return Rot;
}

void NBody::norm(const state_type& R, vector<vector<double> >& norm_r){
   double r1, r2, r3;
   for(unsigned i=0;i<R.size();++i){
        for(unsigned j=0;j<i;++j)
            norm_r(i)(j)=norm_r(j)(i);
         for(unsigned j=i+1;j<R.size();++j){
            r1=R[i](0)-R[j](0); r2=R[i](1)-R[j](1); 
            r3=R[i](2)-R[j](2);
            norm_r(i)(j)=sqrt(r1*r1+r2*r2+r3*r3);
        }
    }
}

void NBody::Fvel(vector<double>& F, const vector<double>& R, const vector<double>& norm_r){
    for(int i=0;i<3;i++) F(i)=R(i+3);
}
void NBody::Faccel(vector<double>& F, const state_type& R, const vector<double>& norm_r,
 unsigned j){
    for(int i=0;i<3;i++){
        F(i+F.size()-3)=0;
        for(unsigned k=0;k<j;++k) 
            F(i+F.size()-3)+=m(k)*(R[k](i)-R[j](i))/(norm_r(k)*norm_r(k)*norm_r(k));
        for(unsigned k=j+1;k<R.size();++k) 
            F(i+F.size()-3)+=m(k)*(R[k](i)-R[j](i))/(norm_r(k)*norm_r(k)*norm_r(k));
        F(i+F.size()-3)*=G;
    }
}




void NBody::f(vector<double>& Derivative, const state_type& R){
    vector<double> a(3);
    vector<double> r1(3), r2(3), F1(3), F2(3),L(3);
    vector<double> norm_R1(R.size()),norm_R2(R.size());
    shared_ptr<matrix<double> > Rot;
    state_type R1(R.copy()),R2(R.copy());
    double cpsi=cos(R[0](6)), spsi=sin(R[0](6)), cphi=cos(R[0](7)), sphi=sin(R[0](7)),
    ctetta=cos(R[0](8)), stetta=sin(R[0](8));


    Rot=Rotate(cpsi, spsi, cphi, sphi, ctetta, stetta);    
    r1=prod(*Rot,Rright); r2=prod(*Rot,Rleft);
    R1[0](0)=R[0](0)+r1(0); R1[0](1)=R[0](1)+r1(1); R1[0](2)=R[0](2)+r1(2);
    R2[0](0)=R[0](0)+r2(0); R2[0](1)=R[0](1)+r2(1); R2[0](2)=R[0](2)+r2(2);
    for(unsigned i=0;i<R.size();++i){
        norm_R1(i)=sqrt((R[i](0)-R1[0](0))*(R[i](0)-R1[0](0))+(R[i](1)-R1[0](1))*(R[i](1)-R1[0](1))+
                                                +(R[i](2)-R1[0](2))*(R[i](2)-R1[0](2)));
        norm_R2(i)=sqrt((R[i](0)-R2[0](0))*(R[i](0)-R2[0](0))+(R[i](1)-R2[0](1))*(R[i](1)-R2[0](1))+
                                                +(R[i](2)-R2[0](2))*(R[i](2)-R2[0](2)));
    }
    Faccel(F1,R1,norm_R1,0);
    Faccel(F2,R2,norm_R2,0);
    
    a=(F1*m1+F2*m2)/(m1+m2);
    L=m1*(*cross_prod(r1,F1))+m2*(*cross_prod(r2,F2));
    
    L=prod(trans(*Rot),L);
    
    Derivative(0)=R[0](3); //x
    Derivative(1)=R[0](4); //y
    Derivative(2)=R[0](5); //z
    
    Derivative(3)=a(0); //Vx
    Derivative(4)=a(1); //Vy
    Derivative(5)=a(2); //Vz

    Derivative(6)=(R[0](9)*sphi+R[0](10)*cphi)/stetta; //psi - прецессия
    Derivative(7)=R[0](11)-(R[0](9)*sphi+R[0](10)*cphi)*ctetta/stetta; //phi - угол собственного вращения
    Derivative(8)=R[0](9)*cphi-R[0](10)*sphi; //tetta - нутация

        
    Derivative(9)=(L(0)+(B-C)*R[0](10)*R[0](11))/A; //p
    Derivative(10)=(L(1)+(C-A)*R[0](9)*R[0](11))/B; //q
    Derivative(11)=(L(2)+(A-B)*R[0](10)*R[0](9))/C; //r


    /*if(norm_R[current](0)(1)<0.07){
    std::cout << a << std::endl;
    std::cout << "norm_R1= " << norm_R1(1) << "; norm_R2= " << norm_R2(1) << std::endl;
    std::cout << "L=" << L << std::endl;
    std::cout << Derivative << std::endl;
    }*/
}


const char* filenamestr(const char* s1, unsigned j, const char* s2){
    static std::string ss;
    ss="";
    ss+= s1; ss+=lexical_cast<std::string>(j); ss+=s2;
    return ss.c_str();
}

double NBody::EnergyIntegral(const state_type& R, const vector<vector<double> >& norm_r){
    double E=0,U;
    E+=m(0)*(A*R[0](9)*R[0](9)+B*R[0](10)*R[0](10)+C*R[0](11)*R[0](11))/2;
    for(unsigned i=0;i<m.size();++i){
        U=0;
        for(unsigned j=i+1;j<m.size();++j)
            U+=m(j)*m(i)/(norm_r(i)(j));
        E-=G*U;
        E+=m(i)*(R[i](3)*R[i](3)+R[i](4)*R[i](4)+R[i](5)*R[i](5))/2;
    }    
    return E;
}

template<typename sysf, typename observer, typename Type>
void RungeKutta(sysf sysF, Type& X0, double t, double h, int n, observer Obs){
    Type k1(X0.size()),k2(X0.size()),k3(X0.size()),k4(X0.size());
    for(unsigned i=0;i<n;++i){
        Obs(X0,t);
        sysF(X0,k1,t);
        sysF(X0+h*k1/2.,k2,t+h/2.);
        sysF(X0+h*k2/2.,k3,t+h/2.);
        sysF(X0+h*k3,k4,t+h);
        X0=X0 + h*(k1+2.*k2+2.*k3+k4)/6.;
        t+=h;
    }
    Obs(X0,t);
}




int main(){
    std::ios_base::sync_with_stdio(false);
    
    std::ofstream out,energy;
    std::ifstream file_in,file_size,file_mass,file_initial,file_RBinitial;
    vector<double> m;
    std::vector<state_type> state_result;
    std::vector<double> time_result, energy_result;
    std::vector<vector<vector<double> > >norm_R;
    double A,B,C,Yleft,Yright;
    unsigned current=0;
    unsigned N,M;
    double T,x3,y3,dx3,dy3,e;
    state_type X0,Y;
    runge_kutta4<state_type> rk;
    
    
    
    file_size.open("file_size.dat",std::ios_base::in);
    file_size >> T >> N >> M;
    file_size.close();

    m.resize(M+1);
    X0.resize(M+1);
    
    file_in.open("RigidBody_in.dat",std::ios_base::in);
    file_in >> m(0) >> Yleft >> Yright >> A >> B >> C;
    file_in.close();
    
    X0[0]=vector<double>(12);
    file_RBinitial.open("RigidBody_initial.dat",std::ios_base::in);
    for(int j=0;j<12;++j)
        file_RBinitial >> X0[0](j);    
    file_RBinitial.close();
    
    file_mass.open("file_mass.dat",std::ios_base::in);
    for(unsigned i=1;i<=M;++i){
        file_mass >> m(i);
    }    
    file_mass.close();
    
    file_initial.open("file_initial.dat",std::ios_base::in);
    for(unsigned i=1;i<=M;++i){
        X0[i]=vector<double>(6);
        for(int j=0;j<6;++j)
            file_initial >> X0[i](j);
    }        
    file_initial.close();

    Y=X0;
    for(unsigned i=0; i< Y.size();++i)
        std::cout << Y[i] << std::endl;
    
    NBody nb(state_result, time_result, energy_result, m,A,B,C,Yleft,Yright,current,norm_R);

    RungeKutta(nb , Y , 0., T/N, N, nb);
    //size_t steps=integrate_n_steps( rk , nb , Y , 0., T/N, N, nb);
    for(unsigned i=0; i< Y.size();++i)
        std::cout << Y[i] << std::endl;

    
    for(unsigned j=0;j<=M;++j){
        out.open(filenamestr("NBodyRK_",j,".dat"),std::ios_base::trunc);
        std::cout << filenamestr("NBodyRK_",j,".dat") << std::endl;
        for(unsigned i=0;i<=current;++i){
            out << time_result[i] << " ";
            if(j==0)
                out << norm_R[i](0)(1) << " ";  
            for(unsigned k=0;k<state_result[i][j].size();++k)
                out << state_result[i][j](k) << " ";
            out << std::endl;
        }        
        out.close();
    }

    energy.open("NBodyRK_energy.dat",std::ios_base::trunc);
    for(unsigned i=0;i<=current;++i){
        energy <<  "h(" << i << ")= " << energy_result[i] << std::endl;
    }
    std::cout << fabs(energy_result[0]-energy_result[current]) << std::endl;
    x3=-0.996404;y3=-0.0027899;
    e=0.3;
    dx3=-sqrt(G)*y3/sqrt(1-e*e)/(1-e*(x3+e));
    dy3=sqrt(G)*sqrt(1-e*e)*(x3+e)/(1-e*(x3+e));

    std::cout << "dx3= " << dx3 << "dy3= " << dy3 << std::endl;


    energy.close();
    
}