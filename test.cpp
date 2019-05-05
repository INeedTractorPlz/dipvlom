#include<iostream>
#include <boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include<boost/shared_ptr.hpp>
#include"Vector.hpp"

using namespace boost::numeric::ublas;
using namespace boost;

typedef vector<vector<double> > state_type;
template<typename T>
std::vector<T> operator+(std::vector<T> a, std::vector<T> b){
    std::vector<T> c(a.size());
    for(unsigned i=0;i<a.size();++i)
        c[i]=a[i]+b[i];
    return c;
}
template<typename T, typename U>
vector<T> operator/(vector<T> a, U b){
    vector<T> c(a.size());
    for(unsigned i=0;i<a.size();++i)
        c(i)=a(i)/b;
    return c;
}
int main(){
    /*Vector<double>(2,1);
    Vector<Vector<double> > array1(2,Vector<double>(2,1)), array2(2,Vector<double>(2)), sum(2,Vector<double>(2));
    for(unsigned i=0;i<2;++i)
        for(unsigned j=0;j<2;++j){
            array1[i][j]=i+j;
            std::cout << "array1["<< i << "][" << j << "]=" << array1[i][j] << std::endl;
            array2[i][j]=2*i;
            std::cout << "array2["<< i << "][" << j << "]=" << array2[i][j] << std::endl;
        }
    
    sum=array1+array2;
    for(unsigned i=0;i<sum.size();++i)
        for(unsigned j=0;j<sum[i].size();++j){
            std::cout << "sum["<< i << "][" << j << "]=" << sum[i][j] << std::endl;
        }
    sum.resize(4);
    std::cout << "sum.size()= " << sum.size() << std::endl;
    for(unsigned i=0;i<sum.size();++i)
        for(unsigned j=0;j<sum[i].size();++j){
            std::cout << "sum["<< i << "][" << j << "]=" << sum[i][j] << std::endl;
        }
    */    
    Vector<double> v1(1),v2(v1.copy());
    v2(0)=3;
    v1(0)=4;
    std::cout << v1 << std::endl;
    std::cout << v2 << std::endl;

}