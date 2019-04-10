#ifndef RUNGE
#define RUNGE

#include<array>
#include<algorithm>
#include<cmath>
#include<utility>
#include<iosfwd>
//#include<type_traits>

namespace rg{
/**
Класс qp_pair создан для разделенных методов Рунге-Кутты. Он определяет операции необходимые
для работы общего алгоритма в случае данного подкласса методов. Класс представляет пару
переменных, к которым операции применяются параллельно и независимо (кроме операций сравнения).
Созданный интерфейс позволяет общему алгоритму воспринимать qp_pair как обычное число.
Операции и их перегрузки определены с превышением для возможного использования в будущем.
**/

template<class Number_t>
class qp_pair: public std::pair<Number_t,Number_t>
{
    using type = qp_pair<Number_t>;
public:
    constexpr qp_pair()=default;
    constexpr qp_pair(const Number_t& number1,const Number_t& number2)
    :std::pair<Number_t,Number_t>(number1,number2)
    {}
    constexpr type operator*(const type& mult) const
    {
        type result=*this;
        result.first*=mult.first;
        result.second*=mult.second;
        return result;
    }
    constexpr type operator*(const Number_t& number) const
    {
        type result=*this;
        result.first*=number;
        result.second*=number;
        return result;
    }
    constexpr type& operator*=(const type& mult)
    {
        this->first*=mult.first;
        this->second*=mult.second;
        return *this;
    }
    constexpr type& operator*=(const Number_t& number)
    {
        this->first*=number;
        this->second*=number;
        return *this;
    }
    constexpr type operator/(const type& div) const
    {
        type result=*this;
        result.first*=div.first;
        result.second*=div.second;
        return result;
    }
    constexpr type operator/(const Number_t& number) const
    {
        type result=*this;
        result.first*=number;
        result.second*=number;
        return result;
    }
    constexpr type& operator/=(const type& div)
    {
        this->first*=div.first;
        this->second*=div.second;
        return *this;
    }
    constexpr type& operator/=(const Number_t& number)
    {
        this->first*=number;
        this->second*=number;
        return *this;
    }
    constexpr type operator+(const type& add) const
    {
        type result=*this;
        result.first+=add.first;
        result.second+=add.second;
        return result;
    }
    constexpr type operator+(const Number_t& number) const
    {
        type result=*this;
        result.first+=number;
        result.second+=number;
        return result;
    }
    constexpr type& operator+=(const type& add)
    {
        this->first+=add.first;
        this->second+=add.second;
        return *this;
    }
    constexpr type& operator+=(const Number_t& number)
    {
        this->first+=number;
        this->second+=number;
        return *this;
    }
    constexpr type operator-(const type& sub) const
    {
        type result=*this;
        result.first-=sub.first;
        result.second-=sub.second;
        return result;
    }
    constexpr type operator-(const Number_t& number) const
    {
        type result=*this;
        result.first-=number;
        result.second-=number;
        return result;
    }
    constexpr type& operator-=(const type& sub)
    {
        this->first-=sub.first;
        this->second-=sub.second;
        return *this;
    }
    constexpr type& operator-=(const Number_t& number)
    {
        this->first-=number;
        this->second-=number;
        return *this;
    }
    constexpr type& operator++()
    {
        this->first++;
        this->second++;
        return *this;
    }
    constexpr type operator++(int)
    {
        type result=*this;
        this->first++;
        this->second++;
        return result;
    }
    constexpr bool operator>(const type& rhs) const
    {
        return this->first>rhs.first||this->second>rhs.second;
    }
    constexpr bool operator<(const type& rhs) const
    {
        return !(*this>rhs);
    }
};

template<class Number_t>
constexpr qp_pair<Number_t> operator*(Number_t number,const qp_pair<Number_t>& mult)
{
    return mult*number;
}

/**
Класс state_vector определяет вектор переменных, используемый в представлений уравнений движения.
Размер определяется параметром system_order во время компиляции. Определенные операции позволяют
воспринимаеть его как вектор евклидового пространства.
**/
template<class Number_t,std::size_t system_order>
class state_vector:public std::array<Number_t,system_order>
{
    using type = state_vector<Number_t,system_order>;
public:
    constexpr type operator*(const Number_t& number) const
    {
        type result=*this;
        for(auto& x:result) x*=number;
        return result;
    }
    constexpr type& operator*=(const Number_t& number)
    {
        for(auto& x:*this) x*=number;
        return *this;
    }
    constexpr type operator+(const type& add) const
    {
        type result=*this;
        for(std::size_t i=0;i<system_order;++i) result[i]+=add[i];
        return result;
    }
    constexpr type& operator+=(const type& add)
    {
        for(std::size_t i=0;i<system_order;++i) this->at(i)+=add[i];
        return *this;
    }
    constexpr type operator-(const type& sub) const
    {
        type result=*this;
        for(std::size_t i=0;i<system_order;++i) result[i]-=sub[i];
        return result;
    }
    constexpr type& operator-=(const type& sub)
    {
        for(std::size_t i=0;i<system_order;++i) this->at(i)-=sub[i];
        return *this;
    }
    constexpr Number_t norm() const
    {
        auto abss=*this;
        for(auto& x:abss) x=std::abs(x);
        return *std::max_element(abss.begin(),abss.end());
    }
};

template<class Number_t,std::size_t system_order>
constexpr state_vector<Number_t,system_order> operator*(Number_t number,const state_vector<Number_t,system_order>& state_vec)
{
    return state_vec*number;
}


/**
Базовый класс методов Рунге-Кутта, определяющий типы с которыми работают конкретные реализации,
вспомогательные функции ks_diff_norm, ks_norm, make_ks для поиска промежуточных
векторов 'ks' (ks_t - тип массива хранящего эти вектора) и основную функцию do_step_impl для
разрешения метода на одном шаге заданной величины. Базовый класс самостоятельно работать не может,
ему требуется наследник реализующий функцию ks_solver.
Встроенные типы:
state_vector_t - используемые интегратором вектор состояния системы;
vector_t - вектор коэффициентов 'a' (узлы сетки на обезразмеренном отрезке времени)
           и вектор коэффициентов 'sigma' (веса квадратурной формулы для соответствующих узлов);
matrix_t - матрица коэффициентов 'b' (коэффициенты для вычисления промежуточных векторов 'ks');
ks_t - массив промежуточных векторов 'ks'.
**/
template<class Number_t,std::size_t system_order,std::size_t method_order,class Solver>
class Runge_Kutta
{
public:
    using state_vector_t = state_vector<Number_t,system_order>;
    using vector_t = std::array<Number_t,method_order>;
    using matrix_t = std::array<vector_t,method_order>;
    using ks_t = std::array<state_vector_t,method_order>;
    Number_t ks_diff_norm(const ks_t& ks1,ks_t& ks2)
    {
        vector_t norms;
        for(std::size_t i=0;i<method_order;i++) ks2[i]-=ks1[i];
        for(std::size_t i=0;i<method_order;i++) norms[i]=ks2[i].norm();
        return *std::max_element(norms.begin(),norms.end());
    }
    Number_t ks_norm(const ks_t& ks)
    {
        vector_t norms;
        for(std::size_t i=0;i<method_order;i++) norms[i]=ks[i].norm();
        return *std::max_element(norms.begin(),norms.end());
    }
    template<class RHS>
    state_vector_t make_ks(const std::size_t limit,const Number_t& ai,const vector_t& bi,const ks_t& ks,const state_vector_t& y,RHS rhs,const Number_t& t,const Number_t& dt)
    {
        auto v=y;
        for(std::size_t o=0;o<limit;o++) v+=ks[o]*(bi[o]*dt);
        return rhs(v,ai*dt+t);
    }
    template<class RHS>
    state_vector_t do_step_impl(const vector_t& a,const matrix_t& b,const vector_t& sigma,const state_vector_t& y,RHS rhs,const Number_t& t,const Number_t& dt,const Number_t& epsilon)
    {
        ks_t ks=static_cast<Solver*>(this)->ks_solver(a,b,y,rhs,t,dt,epsilon);
        auto v=y;
        for(std::size_t i=0;i<method_order;i++) v+=ks[i]*(sigma[i]*dt);
        return v;
    }
};


/**
Наследник базового класса Runge-Kutta, реализующий функцию ks_solver для неявных методов Рунге-Кутты.
**/
template<class Number_t,std::size_t system_order,std::size_t method_order>
class Runge_Kutta_implicit: public Runge_Kutta<Number_t,system_order,method_order,Runge_Kutta_implicit<Number_t,system_order,method_order>>
{
    using type = Runge_Kutta<Number_t,system_order,method_order,Runge_Kutta_implicit<Number_t,system_order,method_order>>;
public:
    template<class RHS>
    typename type::ks_t ks_solver(const typename type::vector_t& a,const typename type::matrix_t& b,const typename type::state_vector_t& y,RHS rhs,const Number_t& t,const Number_t& dt,const Number_t& epsilon)
    {
        typename type::ks_t ks{},ks_prev{};
        ks.fill(rhs(y,t));
        //auto ks_prev=ks;
        //++ks_prev[0][0]; //set initial difference to allow next cycle
        do
        {
            ks_prev=ks;
            for(std::size_t i=0;i<method_order;i++) ks[i]=this->make_ks(method_order,a[i],b[i],ks,y,rhs,t,dt);
        }while(this->ks_diff_norm(ks,ks_prev)>epsilon*this->ks_norm(ks));
        return ks;
    }
};


/**
Наследник базового класса Runge-Kutta, реализующий функцию ks_solver для полуявных методов Рунге-Кутты.
**/
template<class Number_t,std::size_t system_order,std::size_t method_order>
class Runge_Kutta_semiexplicit: public Runge_Kutta<Number_t,system_order,method_order,Runge_Kutta_semiexplicit<Number_t,system_order,method_order>>
{
    using type = Runge_Kutta<Number_t,system_order,method_order,Runge_Kutta_semiexplicit<Number_t,system_order,method_order>>;
public:
    template<class RHS>
    typename type::ks_t ks_solver(const typename type::vector_t& a,const typename type::matrix_t& b,const typename type::state_vector_t& y,RHS rhs,const Number_t& t,const Number_t& dt,const Number_t& epsilon)
    {
        typename type::ks_t ks{};
        ks.fill(rhs(y,t));
        for(std::size_t i=0;i<method_order;i++)
        {
            typename type::state_vector_t ksi{};
            do
            {
                ksi=ks[i];
                ks[i]=this->make_ks(i+1,a[i],b[i],ks,y,rhs,t,dt);
                ksi-=ks[i];
            }while(ksi.norm()>epsilon*ks[i].norm());
        }
        return ks;
    }
};


/**
Наследник базового класса Runge-Kutta, реализующий функцию ks_solver для явных методов Рунге-Кутты.
**/
template<class Number_t,std::size_t system_order,std::size_t method_order>
class Runge_Kutta_explicit: public Runge_Kutta<Number_t,system_order,method_order,Runge_Kutta_explicit<Number_t,system_order,method_order>>
{
    using type = Runge_Kutta<Number_t,system_order,method_order,Runge_Kutta_explicit<Number_t,system_order,method_order>>;
public:
    template<class RHS>
    typename type::ks_t ks_solver(const typename type::vector_t& a,const typename type::matrix_t& b,const typename type::state_vector_t& y,RHS rhs,const Number_t& t,const Number_t& dt,const Number_t& epsilon)
    {
        typename type::ks_t ks{};
        ks[0]=rhs(y,t+a[0]*dt);
        for(std::size_t i=1;i<method_order;i++) ks[i]=this->make_ks(i,a[i],b[i],ks,y,rhs,t,dt);
        return ks;
    }
};


/**
Класс Lobatto3a, определяющий наборы коэффициентов 'a', 'b', 'sigma' для LobattoIIIa
'a' - узлы сетки на обезразмеренном отрезке времени,
'b' - коэффициенты для вычисления промежуточных векторов 'ks' (аналог узлов сетки, только в пространстве, не на отрезке времени),
'sigma' - веса квадратурной формулы для соответствующих узлов.
**/
template<class Number_t,std::size_t system_order,std::size_t method_order>
class Lobatto3a: protected Runge_Kutta_implicit<Number_t,system_order,method_order>
{
    using type = Runge_Kutta_implicit<Number_t,system_order,method_order>;
public:
    constexpr static const typename type::vector_t make_a()
    {
        static_assert(method_order>=2&&method_order<=9,"That method_order is not allowed.");
        if constexpr (method_order==2) return typename type::vector_t{0.0L,1.0L};
        else if constexpr (method_order==3) return typename type::vector_t{0.0L,0.5L,1.0L};
        else if constexpr (method_order==4) return typename type::vector_t{0.0L,2.7639320225002103035908263312687238e-1L,7.2360679774997896964091736687312762e-1L,1.0L};
        else if constexpr (method_order==5) return typename type::vector_t{0.0L,1.7267316464601142810085377187657082e-1L,5.0e-1L,8.2732683535398857189914622812342918e-1L,1.0L};
        else if constexpr (method_order==6) return typename type::vector_t{0.0L,1.1747233803526765357449851302033093e-1L,3.5738424175967745184292450297956046e-1L,6.4261575824032254815707549702043954e-1L,8.8252766196473234642550148697966908e-1L,1.0L};
        else if constexpr (method_order==7) return typename type::vector_t{0.0L,8.488805186071653506398389301626743e-2L,2.6557560326464289309811405904561684e-1L,5.0e-1L,7.3442439673535710690188594095438317e-1L,9.1511194813928346493601610698373257e-1L,1.0L};
        else if constexpr (method_order==8) return typename type::vector_t{0.0L,6.4129925745196692331277119389668282e-2L,2.041499092834288489277446343010234e-1L,3.9535039104876056561567136982732437e-1L,6.0464960895123943438432863017267563e-1L,7.9585009071657115107225536569897659e-1L,9.3587007425480330766872288061033171e-1L,1.0L};
        else if constexpr (method_order==9) return typename type::vector_t{0.0L,5.012100229426992134382737779083102e-2L,1.6140686024463112327705728645432877e-1L,3.1844126808691092064462396564567039e-1L,5.0e-1L,6.815587319130890793553760343543296e-1L,8.3859313975536887672294271354567122e-1L,9.4987899770573007865617262220916898e-1L,1.0L};
    }
    constexpr static const typename type::vector_t make_sigma()
    {
        if constexpr (method_order==2) return typename type::vector_t{0.5L,0.5L};
        else if constexpr (method_order==3) return typename type::vector_t{1.0L/6.0L,2.0L/3.0L,1.0L/6.0L};
        else if constexpr (method_order==4) return typename type::vector_t{8.3333333333333333333333333333333333e-2L,4.1666666666666666666666666666666667e-1L,4.1666666666666666666666666666666667e-1L,8.3333333333333333333333333333333333e-2L};
        else if constexpr (method_order==5) return typename type::vector_t{5.0e-2L,2.7222222222222222222222222222222222e-1L,3.5555555555555555555555555555555555e-1L,2.7222222222222222222222222222222222e-1L,5.0e-2L};
        else if constexpr (method_order==6) return typename type::vector_t{3.3333333333333333333333333333333333e-2L,1.8923747814892349015830640410601233e-1L,2.7742918851774317650836026256065434e-1L,2.7742918851774317650836026256065434e-1L,1.8923747814892349015830640410601233e-1L,3.3333333333333333333333333333333333e-2L};
        else if constexpr (method_order==7) return typename type::vector_t{2.3809523809523809523809523809523809e-2L,1.3841302368078297400535020314503315e-1L,2.1587269060493131170893551114068114e-1L,2.4380952380952380952380952380952381e-1L,2.1587269060493131170893551114068114e-1L,1.3841302368078297400535020314503315e-1L,2.3809523809523809523809523809523809e-2L};
        else if constexpr (method_order==8) return typename type::vector_t{1.7857142857142857142857142857142631e-2L,1.0535211357175301969149603288787757e-1L,1.7056134624175218238212033855387492e-1L,2.0622939732935194078352648570110395e-1L,2.0622939732935194078352648570110581e-1L,1.705613462417521823821203385538733e-1L,1.0535211357175301969149603288787874e-1L,1.7857142857142857142857142857142631e-2L};
        else if constexpr (method_order==9) return typename type::vector_t{1.3888888888888888888888888888888561e-2L,8.2747680780402762523169860014603253e-2L,1.3726935625008086764035280928968679e-1L,1.7321425548652317255756576606985943e-1L,1.8575963718820861678004535147392211e-1L,1.7321425548652317255756576606986015e-1L,1.3726935625008086764035280928968537e-1L,8.2747680780402762523169860014604949e-2L,1.3888888888888888888888888888888561e-2L};
    }
    constexpr static const typename type::matrix_t make_b()
    {
        if constexpr (method_order==2)
            return typename type::matrix_t{
                typename type::vector_t{0.0L,0.0L},
                typename type::vector_t{0.5L,0.5L}};
        else if constexpr (method_order==3)
            return typename type::matrix_t{
                typename type::vector_t{0.0L,0.0L,0.0L},
                typename type::vector_t{5.0L/24.0L, 1.0/3.0L, -1.0L/24.0L},
                typename type::vector_t{1.0L/6.0L, 2.0L/3.0L, 1.0L/6.0L}};
        else if constexpr (method_order==4)
            return typename type::matrix_t{
                typename type::vector_t{0.0L,0.0L,0.0L,0.0L},
                typename type::vector_t{1.1030056647916491413674311390609397e-1L,1.8969943352083508586325688609390603e-1L,-3.3907364229143883777660480779221592e-2L,1.0300566479164914136743113906093969e-2L},
                typename type::vector_t{7.3032766854168419196590219427239366e-2L,4.5057403089581055044432714744588826e-1L,2.2696723314583158080340978057276064e-1L,-2.6967233145831580803409780572760635e-2L},
                typename type::vector_t{8.3333333333333333333333333333333333e-2L,4.1666666666666666666666666666666667e-1L,4.1666666666666666666666666666666667e-1L,8.3333333333333333333333333333333333e-2L}};
        else if constexpr (method_order==5)
            return typename type::matrix_t{
                typename type::vector_t{0.0L,0.0L,0.0L,0.0L,0.0L},
                typename type::vector_t{6.7728432186156897969267419174073483e-2L,1.1974476934341168251615379970493965e-1L,-2.1735721866558113665511351745074293e-2L,1.0635824225415491883105056997129926e-2L,-3.7001392424145306021611522544979462e-3L},
                typename type::vector_t{4.0625e-2L,3.0318418332304277801796699838244475e-1L,1.7777777777777777777777777777777778e-1L,-3.0961961100820555795744776160222533e-2L,9.375e-3L},
                typename type::vector_t{5.3700139242414530602161152254497947e-2L,2.615863979968067303391171652250923e-1L,3.7729127742211366922106690730062985e-1L,1.5247745287881053970606842251728257e-1L,-1.7728432186156897969267419174073482e-2L},
                typename type::vector_t{5.0e-2L,2.7222222222222222222222222222222222e-1L,3.5555555555555555555555555555555555e-1L,2.7222222222222222222222222222222222e-1L,5.0e-2L}};
        else if constexpr (method_order==6)
            return typename type::matrix_t{
                typename type::vector_t{0.0L,0.0L,0.0L,0.0L,0.0L,0.0L},
                typename type::vector_t{4.5679805133755038575653446730922972e-2L,8.1867817008970666864969819153683861e-2L,-1.487460578908983676559396150555719e-2L,7.6276761182509598020429585787853333e-3L,-4.4717804405737092705509645159422079e-3L,1.6434260039545343679772145784381566e-3L},
                typename type::vector_t{2.5908385387879822499353401700632774e-2L,2.1384080863282571965204027280905642e-1L,1.3396073565086083664894428137964585e-1L,-2.4004074733154873937276256033001734e-2L,1.1807696377659694346907243344183639e-2L,-4.1293095563937473670444402209564835e-3L},
                typename type::vector_t{3.746264288972708070037777355428982e-2L,1.7742978177126379581139916076182869e-1L,3.0143326325089805044563651859365607e-1L,1.4346845286688233985941598118100849e-1L,-2.460333048390222949373386870304409e-2L,7.4249479454535108339799316327005587e-3L},
                typename type::vector_t{3.1689907329378798965356118754895177e-2L,1.9370925858949719942885736862195453e-1L,2.6980151239949221670631730398186901e-1L,2.9230379430683301327395422406621153e-1L,1.0736966113995282329333658495232847e-1L,-1.2346471800421705242320113397589639e-2L},
                typename type::vector_t{3.3333333333333333333333333333333333e-2L,1.8923747814892349015830640410601233e-1L,2.7742918851774317650836026256065434e-1L,2.7742918851774317650836026256065434e-1L,1.8923747814892349015830640410601233e-1L,3.3333333333333333333333333333333333e-2L}};
        else if constexpr (method_order==7)
            return typename type::matrix_t{
                typename type::vector_t{0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L},
                typename type::vector_t{3.2846264328292647881547377380331421e-2L,5.9322894027551404504198527596713417e-2L,-1.0768594451189267105573388916275431e-2L,5.5975917805697772306731264515785011e-3L,-3.488929970807462766480046045821455e-3L,2.2170965889145396997313080577989019e-3L,-8.3827044261510438011301150805792508e-4L},
                typename type::vector_t{1.8002223201815165703973460160284566e-2L,1.5770113064168904204778906247437054e-1L,1.0235481204686191521394666173809335e-1L,-1.8478259273459043983140512247516131e-2L,9.5775801007414059542896765174719021e-3L,-5.6818645662243775729731799680753607e-3L,2.0999811132187857342288903709879644e-3L},
                typename type::vector_t{2.7529761904761904761904761904761902e-2L,1.2778825555983746958666847435455748e-1L,2.3748565272164544351033623662864482e-1L,1.2190476190476190476190476190476191e-1L,-2.1612962116714131801400725487963681e-2L,1.0624768120945504418681728790475672e-2L,-3.7202380952380952380952380952380952e-3L},
                typename type::vector_t{2.1709542696305023789580633438535843e-2L,1.4409488824700735157832338311310851e-1L,2.0629511050418990575464583462320924e-1L,2.6228778308298285350695003605703994e-1L,1.1351787855806939649498884940258779e-1L,-1.9288106960906068042438859329337396e-2L,5.8073006077086438198360636492392419e-3L},
                typename type::vector_t{2.4647794252138913903922535317581733e-2L,1.3619592709186843430561889508723425e-1L,2.1936162057573877447541555718650259e-1L,2.3821193202895403229313639735794531e-1L,2.2664128505612057881450890005695657e-1L,7.909012965323156950115167554831973e-2L,-9.0367405187688383577378535708076123e-3L},
                typename type::vector_t{2.3809523809523809523809523809523809e-2L,1.3841302368078297400535020314503315e-1L,2.1587269060493131170893551114068114e-1L,2.4380952380952380952380952380952381e-1L,2.1587269060493131170893551114068114e-1L,1.3841302368078297400535020314503315e-1L,2.3809523809523809523809523809523809e-2L}};
        else if constexpr (method_order==8)
            return typename type::matrix_t{
                typename type::vector_t{0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L},
                typename type::vector_t{2.4737514438875539852864281100249016e-2L,4.4892662602755022208806536433040278e-2L,-8.1407677423897545665601763531105896e-3L,4.2540825485459384954106613515856382e-3L,-2.7042050750158204059997790367644178e-3L,1.845386694481988769056269152924056e-3L,-1.2262209861064896557130967227224978e-3L,4.7147326405026763341242346446679832e-4L},
                typename type::vector_t{1.3258719822130254949486806628934897e-2L,1.206497328240974033374613514561339e-1L,7.9997635786651606350484180603740961e-2L,-1.4480830962625474675528977016578683e-2L,7.6319492068340046025017104515804291e-3L,-4.8326689231347168793347168846303201e-3L,3.1049500159208759803192353768410752e-3L,-1.179578486445104737644956314998851e-3L},
                typename type::vector_t{2.1034356725214087176749624414752258e-2L,9.6280178312081402576316895975299821e-2L,1.8904164984432690362365580574054637e-1L,1.0124595564768955192061451731175456e-1L,-1.8137320211454220198518540068934206e-2L,9.4170098024305847560818239159936158e-3L,-5.6088616425374647708178885822648324e-3L,2.0774225710097205315891311201767933e-3L},
                typename type::vector_t{1.5779720286133136611268011736965683e-2L,1.1096097521429048446231392147014384e-1L,1.6114433643932159762603851463787957e-1L,2.2436671754080616098204502577003992e-1L,1.0498344168166238886291196838934963e-1L,-1.8480303602574721241535467186671696e-2L,9.0719352596716171151791369125779021e-3L,-3.1772138680712300338924815576092221e-3L},
                typename type::vector_t{1.9036721343587961880502099172142061e-2L,1.0224716355583214371117679751103627e-1L,1.7539401516488689926145505543850536e-1L,1.9859744812251793618102477524952352e-1L,2.2071022829197741545905546271768443e-1L,9.0563710455100576031636157950132416e-2L,-1.529761925234438364596531856825521e-2L,4.5984230350126021933703362282077491e-3L},
                typename type::vector_t{1.7385669593092589509444719392675173e-2L,1.0657833455785950934720912961060256e-1L,1.6871595954727019361306406940094817e-1L,2.0893360240436776118952626473787079e-1L,2.0197531478080600228811582434951816e-1L,1.7870211398414193694868051490698548e-1L,6.0459450968997997482689496454837326e-2L,-6.8803715817326827100071382431059395e-3L},
                typename type::vector_t{1.7857142857142857142857142857142631e-2L,1.0535211357175301969149603288787757e-1L,1.7056134624175218238212033855387492e-1L,2.0622939732935194078352648570110395e-1L,2.0622939732935194078352648570110581e-1L,1.705613462417521823821203385538733e-1L,1.0535211357175301969149603288787874e-1L,1.7857142857142857142857142857142631e-2L}};
        else if constexpr (method_order==9)
            return typename type::matrix_t{
                typename type::vector_t{0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L,0.0L},
                typename type::vector_t{1.9293838201043212559591809836576235e-2L,3.5125520977621796835804754698841359e-2L,-6.3641024187047848382867337621813776e-3L,3.3337771969983822838953024622409241e-3L,-2.1368470176082402286362169591461366e-3L,1.49429916272822666158617039699983e-3L,-1.0740606993816592835517151742116807e-3L,7.3377266653928893305941646073068058e-4L,-2.8519577496630157963541016901881371e-4L},
                typename type::vector_t{1.0184080408227649523905913084940606e-2L,9.5086508449362110322827043281952955e-2L,6.3931995628438088310135533623375416e-2L,-1.1585731333842718614601805944510917e-2L,6.1499369005518014276688605000503924e-3L,-3.9808257871826992002764583837745418e-3L,2.7553240894191530283967676671437299e-3L,-1.8475051393836076788723061722057375e-3L,7.1307702904134615787373879735687077e-4L},
                typename type::vector_t{1.6569369843571800049707628059239573e-2L,7.5093517096206969587369572240439695e-2L,1.528820610248834465326923641513242e-1L,8.4085478688913126843291549224452553e-2L,-1.5130673172334494698279439134431745e-2L,8.0006035703992085819150559997474664e-3L,-5.0963258617414847502622958857879155e-3L,3.2896245700396027380943035100630938e-3L,-1.2523876730272542399047725193765277e-3L},
                typename type::vector_t{1.1990017361111111111111111111110974e-2L,8.7869833723063135827338822896489924e-2L,1.2868222230658306372694813717586041e-1L,1.897580803159183818968257919166025e-1L,9.2879818594104308390022675736961138e-2L,-1.6543824829395209339260025846742735e-2L,8.5871339434978039134046721138253774e-3L,-5.1221529426603733041689628818852923e-3L,1.8988715277777777777777777777777104e-3L},
                typename type::vector_t{1.5141276561916143128793661408266021e-2L,7.9458056210363159785075556504539581e-2L,1.4236568211182235239061510517547616e-1L,1.6521365191612396397565071007010965e-1L,2.0089031036054311147832479060835661e-1L,8.9128776797610045714274216845404816e-2L,-1.5612704774802578892339554861636348e-2L,7.6541636841957929358002877741633534e-3L,-2.6804809546829111608187391703502422e-3L},
                typename type::vector_t{1.3175811859847542731015150091531629e-2L,8.459518591978637020204216618681075e-2L,1.3451403216066171461195604162254174e-1L,1.7719508127370587175784222445363445e-1L,1.7960970028765681535237649097387191e-1L,1.8479998682036589117216757201437055e-1L,7.3337360621642779330217275666310528e-2L,-1.2338827668959347799657183267348473e-2L,3.704808480661239364982975803948145e-3L},
                typename type::vector_t{1.417408466385519046852429905790681e-2L,8.2013908113863473590110443553875801e-2L,1.3834341694946252692390452446389464e-1L,1.7171995632379494589597959567286352e-1L,1.8789648420581685700868156843306465e-1L,1.6988047828952479027367046360762227e-1L,1.4363345866878565247863954305186433e-1L,4.7622159802780965687365105315765329e-2L,-5.4049493121543236707029209476883653e-3L},
                typename type::vector_t{1.3888888888888888888888888888889388e-2L,8.2747680780402762523169860014603253e-2L,1.3726935625008086764035280928968679e-1L,1.7321425548652317255756576606985943e-1L,1.8575963718820861678004535147392211e-1L,1.7321425548652317255756576606986015e-1L,1.3726935625008086764035280928968537e-1L,8.2747680780402762523169860014604949e-2L,1.3888888888888888888888888888888561e-2L}};
    }
    using type::state_vector_t;
    constexpr static const typename type::vector_t a=make_a();
    constexpr static const typename type::vector_t sigma=make_sigma();
    constexpr static const typename type::matrix_t b=make_b();
    template<class RHS>
    auto do_step(const typename type::state_vector_t& y,RHS rhs,const Number_t& t,const Number_t& dt,const Number_t& epsilon)
    {
        return this->do_step_impl(a,b,sigma,y,rhs,t,dt,epsilon);
    }
};


/**
Класс Lobatto3b, определяющий наборы коэффициентов 'a', 'b', 'sigma' для LobattoIIIb
'a' - узлы сетки на обезразмеренном отрезке времени,
'b' - коэффициенты для вычисления промежуточных векторов 'ks' (аналог узлов сетки, только в пространстве, не на отрезке времени),
'sigma' - веса квадратурной формулы для соответствующих узлов.
**/
template<class Number_t,std::size_t system_order,std::size_t method_order>
class Lobatto3b: protected Runge_Kutta_implicit<Number_t,system_order,method_order>
{
    using type = Runge_Kutta_implicit<Number_t,system_order,method_order>;
public:
    constexpr static const typename type::matrix_t make_b()
    {
        if constexpr (method_order==2)
            return typename type::matrix_t{
                typename type::vector_t{0.5L,0.0L},
                typename type::vector_t{0.5L,0.0L}};
        else if constexpr (method_order==3)
            return typename type::matrix_t{
                typename type::vector_t{1.0L/6.0L,-1.0L/6.0L,0.0L},
                typename type::vector_t{1.0L/6.0L,1.0L/3.0L,0.0L},
                typename type::vector_t{1.0L/6.0L,5.0L/6.0L,0.0L}};
        else if constexpr (method_order==4)
            return typename type::matrix_t{
                typename type::vector_t{8.3333333333333333333333333333333333e-2L,-1.3483616572915790401704890286380318e-1L,5.1502832395824570683715569530469843e-2L,0.0L},
                typename type::vector_t{8.3333333333333333333333333333333333e-2L,2.2696723314583158080340978057276064e-1L,-3.3907364229143883777660480779221592e-2L,0.0L},
                typename type::vector_t{8.3333333333333333333333333333333333e-2L,4.5057403089581055044432714744588826e-1L,1.8969943352083508586325688609390603e-1L,0.0L},
                typename type::vector_t{8.3333333333333333333333333333333333e-2L,3.6516383427084209598295109713619682e-1L,5.5150283239582457068371556953046984e-1L,0.0L}};
        else if constexpr (method_order==5)
            return typename type::matrix_t{
                typename type::vector_t{5.0e-2L,-9.6521464124632000054900393281066737e-2L,6.6666666666666666666666666666666667e-2L,-2.0145202542034666611766273385599929e-2L,0.0L},
                typename type::vector_t{5.0e-2L,1.5247745287881053970606842251728257e-1L,-4.0440112458214603488319707637841675e-2L,1.0635824225415491883105056997129926e-2L,0.0L},
                typename type::vector_t{5.0e-2L,2.8886363427630577799737935090204473e-1L,1.7777777777777777777777777777777778e-1L,-1.6641412054083555775157128679822506e-2L,0.0L},
                typename type::vector_t{5.0e-2L,2.615863979968067303391171652250923e-1L,3.9599566801377015904387526319339723e-1L,1.1974476934341168251615379970493965e-1L,0.0L},
                typename type::vector_t{5.0e-2L,2.9236742476425688883398849560782215e-1L,2.8888888888888888888888888888888889e-1L,3.6874368634685422227712261550328896e-1L,0.0L}};
        else if constexpr (method_order==6)
            return typename type::matrix_t{
                typename type::vector_t{3.3333333333333333333333333333333333e-2L,-7.0092455626458075245278394957775769e-2L,6.1796918498809558113257375632222398e-2L,-3.4367729981066381604897113071568057e-2L,9.3299337753815654035847990637880952e-3L,0.0L},
                typename type::vector_t{3.3333333333333333333333333333333333e-2L,1.0736966113995282329333658495232847e-1L,-3.6069398502611984526326140173523631e-2L,1.7310522505167190744705699424134963e-2L,-4.4717804405737092705509645159422081e-3L,0.0L},
                typename type::vector_t{3.3333333333333333333333333333333333e-2L,1.9938360914193799858198254230684248e-1L,1.4346845286688233985941598118100849e-1L,-2.4004074733154873937276256033001733e-2L,5.2029211506786540054689021913778951e-3L,0.0L},
                typename type::vector_t{3.3333333333333333333333333333333333e-2L,1.8403455699824483615283750191463443e-1L,3.0143326325089805044563651859365607e-1L,1.3396073565086083664894428137964585e-1L,-1.0146130993014508423676138200830156e-2L,0.0L},
                typename type::vector_t{3.3333333333333333333333333333333333e-2L,1.9370925858949719942885736862195454e-1L,2.6011866601257598576365456313651938e-1L,3.1349858702035516103468640273417797e-1L,8.1867817008970666864969819153683861e-2L,0.0L},
                typename type::vector_t{3.3333333333333333333333333333333333e-2L,1.7990754437354192475472160504222423e-1L,3.117969184988095581132573756322224e-1L,2.1563227001893361839510288692843194e-1L,2.593299337753815654035847990637881e-1L,0.0L}};
        else if constexpr (method_order==7)
            return typename type::matrix_t{
                typename type::vector_t{2.3809523809523809523809523809523809e-2L,-5.2533708335700574120948015395724229e-2L,5.2652779508184141373465335535682369e-2L,-3.8095238095238095238095238095238095e-2L,1.9039800071463282903331944936881433e-2L,-4.8731569582325644415635507911252871e-3L,0.0L},
                typename type::vector_t{2.3809523809523809523809523809523809e-2L,7.909012965323156950115167554831973e-2L,-3.0082252634905692192822699037316039e-2L,1.8715143902415751594872873751978001e-2L,-8.8615894584634430627587891140369716e-3L,2.2170965889145396997313080577989018e-3L,0.0L},
                typename type::vector_t{2.3809523809523809523809523809523809e-2L,1.4531761969338837508155554931940708e-1L,1.1351787855806939649498884940258779e-1L,-2.4409970464642816920111989492718176e-2L,9.5775801007414059542896765174719013e-3L,-2.2370284324372770364175505106555656e-3L,0.0L},
                typename type::vector_t{2.3809523809523809523809523809523809e-2L,1.3523521671256684609907229694902606e-1L,2.3223362468465780937437789434758218e-1L,1.219047619047619047619047619047619e-1L,-1.6360934079726497665442383206901035e-2L,3.177806968216127906277906196007087e-3L,0.0L},
                typename type::vector_t{2.3809523809523809523809523809523809e-2L,1.4065005211322025104176775365568871e-1L,2.0629511050418990575464583462320924e-1L,2.6821949427416662644392151330224198e-1L,1.0235481204686191521394666173809335e-1L,-6.9045960126054010762053461743739305e-3L,0.0L},
                typename type::vector_t{2.3809523809523809523809523809523809e-2L,1.3619592709186843430561889508723425e-1L,2.2473428006339475477169430025471811e-1L,2.2509437990710805792893665005754581e-1L,2.4595494323983700390175821017799718e-1L,5.9322894027551404504198527596713417e-2L,0.0L},
                typename type::vector_t{2.3809523809523809523809523809523809e-2L,1.4328618063901553844691375393615843e-1L,1.9683289053346802880560356620379971e-1L,2.819047619047619047619047619047619e-1L,1.6321991109674717033547017560499877e-1L,1.9094673201648354812629821854075738e-1L,0.0L}};
        else if constexpr (method_order==8)
            return typename type::matrix_t{
                typename type::vector_t{1.7857142857142857142857142857142631e-2L,-4.0592254544495558974329932166468937e-2L,4.3921540512686668232481179022219043e-2L,-3.6693154467132177093605731558013309e-2L,2.3991833869792740109631243959109918e-2L,-1.1266667700169587582956970426641225e-2L,2.7815594721750581659230683126518787e-3L,0.0L},
                typename type::vector_t{1.7857142857142857142857142857142631e-2L,6.0459450968997997482689496454838378e-2L,-2.4766304590524872846824351738897966e-2L,1.7758540173365830172565864838734536e-2L,-1.097948695121725330870835577443655e-2L,5.0268042735386233444104194750097963e-3L,-1.2262209861064896557130967227225442e-3L,0.0L},
                typename type::vector_t{1.7857142857142857142857142857142631e-2L,1.1038049254890935182366073474172713e-1L,9.0563710455100576031636157950132504e-2L,-2.2344933118788238926140722175395304e-2L,1.1386309377783603689753738926155842e-2L,-4.8326689231347168793347168846302215e-3L,1.1398560864154160453122988858908222e-3L,0.0L},
                typename type::vector_t{1.7857142857142857142857142857142631e-2L,1.0317891909219839818588095699836862e-1L,1.8253766996508743356716923571366439e-1L,1.0498344168166238886291196838935069e-1L,-1.8137320211454220198518540068934606e-2L,6.3119785443945133407638687915642126e-3L,-1.3814408802708052853932628538315608e-3L,0.0L},
                typename type::vector_t{1.7857142857142857142857142857142631e-2L,1.0673355445202382497688929574171017e-1L,1.6424936769735766904135646976230958e-1L,2.243667175408061609820450257700395e-1L,1.012459556476895519206145173117545e-1L,-1.1976323723335251185048897159790793e-2L,2.1731944795546215056150758895100477e-3L,0.0L},
                typename type::vector_t{1.7857142857142857142857142857142631e-2L,1.0421225748533760364618373400198783e-1L,1.7539401516488689926145505543850387e-1L,1.9484308795156833709377274677494924e-1L,2.2857433044814017970966720787650033e-1L,7.9997635786651606350484180603741184e-2L,-5.028378977156332132164701853848496e-3L,0.0L},
                typename type::vector_t{1.7857142857142857142857142857142631e-2L,1.0657833455785950934720912961060104e-1L,1.655345419682135590377099190788643e-1L,2.1720888428056919409223484147554102e-1L,1.8847085715598611061096062086237102e-1L,1.9532765083227705522894469029277136e-1L,4.4892662602755022208806536433040335e-2L,0.0L},
                typename type::vector_t{1.7857142857142857142857142857142631e-2L,1.025705540995779615255729645752267e-1L,1.818280139419217699650773089805151e-1L,1.8223756345955920067389524174199486e-1L,2.4292255179648411787713221725911861e-1L,1.2663980572906551414963915953165449e-1L,1.4594436811624857866582596505434762e-1L,0.0L}};
        else if constexpr (method_order==9)
            return typename type::matrix_t{
                typename type::vector_t{1.3888888888888888888888888888888561e-2L,-3.2201785462781049136033796318078201e-2L,3.6616080612255012077627666284754148e-2L,-3.3429420929606782752350330725495349e-2L,2.5396825396825396825396825396824787e-2L,-1.5619060682201884450601052753486368e-2L,7.0476209807922580738656498350504659e-3L,-1.6991488041718395267938506084580441e-3L,0.0L},
                typename type::vector_t{1.3888888888888888888888888888888561e-2L,4.7622159802780965687365105315763513e-2L,-2.0468766194107845784129194852202531e-2L,1.6022325356144010036764002880109237e-2L,-1.1498682057037766338016808097283472e-2L,6.8861128838372528815426198624926126e-3L,-3.0648090527748729616466526676676413e-3L,7.3377266653928893305941646073073998e-4L,0.0L},
                typename type::vector_t{1.3888888888888888888888888888888561e-2L,8.6584041129010884020032117616147375e-2L,7.3337360621642779330217275666309971e-2L,-1.9700995965708996038950828511862417e-2L,1.1620531554942463292097179597816928e-2L,-6.4308328819582765688938302454615701e-3L,2.7553240894191530283967676671436319e-3L,-6.4745719160577267473028422465370608e-4L,0.0L},
                typename type::vector_t{1.3888888888888888888888888888888561e-2L,8.1155073229047548984911994365987006e-2L,1.464508516034279145424571288938545e-1L,8.9128776797610045714274216845407126e-2L,-1.774204374450483878223640672478506e-2L,8.0006035703992085819150559997479908e-3L,-3.1547368408878239492149631474275928e-3L,7.1385458292997666362805052399786077e-4L,0.0L},
                typename type::vector_t{1.3888888888888888888888888888888561e-2L,8.3699551367393907951075549069003135e-2L,1.3272478514682211212589997536234294e-1L,1.8732307012584349345713350781886726e-1L,9.2879818594104308390022675736960873e-2L,-1.4108814639320320899567741749006992e-2L,4.5445711032587555144528339273424215e-3L,-9.5187058699114542790568905439819159e-4L,0.0L},
                typename type::vector_t{1.3888888888888888888888888888888561e-2L,8.2033826197472785859541809490607058e-2L,1.4042409309096869158956777243711303e-1L,1.6521365191612396397565071007011216e-1L,2.035016809327134555622817581987069e-1L,8.4085478688913126843291549224453077e-2L,-9.181495353347046902104319604169114e-3L,1.5926075513552135382578656486179302e-3L,0.0L},
                typename type::vector_t{1.3888888888888888888888888888888561e-2L,8.339513797200853519790014423925861e-2L,1.3451403216066171461195604162254186e-1L,1.7964508836848144912645959631532162e-1L,1.7413910563326615348794817187610505e-1L,1.9291525145223216859651659458172252e-1L,6.3931995628438088310135533623375463e-2L,-3.8363603486081214968622576015424479e-3L,0.0L},
                typename type::vector_t{1.3888888888888888888888888888888561e-2L,8.201390811386347359011044355387433e-2L,1.403341653028557406019994619573525e-1L,1.6632814260268591967602314620736849e-1L,1.9725831924524638311806215957120439e-1L,1.5719193013037916252080176318975157e-1L,1.5773812244418871342448200414188764e-1L,3.512552097762179683580475469884149e-2L,0.0L},
                typename type::vector_t{1.3888888888888888888888888888888561e-2L,8.4446829584574602049963710623063254e-2L,1.30221735269288609566487159454634e-1L,1.8883331616872505700816681882334804e-1L,1.6036281179138321995464852607709553e-1L,2.0664367641612995530991609679535665e-1L,1.0065327563782585556272514300493069e-1L,1.1494946624318381165920365633268328e-1L,0.0L}};
    }
    using type::state_vector_t;
    constexpr static const typename type::vector_t a=Lobatto3a<Number_t,system_order,method_order>::make_a();
    constexpr static const typename type::vector_t sigma=Lobatto3a<Number_t,system_order,method_order>::make_sigma();
    constexpr static const typename type::matrix_t b=make_b();
    template<class RHS>
    auto do_step(const typename type::state_vector_t& y,RHS rhs,const Number_t& t,const Number_t& dt,const Number_t& epsilon)
    {
        return this->do_step_impl(a,b,sigma,y,rhs,t,dt,epsilon);
    }
};


/**
Класс Lobatto3_pair, определяющий наборы коэффициентов 'a', 'b', 'sigma' для пары LobattoIIIa-LobattoIIIb
'a' - узлы сетки на обезразмеренном отрезке времени,
'b' - коэффициенты для вычисления промежуточных векторов 'ks' (аналог узлов сетки, только в пространстве, не на отрезке времени),
'sigma' - веса квадратурной формулы для соответствующих узлов.

Все числа Number_t становятся парами qp_pair.
**/
template<class Number_t,std::size_t system_order,std::size_t method_order>
class Lobatto3_pair: protected Runge_Kutta_implicit<qp_pair<Number_t>,system_order,method_order>
{
    using type = Runge_Kutta_implicit<qp_pair<Number_t>,system_order,method_order>;
public:
    constexpr static const typename type::vector_t make_a()
    {
        auto a1=Lobatto3a<Number_t,system_order,method_order>::make_a();
        typename type::vector_t a{};
        for(std::size_t i=0;i<method_order;i++) a[i]=qp_pair<Number_t>{a1[i],a1[i]};
        return a;
    }
    constexpr static const typename type::vector_t make_sigma()
    {
        auto sigma1=Lobatto3a<Number_t,system_order,method_order>::make_sigma();
        typename type::vector_t sigma{};
        for(std::size_t i=0;i<method_order;i++) sigma[i]=qp_pair<Number_t>{sigma1[i],sigma1[i]};
        return sigma;
    }
    constexpr static const typename type::matrix_t make_b()
    {
        auto b1=Lobatto3a<Number_t,system_order,method_order>::make_b();
        auto b2=Lobatto3a<Number_t,system_order,method_order>::make_b();
        typename type::matrix_t b{};
        for(std::size_t i=0;i<method_order;i++) for(std::size_t o=0;o<method_order;o++) b[i][o]=qp_pair<Number_t>{b1[i][o],b2[i][o]};
        return b;
    }
    using type::state_vector_t;
    const typename type::vector_t a=make_a();
    const typename type::vector_t sigma=make_sigma();
    const typename type::matrix_t b=make_b();
    template<class RHS>
    auto do_step(const typename type::state_vector_t& y,RHS rhs,const Number_t& t,const Number_t& dt,const Number_t& epsilon)
    {
        return this->do_step_impl(a,b,sigma,y,[&rhs](const typename type::state_vector_t& x,const qp_pair<Number_t>& t){return rhs(x,t.first);},qp_pair<Number_t>{t,t},qp_pair<Number_t>{dt,dt},qp_pair<Number_t>{epsilon,epsilon});
    }
};
}

template<class Number_t>
constexpr rg::qp_pair<Number_t> std::abs(const rg::qp_pair<Number_t>& qppair)
{
    return rg::qp_pair<Number_t>{std::abs(qppair.first),std::abs(qppair.second)};
}

template<class Number_t>
std::ostream& operator<<(std::ostream& s,const rg::qp_pair<Number_t>& qppair)
{
    return s<<"<"<<qppair.first<<","<<qppair.second<<">";
}
#endif