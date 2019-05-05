#include <cstddef>
#include<iostream>

#include <boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include<boost/shared_ptr.hpp>

using namespace boost::numeric::ublas;
using namespace boost;

/*В данном hpp файле реализован класс массивов Array со следующими особенностями.
1)Array<T> A(B); Происходит копирование Умного Указателя на память, хранящую элементы B в A.
Т.е. изменение элементов А изменит элементы B и наоборот.
Это полезно при возвращении больших массивов из функции, если массив создан внутри функции.
Т.к. копирование указателя, а не элементов приведёт к серьезной экономии времени, а умный указатель
предотвратит удаление памяти, на которую указывал массив, созданный в функции.
Конечно, в обычных функциях можно подать ссылку на массив std::vector<T>& в качестве аргумента, но что
делать, если хочется переопределить оператор?
template<typename T>
std::vector<T> operator+(std::vector<T> a, std::vector<T> b){
    std::vector<T> c(a.size());
    for(unsigned i=0;i<a.size();++i)
        c[i]=a[i]+b[i];
    return c;
}
Такой код будет работать за 3*N операций(N сложений + 2*N копирований), не считая выделения памяти.
template<typename T>
Array<T> operator+(Array<T> left, Array<T> right){
    Array<T> result(left.size());
    for(unsigned i=0;i<result.size();++i)
        result[i]=left[i]+right[i];
    return result;
}
Такой код будет работать за 2*N операций, не считая выделения памяти.
Также здесь я не учитываю разные операторы [] у Array и std::vector.
Нужно проверить верность данного утверждения на практике.
Есть также boost::numeric::ublas::vector, но у него почему-то возникают проблемы в случае
vector<vector<double> > и также приходится прописывать перегрузку:
template<typename T>
vector<T> operator+(vector<T> a,vector<T> b){
    vector<T> c(a.size());
    for(unsigned i=0;i<a.size();++i)
        c(i)=a(i)+b(i);
    return c;
}
и template<typename T, typename U>
vector<T> operator*(U b, vector<T> a){...}

При использовании оператора копирования следует также помнить, что память выделенная под массив затрётся
только после окончания времени существования последней копии этого массива(не поэлементной копии).
2)Array<T> A(B.copy()); Копирование поэлементное.
3)A=B; Также поэлементное копирование и затирание памяти из A(при условии, что это последняя копия).
4)A.subjoin(B); Обычное копирование(умного указателя) массива B в следующую часть массива A.
4.1)Array<Array<T> > A;
  Array<T> B;
    -//- Массивы A и B как-то заполнены.
    B.subjoin(A[0]); Не поэлементное копирование.
    В таких более сложных случаях можно забыть, что копирование не поэлементное и массив B остаётся связан
с массивом A[0].
5)A=Array<T>(B); Здесь копирование поэлементное, т.к. используется оператор присваивания.
    Если при задании переменной A не сделать её непоэлементной копией массива B,
    то потом сделать это с помощью конструктора копирования уже нельзя.
  Для этого используется:
  A.mirror(B);
6)Array<T> A(M),B;
B.mirror(A);    Указатель на Arr<T> в B совпадает с указателем в A. Т.е. если менять элементы
массива B, будут меняться элементы массива A,и наоборот.
B=Array<T>(N); Здесь указатель на Arr<T> в B заменяется на новый, т.е. элементы старого не меняются,
а значит не меняются и элементы A.
7)A.resize(N); A.size() становится равным N. Если размер был больше N, то удаляются части массива с
конца так, чтобы размер массива стал меньше N и вызывается resize.
Если размер массива был меньше N, то добавялется новая часть массива по размеру равная N-A.size().
Т.о. сохранить все элементы массива после этой операции можно только в случае N>=A.size(). 
*/
template<typename T>
struct Array;

template<typename U>
struct Arr
{
    // Список операций:
    //
    explicit Arr(size_t size, const U& value = U())
    //   конструктор класса, который создает
    //   Arr размера size, заполненный значениями
    //   value типа U. Считаем, что у типа U есть
    //   конструктор, который можно вызвать без
    //   без параметров, либо он ему не нужен.
    : size_(size)    {
        data_= static_cast<U*>(operator new[] (size * sizeof(U)));
        for(unsigned i=0; i< size_; i++)
            new (data_+i)U(value);
    }
    Arr(const Arr<U> &A)
    //   Конструктор копирования, который копирует поэлементно первую часть массива А и все последующие,
    //   если они есть, в первую и последующие части нашего массива соответственно.
    :    size_(A.size_) {
        data_= static_cast<U*>(operator new[] (size_ * sizeof(U)));
        for(unsigned i=0; i< size_; i++)
            new (data_+i)U(A.data_[i]);
        data_next=A.data_next;
    }
    Arr(){
        data_=0;
        size_=0;
    }   
    ~Arr(){    
    //   деструктор
        for (int i = 0; i < size_; i++)
            data_[i].~U();
        operator delete [] (data_);    
    }    
    Arr& operator=(const Arr<U> &A){
    //   Оператор присваивания.
    //  Удаляем содержимое нашего массива, выделяем новую память и заполняем значениями из массива А.
    //  Тоже самое делаем для следующей части массива, используя конструктур копирования
    //   и функцию copy, определённые в классе Array. 
        if(this->data_!=A.data_){
            this->~Arr();
            this->size_=A.size_;
            data_= static_cast<U*>(operator new[] (size_ * sizeof(U)));
            for(unsigned i=0; i< size_; i++)
                new (data_+i)U(A.data_[i]);     
        }
        data_next=A.data_next;
        return *this;
    }    
    void subjoin(size_t size, const U& value = U()){
    //  Добавляем новую часть массива размера size, заполненную значениями value.
        if(data_next.array==0)
            data_next.mirror(Array<U>(size, value));
        else
            data_next.subjoin(size, value);
    }
    void subjoin(const Array<U>& A){
    //  Добавляем новую часть массива и копируем в неё Array<U> A,
    //  стоит обратить внимание на то, что копирование происходит без использования функции copy,
    // а значит копируется лишь указатель, и при изменении элементов массива A будут меняться
    // и элементы нашего массива, и наоборот. Чтобы избежать подобного необходимо использовать
    //  функцию copy при передаче аргумента.
        if(data_next.array==0)
            data_next.mirror(A);
        else
            data_next.subjoin(A);
    }
    size_t size() const{
    //  Возвращает полный размер нашего массива, т.е. сумму размеров всех его частей.
        if(data_next.array==0)    
            return this->size_;
        else
            return size_+data_next.size();
    }
    U& operator[](size_t i){
    //  Возвращает элемент массива по номеру, игнорируя разбиение на части.
        if(i<size_)
            return data_[i];
        else
            return data_next[i-size_];
    }
    const U& operator[](size_t i) const{
    //  Возвращает элемент массива по номеру, игнорируя разбиение на части.
        if(i<size_)
            return data_[i];
        else
            return data_next[i-size_];
    }
    
    U *data_;   //  Указатель на память, в которой хранится первая часть данного массива.
    Array<U> data_next;    // Вторая часть данного массива.
    size_t size_;   // Размер первой части данного массива.
};


// Обёртка над Arr, позволяющая копировать только указатели на память, а не сами элементы массива,
//  без явного использования указателей.
template<typename T>
struct Array{
    T& operator[](size_t i){return (*array)[i];}
    const T& operator[](size_t i) const{return (*array)[i];}
    size_t size() const{
        if(array!=0)
            return (this->array)->size();
        else
            return 0;
    }
    explicit Array(size_t size, const T& value = T()) : array(new Arr<T>(size,value)){}
    Array(Arr<T>* copy) : array(copy){ } // Конструктор класса, получающий на вход обычный указатель на Arr,
    // и конвертирующий его в умный, записывая в array. Очевидно, изменение элементов copy приведет к
    // изменению элементов текущего массива и наоборот.
    Array(const Array<T>& copy) : array(copy.array){ } //Конструктор, копирующий умный указатель на Arr
    // в array. Очевидно, изменение элементов copy приведет к
    // изменению элементов текущего массива и наоборот.    
    Array() : array(0) {} //Конструктор, создающий пустой указатель.

    Array& operator=(const Array<T> &A){
    // Оператор присваивания.
    // Копирование массива A в текущий с использованием констуктора Arr, т.е. поэлементно.
    // Освобождать память не нужно, умный указатель справится с этим сам.
        if(A.array==0)
            this->array=0;
        else
            this->array=shared_ptr<Arr<T> >(new Arr<T>(*(A.array)));
        return *this;
    }
    Array<T> copy() const{
    // Копирование массива A в возвращаемый с использованием констуктора Arr, т.е. поэлементно.
    // Освобождать память не нужно, умный указатель справится с этим сам.
        if(this->array==0)
            return Array<T>();
        else
            return Array<T>(new Arr<T>(*(this->array))); 
    }
    void mirror(const Array<T>& A){
    //  Копирование указателей. Очевидно, изменение элементов А приведет к
    // изменению элементов текущего массива и наоборот.
        this->array=A.array;
    }
    void subjoin(size_t size, const T& value = T()){
        if(array!=0)
            array->subjoin(size, value);
        else  
            this->mirror(Array<T>(size, value));
    }
    void subjoin(const Array<T>& A){
        if(array!=0)
            array->subjoin(A);
        else  
            this->mirror(A); 
    }
    void resize(size_t n){
        if(n==this->size())
                return;
        if(n>this->size()){
            this->subjoin(n-this->size());
        }else{
            Array<T> last(*this),penult;
            if(this->array->data_next.array==0){
                this->array=0;
            }else{    
                while(last.array->data_next.array!=0){
                    penult.mirror(last);
                    last.mirror(last.array->data_next);
                    if(this->size()-last.array->data_next.size()>=n){
                        penult.array->data_next.array=0;
                        break;
                    }       
                }    
                if(last.array->data_next.array==0)
                    penult.array->data_next.array=0;
            }    
            this->resize(n);
        }
    }

    shared_ptr<Arr<T> > array;  //  Умный указатель на массив Arr.
};

//  В случае массива массива конструктор класса работает некорректно. А именно создаёт массив
//  указателей на один и тот же массив. Проблема в том, что подавая в value массив Array,
//  мы будем использовать конструктор копирования Array, который копирует только указатель, а не элементы.

template<typename T>
struct Array<Array<T> >{
    explicit Array(size_t size, Array<T> value=Array<T>(size_t(0))) : array(new Arr<Array<T> >(size,value)){
        for(unsigned i=1;i<size;++i)
            (*this)[i]=value.copy();// Поэлементно копируем массив value в элементы нашего массива.
    }
//  Остальное всё тоже самое. Скопировать пришлось из-за особенностей частичной специализации.
    Array(Arr<Array<T> >* copy) : array(copy){ }
    Array(const Array<Array<T> >& copy) : array(copy.array){ }    
    Array() : array(0) {}


    void resize(size_t n){
        if(n==this->size())
                return;
        if(n>this->size()){
            this->subjoin(n-this->size());
        }else{
            Array<Array<T> > last(*this),penult;
            if(this->array->data_next.array==0){
                this->array=0;
            }else{    
                while(last.array->data_next.array!=0){
                    penult.mirror(last);
                    last.mirror(last.array->data_next);
                    if(this->size()-last.array->data_next.size()>=n){
                        penult.array->data_next.array=0;
                        break;
                    }       
                }    
                if(last.array->data_next.array==0)
                    penult.array->data_next.array=0;
            }    
            this->resize(n);
        }
    }

    Array& operator=(const Array<Array<T> > &A){
        if(A.array==0)
            this->array=0;
        else{
            this->array=0;
            this->resize(A.size());
            for(unsigned i=0;i<this->size();++i)
                (*this)[i]=A[i];
        }
        return *this;
    }

    Array<Array<T> > copy(){
        if(this->array==0)
            return Array<Array<T> >();
        else{
            Array<Array<T> > result(this->size());
            for(unsigned i=0;i<result.size();++i)
                result[i]=(*this)[i];
            return result;
        }
    }


    void mirror(const Array<Array<T> >& A){
        this->array=A.array;
    }

    void subjoin(const Array<Array<T> >& A){
        if(array!=0)
            array->subjoin(A);
        else  
            this->mirror(A); 
    }

    void subjoin(size_t size, const Array<T>& value=Array<T>(size_t(0))){
        if(array!=0)
            array->subjoin(size, value);
        else 
            this->mirror(Array<Array<T> >(size, value));
    }


    Array<T>& operator[](size_t i){return (*array)[i];}
    const Array<T>& operator[](size_t i) const{return (*array)[i];}
    size_t size() const{
        if(array!=0)
            return (this->array)->size();
        else
            return 0;
    }
    shared_ptr<Arr<Array<T> > > array;

};

//  Дальше идёт реализация различных операций с массивами.

template<typename T>
Array<T> operator+(Array<T> left, Array<T> right){
    if(left.size()!=right.size()){
        std::cout << "Разные размеры складываемых массивов!!!" << std::endl;
        exit(0);
    }
    Array<T> result(left.size());
    for(unsigned i=0;i<result.size();++i)
        result[i]=std::move(left[i]+right[i]);
    return result;
}


template<typename T>
Array<Array<T> > operator+(Array<Array<T> > left, Array<Array<T> > right){
    if(left.size()!=right.size()){
        std::cout << "Разные размеры складываемых массивов!!!" << std::endl;
        exit(0);
    }
    Array<Array<T> > result(left.size());
    for(unsigned i=0;i<result.size();++i)
        result[i].mirror(left[i]+right[i]);
    return result;
}

template<typename T>
Array<T> operator+=(Array<T> left, Array<T> right){
    left.mirror(left+right);
    return left;
}


template<typename T, typename U>
Array<T> operator*(Array<T> left, const U &right){
    Array<T> result(left.size());
    for(unsigned i=0;i<left.size();++i)
        result[i]=std::move(left[i]*right);
    return result;
}


template<typename T, typename U>
Array<Array<T> > operator*(Array<Array<T> > left, const U &right){
    Array<Array<T> > result(left.size());
    for(unsigned i=0;i<left.size();++i)
        result[i].mirror(left[i]*right);
    return result;
}

template<typename T, typename U>
Array<T> operator*(const U &right, Array<T> left){
    Array<T> result(left.size());
    for(unsigned i=0;i<left.size();++i)
        result[i]=std::move(left[i]*right);
    return result;
}


template<typename T, typename U>
Array<Array<T> > operator*(const U &right,Array<Array<T> > left){
    Array<Array<T> > result(left.size());
    for(unsigned i=0;i<left.size();++i)
        result[i].mirror(left[i]*right);
    return result;
}

template<typename T, typename U>
Array<T> operator/(Array<T> left, const U &right){
    Array<T> result(left.size());
    for(unsigned i=0;i<left.size();++i)
        result[i]=std::move(left[i]/right);
    return result;
}


template<typename T, typename U>
Array<Array<T> > operator/(Array<Array<T> > left, const U &right){
    Array<Array<T> > result(left.size());
    for(unsigned i=0;i<left.size();++i)
        result[i].mirror(left[i]/right);
    return result;
}

template<typename T, typename U>
Array<T> operator/(const U &right, Array<T> left){
    Array<T> result(left.size());
    for(unsigned i=0;i<left.size();++i)
        result[i]=right/left[i];
    return result;
}

template<typename T, typename U>
Array<Array<T> > operator/(const U &right,Array<Array<T> > left){
    Array<Array<T> > result(left.size());
    for(unsigned i=0;i<left.size();++i)
        result[i].mirror(right/left[i]);
    return result;
}

template<typename T, typename U>
Array<T> operator*=(Array<T> left, const U &right){
    left.mirror(left*right);
    return left;
}
template<typename T, typename U>
Array<T> operator/=(Array<T> left, const U &right){
    left.mirror(left/right);
    return left;
}