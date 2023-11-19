#ifndef POLYNOMIAL_H_INCLUDED
#define POLYNOMIAL_H_INCLUDED

#include <iostream>
#include <vector>


template<class T>
class Polynomial
{
private:

    std::vector<T> m_coeffs = {T()};

    Polynomial(std::vector<T>, bool);
    void delete_zeroes();
    Polynomial base_multiply(size_t) const;

public:

    Polynomial();
    Polynomial(std::vector<T>);
    Polynomial(std::initializer_list<T>);
    Polynomial(T);


    size_t deg() const;

    T operator()(T) const;
    T operator[](size_t) const;
    T& operator[](size_t);
    Polynomial operator-() const;
    Polynomial& operator/=(T);
    Polynomial operator/(T) const;

    Polynomial& operator+=(const Polynomial&);

    Polynomial& operator-=(const Polynomial&);

    Polynomial operator*(const Polynomial&) const;
    Polynomial& operator*=(const Polynomial&);
    Polynomial& operator*=(T);

    Polynomial& operator/=(const Polynomial&);
    Polynomial  operator/(const Polynomial&) const;
    Polynomial  operator%(const Polynomial&) const;

    Polynomial  operator^(size_t) const;

    //friend std::ostream& operator<<<>(std::ostream&, const Polynomial&);
};



// --------------------------------- CONSTRUCTORS----------------------------------------

template<class T>
Polynomial<T>::Polynomial() = default;


template<class T>
Polynomial<T>::Polynomial(std::vector<T> coeffs, bool) : m_coeffs(coeffs.size())
{
    for (size_t i = 0; i  < coeffs.size(); i++)
        m_coeffs[i] = coeffs[coeffs.size()-i-1];
}


template<class T>
Polynomial<T>::Polynomial(std::vector<T> coeffs) : Polynomial(coeffs, 1)
{
    delete_zeroes();
}


template<class T>
Polynomial<T>::Polynomial(std::initializer_list<T> lst) : m_coeffs(lst.size())
{
    auto iter = lst.end()-1;
    for (size_t i = 0; i < lst.size(); i++)
    {
        m_coeffs[i] = *iter;
        iter--;
    }
    delete_zeroes();
}


template<class T>
Polynomial<T>::Polynomial(T value) : m_coeffs({value}) {}


//-------------------------------------MEMBER-FUNTCTIONS------------------------------------------


template<class T>
void Polynomial<T>::delete_zeroes()
{
    while (m_coeffs.size() > 1 && m_coeffs[m_coeffs.size()-1] == T(0))
        m_coeffs.pop_back();
}


template<class T>
Polynomial<T> Polynomial<T>::base_multiply(size_t deg) const
{
    Polynomial res(std::vector<T>(m_coeffs.size() + deg, T(0)), 1);
    for (size_t i = 0; i < m_coeffs.size(); i++)
        res.m_coeffs[i+deg] = m_coeffs[i];
    return res;
}


template<class T>
size_t Polynomial<T>::deg() const
{
    return m_coeffs.size()-1;
}


// --------------------------------------------- OPERATORS --------------------------------------------------------------------


template<class T>
T Polynomial<T>::operator()(T value) const
{
    T res(0);
    for (size_t i = m_coeffs.size()-1; i  < m_coeffs.size(); i--)
        res = res*value + m_coeffs[i];
    return res;
}


template<class T>
T Polynomial<T>::operator[](size_t i) const
{
    return m_coeffs[i];
}


template<class T>
T& Polynomial<T>::operator[](size_t i)
{
    return m_coeffs[i];
}


template<class T>
Polynomial<T> Polynomial<T>::operator-() const
{
    Polynomial res(*this);
    for (auto& cf : res.m_coeffs)
        cf = -cf;
    return res;
}


template<class T>
Polynomial<T>& Polynomial<T>::operator/=(T value)
{
    if (value == T(0)) return *this;

    for (auto& cf : m_coeffs)
        cf /= value;

    return *this;
}


template<class T>
Polynomial<T> Polynomial<T>::operator/(T value) const
{
    Polynomial cpy(*this);
    return cpy/=value;
}


template<class T, class U>
Polynomial<T> operator*(U value, const Polynomial<T>& poly)
{

    return poly*T(value);
}

// --------------------------------------------------------------------------------------------------------

template<class T>
Polynomial<T>& Polynomial<T>::operator+=(const Polynomial<T>& other)
{
    if (this == &other)
        return *this *= 2;

    if (other.deg() > deg())
        m_coeffs.resize(other.deg()+1);

    size_t m = std::min(other.deg(), deg())+1;
    for (size_t i = 0; i < m; i++)
        m_coeffs[i] += other.m_coeffs[i];

    delete_zeroes();

    return *this;
}


template<class T>
Polynomial<T>& Polynomial<T>::operator-=(const Polynomial<T>& other)
{
    Polynomial tmp = -other;
    return *this += tmp;
}


template<class T>
Polynomial<T> operator+(const Polynomial<T>& p1, const Polynomial<T>& p2)
{
    Polynomial cpy(p1);
    return cpy += p2;
}


template<class T>
Polynomial<T> operator-(const Polynomial<T>& p1, const Polynomial<T>& p2)
{
    Polynomial cpy(p1);
    return cpy -= p2;
}


template<class T>
Polynomial<T>& Polynomial<T>::operator*=(const Polynomial<T>& other)
{
    *this = *this * other;
    return *this;
}


template<class T>
Polynomial<T>& Polynomial<T>::operator*=(T num)
{
    if (num == 0)
        m_coeffs = {0};
    for (T& t : m_coeffs)
        t *= num; 
    return *this;
}


template<class T>
Polynomial<T> Polynomial<T>::operator*(const Polynomial<T>& other) const
{
    Polynomial res(std::vector<T> (other.deg()+deg()+1, T(0)), 1);

    for (size_t i = 0; i <= deg(); i++)
    {
        for (size_t j = 0; j <= other.deg(); j++)
        {
            res.m_coeffs[i+j] += m_coeffs[i] * other.m_coeffs[j];
        }
    }
    return res;
}


template<class T>
Polynomial<T>& Polynomial<T>::operator/=(const Polynomial& other)
{
    *this = *this / other;
    return *this;
}


template<class T>
Polynomial<T> Polynomial<T>::operator/(const Polynomial<T>& other) const
{
    if (other.m_coeffs.size() == 1 && other.m_coeffs[0] == T(0)) return *this;
    if (other.deg() > deg()) return Polynomial(T(0));


    Polynomial res(std::vector<T>(deg()-other.deg()+1), 1);
    if (other.deg() == 1)
    {
        T x = -other.m_coeffs[0] / other.m_coeffs[1], tmp(0);
        for (size_t i = 0; i < deg(); i++)
        {
            tmp = tmp*x + m_coeffs[deg()-i];
            res.m_coeffs[deg()-i-1] = tmp;
        }
        return res;
    }

    Polynomial cpy(*this);

    for (size_t i = res.m_coeffs.size()-1; i < res.m_coeffs.size(); i--)
    {
        res.m_coeffs[i] = cpy.m_coeffs[i+other.deg()] / other.m_coeffs[other.deg()];
        cpy -= (res.m_coeffs[i] * other).base_multiply(i);
    }

    return res;
}




template<class T>
Polynomial<T> Polynomial<T>::operator%(const Polynomial<T>& other) const
{
    if (other.deg() == 1)
        return {(*this)(-other[0]/other[1])};

    Polynomial<T> cpy(*this);
    while (cpy.deg() >= other.deg())
    {
        T cf = cpy[cpy.deg()] / other[other.deg()];
        size_t prev_deg = cpy.deg();
        cpy -= cf * other.base_multiply(cpy.deg()-other.deg());
        if (cpy.deg() == prev_deg)
            cpy.m_coeffs.pop_back();

        //std::cout << "wh " << cpy << " ww\n";
    }
    return cpy;
}



template<class T>
Polynomial<T> Polynomial<T>::operator^(size_t power) const
{
    if (power == 0) return T(1);
    if (power == 1) return *this;
    Polynomial res = *this^(power/2);
    return (power%2 ? *this : T(1))*res*res;
}


template<class T>
std::ostream& operator<<(std::ostream& os, const Polynomial<T>& poly)
{
    size_t n = poly.deg();

    os << poly[n];
    if (n > 1) os << "*x^" << n;
    else if (n == 1) os << "*x";

    for (size_t i = n-1; i < n ; i--)
        if (poly[i] != T(0))
        {
            if (poly[i] < T(0)) os << " - ";
            else os << " + ";

            os << (poly[i] > T(0) ? poly[i] : -poly[i]);
            if (i > 1) os << "*x^" << i ;
            else if(i == 1) os << "*x";
        }

    return os;
}


template<class T>
Polynomial<T> D(const Polynomial<T>& poly, size_t n = 1)
{
    if (n == 0) return poly;
    if (n > poly.deg()) return T(0);

    Polynomial res(std::vector<T>(poly.deg()+1-n, 1));
    for (size_t i = 0; i < poly.deg()+1-n; i++)
    {
        res[i] = poly[i+n];
        for (size_t j = 0; j < n; j++)
            res[i] *= i+n-j;
    }
    return res;
}


#endif // POLYNOMIAL_H_INCLUDED
