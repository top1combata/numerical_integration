#ifndef ROW_H_INCLUDED
#define ROW_H_INCLUDED

#include <vector>
#include <cmath>
#include <iostream>
#include <complex>

template<class Number>
class Row
{
private:

    std::vector<Number> m_data = {Number(0)};

public:

    Row() = default;
    Row(int);
    Row(const std::vector<Number>&);
    Row(const std::initializer_list<Number>&);
    Row(int, Number);
    Row(const Row&);
    Row& operator=(const Row&);
    ~Row(void){}

    std::vector<Number> getData(void) const;
    void setData(const std::vector<Number>&);
    int size(void) const;

    Number& operator[](int);
    const Number& operator[](int) const;

    Row operator*(Number) const;

    Row operator-() const;
    Row operator+(const Row&) const;
    Row operator-(const Row&) const;
    Row operator/(Number) const;
    Number operator*(const Row&) const;

    void operator+=(const Row&);
    void operator-=(const Row&);
    void operator*=(Number);
    void operator/=(Number);

};

template<class Number>
Row<Number>::Row(int _size) : m_data(_size) {}


template<class Number>
Row<Number>::Row(const std::vector<Number>& data)
{
    m_data = data;
}

template<class Number>
Row<Number>::Row(const std::initializer_list<Number>& lst)
{
    m_data.clear();
    for (const Number& num : lst)
        m_data.push_back(num);
}


template<class Number>
Row<Number>::Row(int length, Number value)
{
    m_data = std::vector<Number>(length, value);
}


template<class Number>
Row<Number>::Row(const Row& row)
{
    m_data = row.m_data;
}


template<class Number>
Row<Number>& Row<Number>::operator=(const Row& row)
{
    m_data = row.m_data;
    return *this;
}


template<class Number>
std::ostream& operator<<(std::ostream& os, const Row<Number>& row)
{
    os << '\n';
    for (int i = 0; i < row.size(); i++)
        os << row[i] << ' ';
    os << '\n';

    return os;
}


template<class Number>
std::vector<Number> Row<Number>::getData(void) const
{
    return m_data;
}


template<class Number>
void Row<Number>::setData(const std::vector<Number>& data)
{
    m_data = data;
}


template<class Number>
int Row<Number>::size(void) const
{
    return m_data.size();
}


template<class Number>
Number& Row<Number>::operator[](int index)
{
    return m_data[index];
}


template<class Number>
const Number& Row<Number>::operator[](int index) const
{
    return m_data[index];
}


template<class Number>
Row<Number> Row<Number>::operator*(Number alpha) const
{
    Row res(m_data);

    for (Number& num : res.m_data)
        num *= alpha;

    return res;
}


template<class Number>
Row<Number> operator*(Number alpha, const Row<Number>& row)
{
    return row*alpha;
}

template<class Number>
Row<Number> Row<Number>::operator-() const
{
    return Row(size(),0)-*this;
}

template<class Number>
Row<Number> Row<Number>::operator+(const Row& row) const
{
    Row res(*this);

    if (m_data.size() != row.m_data.size())
    {
        std::cout << "\nWarning ! Rows have different size\n";
        return res;
    }

    for (int i = 0; i < m_data.size(); i++)
        res.m_data[i] += row.m_data[i];

    return res;
}


template<class Number>
Row<Number> Row<Number>::operator-(const Row<Number>& row) const
{
    Row res(*this);

    if (m_data.size() != row.m_data.size())
    {
        std::cout << "\nWarning ! Rows have different size\n";
        return res;
    }

    for (size_t i = 0; i < m_data.size(); i++)
        res.m_data[i] -= row.m_data[i];

    return res;
}


template<class Number>
Row<Number> Row<Number>::operator/(Number alpha) const
{
    Row res = *this;
    if (alpha == Number())
    {
        std::cout << "\nWarning ! deviding by zero\n";
        return res;
    }

    for (Number& num : res.m_data)
        num /= alpha;

    return res;
}


template<class Number>
Number Row<Number>::operator*(const Row<Number>& row) const
{
    if (m_data.size() != row.m_data.size())
    {
        std::cout << "\nWarning ! Rows have different size\n";
        return 0;
    }

    Number res = 0;

    for (int i = 0; i < m_data.size(); i++)
        res += m_data[i] * row.m_data[i];

    return res;
}


template<class Number>
void Row<Number>::operator+=(const Row<Number>& row)
{
    if (m_data.size() != row.m_data.size())
    {
        std::cout << "\nWarning ! Rows have different size\n";
    }

    for (int i = 0; i < m_data.size(); i++)
        m_data[i] += row.m_data[i];

}


template<class Number>
void Row<Number>::operator-=(const Row<Number>& row)
{
    if (m_data.size() != row.m_data.size())
    {
        std::cout << "\nWarning ! Rows have different size\n";
    }

    for (size_t i = 0; i < m_data.size(); i++)
        m_data[i] -= row.m_data[i];

}

template<class Number>
void Row<Number>::operator*=(Number alpha)
{
    for (Number& num : m_data)
        num *= alpha;
}


template<class Number>
void Row<Number>::operator/=(Number alpha)
{
    if (alpha == Number())
    {
        std::cout << "\nWarning ! deviding by zero\n";
        return;
    }

    for (Number& num : m_data)
        num /= alpha;
}
#endif // ROW_H_INCLUDED
