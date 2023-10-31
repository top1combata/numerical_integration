#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include "row.h"

template<class Number>
class Matrix
{

private:

    int m_cols;
    int m_rows;
    std::vector<Row<Number>> m_rows_data;
    void swap_rows(int i, int j);

public:
    Matrix();
    Matrix(int); // unit matrix n x n
    Matrix(int, int, Number = Number());
    Matrix(const std::vector<std::vector<Number>>&);
    Matrix(const std::initializer_list<std::initializer_list<Number>>&);
    Matrix(const Matrix<Number>&);

    std::vector<std::vector<Number>> getData(void) const;

    Matrix<Number> operator+(const Matrix<Number>&) const;
    Matrix<Number> operator-(const Matrix<Number>&) const;
    Matrix<Number> operator*(const Matrix<Number>&) const;
    Row<Number> operator*(const Row<Number>&) const;
    Matrix<Number> operator*(Number) const;
    Matrix<Number> operator^(int) const;

    size_t rows() const;
    size_t cols() const;
    Row<Number>& operator[](int);
    const Row<Number>& operator[](int) const;

    Row<Number> getColumn(int) const;
};



template<class Number>
Matrix<Number>::Matrix()
{
    m_rows = 1;
    m_cols = 1;
    m_rows_data = {Row<Number>()};
}


template<class Number>
Matrix<Number>::Matrix(int _size)
{
    if (_size <= 0) return;

    m_cols = _size;
    m_rows = _size;

    m_rows_data.resize(_size);

    for (int i = 0; i < _size; i++)
    {
        Row<Number> tmp(_size, 0);
        tmp[i] = 1;
        m_rows_data[i] = tmp;
    }
}


template<class Number>
Matrix<Number>::Matrix(int rows, int cols, Number num)
{
    if (rows <= 0 || cols <= 0) return;

    m_rows = rows;
    m_cols = cols;

    m_rows_data = std::vector<Row<Number>>(rows, Row<Number>(cols, num));
}


template<class Number>
Matrix<Number>::Matrix(const std::vector<std::vector<Number>>& data)
{
    int cols = data[0].size();
    bool f = 1;
    for (std::vector<Number> row : data)
        if (row.size() != cols) f = 0;
    if (!f)
    {
        std::cout << "\nWarning ! rows have different size\n";
        m_rows = 1;
        m_cols = 1;
        m_rows_data = {Row<Number>()};
        return;
    }

    m_rows = data.size();
    m_cols = cols;
    m_rows_data.resize(m_rows);


    for (int i = 0; i < m_rows; i++)
    {
        m_rows_data[i] = Row<Number>(data[i]);
    }

}


template<class Number>
Matrix<Number>::Matrix(const std::initializer_list<std::initializer_list<Number>>& data)
{
    size_t cols = (*data.begin()).size();
    bool f = 1;
    for (std::initializer_list<Number> row : data)
        if (row.size() != cols) f = 0;
    if (!f)
    {
        std::cout << "\nWarning ! rows have different size\n";
        m_rows = 1;
        m_cols = 1;
        m_rows_data = {Row<Number>()};
        return;
    }

    m_rows = data.size();
    m_cols = cols;
    //m_rows_data.resize(m_rows);


    for (const std::initializer_list<Number>& row : data)
    {
        m_rows_data.push_back(Row<Number>(row));
    }

}


template<class Number>
Matrix<Number>::Matrix(const Matrix<Number>& A)
{
    m_rows = A.m_rows;
    m_cols = A.m_cols;
    m_rows_data = A.m_rows_data;
}


template<class Number>
void Matrix<Number>::swap_rows(int i, int j)
{
    std::swap(m_rows_data[i], m_rows_data[j]);
}


/*
std::ostream& operator<<(std::ostream& os, const Matrix& A)
{
    os << '\n';
    int width = 0;

    for (int i = 0; i < A.m_rows; i++)
    {
        width = std::max(width, A.m_rows_data[i].getMaxNumLength() + 1);
    }


    for (int i = 0; i < A.m_rows; i++)
    {
        os << ' ';
        for (int j = 0; j < A.m_cols; j++)
        {
            int l = len(A[i][j]);
            int k = 0;

            for (k; k < (width-l)/2; k++)
                os << ' ';

            os << A[i][j];

            k += l;

            for (k; k <= width; k++)
                os << ' ';
        }
        os << '\n';
    }

    return os << '\n';
}
*/

template<class Number>
std::ostream& operator<<(std::ostream& os, const Matrix<Number>& A)
{
    os << '\n';
    for (size_t i = 0; i < A.rows(); i++)
    {
        for (size_t j = 0; j < A.cols(); j++)
            os << A[i][j] << ' ';
        os << '\n';
    }
    return os << '\n';
}


template<class Number>
Matrix<Number> gram(const std::vector<Row<Number>>& data)
{
    int n = data.size();
    Matrix res(n,n,0);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res[i][j] = data[i] * data[j];
        }
    }

    return res;
}


template<class Number>
std::vector<std::vector<Number>> Matrix<Number>::getData(void) const
{
    int rows = m_rows, cols = m_cols;
    std::vector<std::vector<Number>> res(rows, std::vector<Number>(cols));
    for (int i = 0; i < rows; i++)
    {
        res[i] = m_rows_data[i].getData();
    }
    return res;
}

template<class Number>
size_t Matrix<Number>::rows() const
{
    return m_rows;
}

template<class Number>
size_t Matrix<Number>::cols() const
{
    return m_cols;
}

template<class Number>
Row<Number>& Matrix<Number>::operator[](int index)
{
    if (index >= m_rows)
    {
        std::cout << "\nWarning ! index exceed rows count\n";
        return m_rows_data[0];
    }

    return m_rows_data[index];
}


template<class Number>
const Row<Number>& Matrix<Number>::operator[](int index) const
{
    if (index >= m_rows)
    {
        std::cout << "\nWarning ! index exceed rows count\n";
        return m_rows_data[0];
    }

    return m_rows_data[index];
}


template<class Number>
Number det(const Matrix<Number>& A)
{
    if (A.cols() != A.rows())
    {
        std::cout << "\nWarning ! matrix must be square\n";
        return 0;
    }
    Matrix<Number> tmp = A;

    Number res(1);
    int n = A.rows();
    Number a;

    for (int col = n-1; col > 0; col--)
    {
        for (int row = 0; row <= col; row++)
        {
            if (row == col && tmp[row][col] == Number())
                return 0;

            if (tmp[row][col] != Number())
            {
                if (row != col)
                {
                    res = -res;
                    std::swap(tmp[row], tmp[col]);//tmp.swap_rows(row, col);
                }
                break;
            }
        }

        a = tmp[col][col];

        for (int row = 0; row < col; row++)
        {
            tmp[row] -= tmp[col] * (tmp[row][col]/a);
        }
    }

    for (int i = 0; i < n; i++)
        res *= tmp[i][i];

    return res;
}


template<class Number>
Matrix<Number> inverse(const Matrix<Number>& A)
{
    if (A.cols() != A.rows() || det(A) == Number())
    {
        std::cout << "\nWarning ! inverse matrix doesn't exist\n";
        return A;
    }

    int n = A.rows();

    Matrix<Number> tmp(n, 2*n, 0);
    Matrix<Number> unit(n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            tmp[i][j] = A[i][j];
            tmp[i][j+n] = unit[i][j];
        }
    }

    Number a;

    //------------TRANSFORMATE MATRIX TO TRIANGLE--------

    for (int col = n-1; col >= 0; col--)
    {
        for (int row = 0; row <= col; row++)
        {
            if (tmp[row][col] != Number())
            {
                std::swap(tmp[row], tmp[col]);//tmp.swap_rows(row, col);
                break;
            }
        }

        a = tmp[col][col];

        for (int row = 0; row < col; row++)
        {
            tmp[row] -= tmp[col] * (tmp[row][col]/a);
        }

        tmp[col] /= a;
    }

    //---------TRANSFORMATE MATRIX TO UNIT------

    for (int row = 1; row < n; row++)
    {
        for (int col = 0; col < row; col++)
        {
            tmp[row] -= tmp[col] * tmp[row][col];
        }
    }

    Matrix<Number> res(n,n,0);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res[i][j] = tmp[i][n+j];
        }
    }

    return res;
}


template<class Number>
Matrix<Number> transpose(const Matrix<Number>& A)
{
    Matrix<Number> res(A.cols(), A.rows(), 0);

    for (size_t i = 0; i < A.cols(); i++)
    {
        for (size_t j = 0; j < A.rows(); j++)
        {
            res[i][j] = A[j][i];
        }
    }

    return res;
}


template<class Number>
Matrix<Number> Matrix<Number>::operator+(const Matrix<Number>& A) const
{
    if (m_rows != A.m_rows || m_cols != A.m_cols)
    {
        std::cout << "\nWarning ! matrices have different size\n";
        return *this;
    }
    int cols = m_cols, rows = m_rows;

    Matrix<Number> res(rows, cols, 0);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            res[i][j] = m_rows_data[i][j] + A[i][j];
        }
    }
    return res;
}


template<class Number>
Matrix<Number> Matrix<Number>::operator-(const Matrix<Number>& A) const
{
    if (m_rows != A.m_rows || m_cols != A.m_cols)
    {
        std::cout << "\nWarning ! matrices have different size\n";
        return *this;
    }
    int cols = m_cols, rows = m_rows;

    Matrix<Number> res(rows, cols, 0);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            res[i][j] = m_rows_data[i][j] - A[i][j];
        }
    }
    return res;
}


template<class Number>
Matrix<Number> Matrix<Number>::operator*(const Matrix<Number>& A) const
{
    if (m_cols != A.m_rows)
    {
        std::cout << "Warning ! matrices are not compatible\n";
        return *this;
    }

    int rows = m_rows, cols = A.m_cols;

    Matrix<Number> res(rows, cols, 0);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            res[i][j] = m_rows_data[i] * A.getColumn(j);
        }
    }

    return res;
}

template<class Number>
Row<Number> Matrix<Number>::operator*(const Row<Number>& row) const
{
    if (m_cols != row.size())
    {
        std::cout << "cols != row size\n";
        return Row<Number>();
    }

    Row<Number> res(m_rows, 0);

    for (int i = 0; i < m_rows; i++)
        for (int j = 0; j < m_cols; j++)
            res[i] += m_rows_data[i][j] * row[j];
    return res;
}


template<class Number>
Matrix<Number> Matrix<Number>::operator*(Number num) const
{
    int rows = m_rows, cols = m_cols;
    Matrix<Number> res(rows, cols, 0);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
            res[i][j] = m_rows_data[i][j] * num;
    }

    return res;
}

template<class Number>
Matrix<Number> Matrix<Number>::operator^(int power) const
{
    if (m_cols != m_rows)
    {
        std::cout << "Warning ! matrix cannot be powered up\n";
        return *this;
    }

    Matrix<Number> mult;
    if (power < 0)
    {
        power = -power;
        mult = inverse(*this);
    }
    else
        mult = *this;

    Matrix<Number> res(m_rows);
    for (int i = 0; i < power; ++i)
        res = res * mult;

    return res;
}


template<class Number>
Matrix<Number> operator*(Number num, const Matrix<Number>& A)
{
    return A*num;
}



template<class Number>
Row<Number> Matrix<Number>::getColumn(int index) const
{
    Row<Number> res(m_rows, 0);
    for (int i = 0; i < m_rows; i++)
        res[i] = m_rows_data[i][index];

    return res;
}


template<class Number>
Row<Number> linear_system_solve(const Matrix<Number>& A, const Row<Number>& b)
{
    int n = b.size();

    // checking
    /*
    if (det(A) == 0)
    {
        std::cout << "det == 0\n";
        return Row<Number>();
    }
    */

    Matrix<Number> W(n, n+1);
    for (int i = 0; i < n; i++)
    {
        std::copy(&A[i][0], &A[i][n-1]+1,&W[i][0]);
        W[i][n] = b[i];
    }

    // forward
    for (int j = 0; j < n; j++)
    {
        int maxx = j;
        for (int i = j+1; i < n; i++)
            maxx = abs(W[i][j]) > abs(W[i][maxx]) ? j : maxx;

        std::swap(W[maxx], W[j]);
        Number tmp = W[j][j];
        for (int s = 0; s <= n; s++)
            W[j][s] /= tmp;

        for (int i = j+1; i < n; i++)
        {
            tmp = W[i][j];
            for (int s = 0; s <= n; s++)
                W[i][s] -= W[j][s]*tmp;
        }
    }

    // inverse
    Row<Number> res(n);
    for (int i = n-1; i >= 0; i--)
    {
        for (int j = i+1; j < n; j++)
            W[i][n] -= res[j]*W[i][j];
        res[i] = W[i][n]/W[i][i];
    }

    return res;
}



template<class Number>
Row<Number> LU_solve(const Matrix<Number>& A, const Row<Number>& b)
{
    int n = b.size();

    Matrix<Number> L(n,n,0), U(n,n,0);
    for (int i = 0; i < n; i++)
        L[i][i] = 1;
    
    // finding L and U
    for (int i = 0; i < n; i++)
    {
        // finding L part
        for (int j = 0; j < i; j++)
        {
            L[i][j] = A[i][j];
            for (int s = 0; s < j; s++)
                L[i][j] -= L[i][s]*U[s][j];
            L[i][j] /= U[j][j];
        }
        // finding U part
        for (int j = i; j < n; j++)
        {
            U[i][j] = A[i][j];
            for (int s = 0; s < i; s++)
                U[i][j] -= L[i][s]*U[s][j];
            U[i][j] /= L[i][i];
        }
    }

    Row<Number> y(n);
    for (int i = 0; i < n; i++)
    {
        y[i] = b[i];
        for (int j = 0; j < i; j++)
            y[i] -= L[i][j]*y[j];
        y[i] /= L[i][i];
    }

    Row<Number> res(n);
    for (int i = n-1; i >= 0; i--)
    {
        res[i] = y[i];
        for (int j = n-1; j > i; j--)
            res[i] -= U[i][j]*res[j];
        res[i] /= U[i][i];
    }
  
    return res;
}





template<class Number>
Matrix<Number> vandermonde(const std::vector<Number>& _values, size_t n, bool _transpose = 0)
{
    size_t m = _values.size();

    if (_transpose)
    {
        Matrix<Number> res(n+1, m);
        for (size_t j = 0; j < m; j++)
            res[0][j] = 1;
        for (size_t i = 1; i < n+1; i++)
            for (size_t j = 0; j < m; j++)
                res[i][j] = res[i-1][j] * _values[j];
        return res;
    }
    else
    {
        Matrix<Number> res(m, n+1);
        for (size_t i = 0; i < m; i++)
            res[i][0] = 1;
        for (size_t i = 0; i < m; i++)
            for (size_t j = 1; j < n+1; j++)
                res[i][j] = res[i][j-1] * _values[i];
        return res;
    }
}

#endif // MATRIX_H_INCLUDED
