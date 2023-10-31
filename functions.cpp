#include "functions.h"


val_t fp(val_t x)
{
    return f(x)*p(x);
}


val_t f(val_t x)
{
    return 4*my_cos(x/2)*my_exp(-2*x/3) + val_t(12)/5 * my_sin(9*x/2)*my_exp(x/8)+2;
}


val_t p(val_t x)
{
	return my_pow(x-A, -ALPHA) * my_pow(B-x, -BETA);
}




val_t integrate_riemann(const std::function<val_t(val_t)>& f, val_t l, val_t r, size_t n)
{
	val_t sum = 0, step = (r-l)/(n+1);
	for (size_t i = 0; i < n-1; i++)
		sum += f(l += step);
	return sum*step;
}


static val_t moment(size_t i, val_t l, val_t r)
{
    val_t power = i-BETA+1;
    return (my_pow(B-l, power) - my_pow(B-r, power)) / (power);
}


static std::vector<val_t> poly_interpol_quadrature(const std::vector<val_t>& nodes, val_t l, val_t r)
{
    int n = nodes.size();

    Matrix<val_t> A = vandermonde(nodes, n-1, true);

    Row<val_t> b(n);
    for (int i = 0; i < n; i++)
        b[i] = moment(i,l,r);

    //return linear_system_solve(A,b).getData();
    return LU_solve(A,b).getData();
}


static std::vector<val_t> find_gaussian_nodes(val_t l, val_t r, size_t n)
{
    std::vector<val_t> moments(2*n);
    for (size_t i = 0; i < 2*n; i++)
        moments[i] = moment(i,l,r);

    Matrix<val_t> M(n,n);
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            M[i][j] = moments[i+j];

    std::vector<val_t> tmp = linear_system_solve(M,-Row(std::vector(moments.begin()+n, moments.begin()+2*n))).getData();

    tmp.push_back(1);
    std::reverse(tmp.begin(), tmp.end());
    Polynomial<val_t> node_poly = tmp;

    return find_roots(node_poly, l, r);
}


// 

// n nodes
Solve integrate_newton(const std::function<val_t(val_t)>& f, val_t l, val_t r, size_t n)
{
        std::vector<val_t> nodes(n);
        val_t gap = (r-l)/(n-1);
        for (size_t i = 0; i < n; i++)
                nodes[i] = B-r + gap*i;
        
        std::vector<val_t> coeffs = poly_interpol_quadrature(nodes, l, r);

        val_t res = 0;
        for (size_t i = 0; i < n; i++)
                res += coeffs[i] * f(B-nodes[i]);

        val_t abs_sum = 0;
        for (size_t i = 0; i < n; i++)
                abs_sum += my_abs(coeffs[i]);
        return {res, abs_sum};
}




Solve integrate_gauss(const std::function<val_t(val_t)>& f, val_t l, val_t r, size_t n)
{
    std::vector<val_t> nodes = find_gaussian_nodes(l,r,n);

    if (nodes.size() < n)
    {
        std::cout << "\nNO ROOTS\n";
        return {0,0};
    }
    for (size_t i = 0; i < n; i++)
    {
        if (nodes[i] < B-r || nodes[i] > B-l)
        {
            std::cout << "\nROOTS OUT OF BOUNDS\n";
            return {0,0};
        }
    }
    std::vector<val_t> coeffs = poly_interpol_quadrature(nodes, l, r);
    for (auto cf : coeffs)
        if (cf < 0)
        {
            std::cout << "\nNEGATIVE \n";
            return {0,0};
        }

    val_t res = 0;
    for (size_t i = 0; i < n; i++)
            res += coeffs[i] * f(B-nodes[i]);

    val_t abs_sum = 0;
    for (size_t i = 0; i < n; i++)
            abs_sum += my_abs(coeffs[i]);

    return {res, abs_sum};
}




val_t integrate(const std::function<val_t(val_t)>& f, val_t l, val_t r, val_t h, bool gauss, size_t k)
{
    size_t n = ceil(double((r-l)/h));
    h = (r-l)/n;

    val_t sum = 0;
    for (size_t i = 0; i < n; i++)
        sum += (gauss ? integrate_gauss(f, l+i*h, l+(i+1)*h, k).res : integrate_newton(f, l+i*h, l+(i+1)*h, k).res);
        
    return sum;
}