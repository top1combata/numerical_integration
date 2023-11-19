#include "functions.h"



/*
   b
   /
   |   x^k
   | --------- dx
   /  x^alpha
   a
*/
val_t moment(size_t k, val_t a, val_t b, val_t alpha)
{
    val_t power = k-alpha+1;
    return (pow(b, power) - pow(a, power)) / power;
}


std::vector<val_t> poly_interpol_quadrature(const std::vector<val_t>& nodes, val_t a, val_t b, val_t alpha)
{
    int n = nodes.size();

    Matrix<val_t> A = vandermonde(nodes, n-1, true);

    Row<val_t> c(n);
    for (int i = 0; i < n; i++)
        c[i] = moment(i,a,b,alpha);

    return LU_solve(A,c).getData();
}


std::vector<val_t> find_gaussian_nodes(val_t a, val_t b, val_t alpha, size_t n)
{
    std::vector<val_t> moments(2*n);
    for (size_t i = 0; i < 2*n; i++)
        moments[i] = moment(i,a,b,alpha);

    Matrix<val_t> M(n,n);
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            M[i][j] = moments[i+j];

    std::vector<val_t> tmp = LU_solve(M,-Row(std::vector(moments.begin()+n, moments.begin()+2*n))).getData();

    tmp.push_back(1);
    std::reverse(tmp.begin(), tmp.end());
    Polynomial<val_t> node_poly = tmp;

    auto res = find_roots(node_poly, a, b);
    std::sort(res.begin(), res.end());
    return res;
}



val_t integrate_riemann(const std::function<val_t(val_t)>& f, val_t l, val_t r, size_t n)
{
	val_t sum = 0, step = (r-l)/(n+1);
	for (size_t i = 0; i < n-1; i++)
		sum += f(l += step);
	return sum*step;
}


// n nodes
val_t integrate_newton(const std::function<val_t(val_t)>& f, val_t a, val_t b, val_t alpha, size_t n)
{
    if (a > b)
        return 0;
    std::vector<val_t> nodes(n);
    val_t h = (b-a)/(n-1);
    for (size_t i = 0; i < n; i++)
            nodes[i] = a + h*i;
    
    std::vector<val_t> coeffs = poly_interpol_quadrature(nodes, a, b, alpha);
    for (val_t cf : coeffs)
        if (cf < 0)
        {
            std::cout << "NEWTON NEGATIVE WEIGHT\n";
            exit(-3);
        }

    val_t res = 0;
    for (size_t i = 0; i < n; i++)
            res += coeffs[i] * f(nodes[i]);

    return res;
}




val_t integrate_gauss(const std::function<val_t(val_t)>& f, val_t a, val_t b, val_t alpha, size_t n)
{
    std::vector<val_t> nodes = find_gaussian_nodes(a, b, alpha, n);

    if (nodes.size() != n)
    {
        std::cout << "GAUSS NO ROOTS\n";
        exit(-1);
    }

    for (size_t i = 0; i < n; i++)
        if (nodes[i] < a || nodes[i] > b)
        {
            std::cout << "GAUSS ROOTS OUT OF BOUNDS\n";
            exit(-2);
        }
        

    std::vector<val_t> coeffs = poly_interpol_quadrature(nodes, a, b, alpha);

    for (val_t cf : coeffs)
        if (cf < 0)
        {
            std::cout << "GAUSS NEGATIVE WEIGHT\n";
            exit(-3);
        }

    val_t res = 0;
    for (size_t i = 0; i < n; i++)
            res += coeffs[i] * f(nodes[i]);

    return res;
}



/*
   b
   /
   |   f(x)
   | --------- dx
   /  x^alpha
   0
*/
static val_t integrate_qaws_wrapper(const std::function<val_t(val_t)>& f, val_t b, val_t alpha, val_t h, bool gauss, size_t k)
{
    size_t n = ceil(double(b/h));
    h = b/n;

#ifdef PARALLEL
    omp_set_num_threads(NUM_THREADS);
    
    val_t sums[NUM_THREADS];
    for (val_t& sum : sums)
        sum = 0;

    #pragma omp parallel
    {
        size_t i = omp_get_thread_num();

        val_t a = h*(n/NUM_THREADS * i + std::min(i, n%NUM_THREADS));
        size_t intervals = n/NUM_THREADS + (i < n%NUM_THREADS);

        for (size_t j = 0; j < intervals; j++)
        {
            sums[i] += (gauss ? integrate_gauss(f, a, a+h, alpha, k) : integrate_newton(f, a, a+h, alpha, k));
            a += h;
        }
    }
    val_t sum = 0;
    for (int i = 0; i < NUM_THREADS; i++)
        sum += sums[i];
    return sum;
#else

    val_t sum = 0, a = 0;
    for (size_t i = 0; i < n; i++)
    {
        sum += (gauss ? integrate_gauss(f, a, a+h, alpha, k) : integrate_newton(f, a, a+h, alpha, k));
        a += h;
    }
    return sum;
#endif
}



val_t integrate_qaws(const std::function<val_t(val_t)>& f, val_t a, val_t b, val_t alpha, val_t beta, val_t h, bool gauss, size_t k)
{
    if (alpha == 0)
        return integrate_qaws_wrapper([f,b](val_t x){return f(b-x);}, b-a, beta, h, gauss, k);
    
    if (beta == 0)
        return integrate_qaws_wrapper([f,a](val_t x){return f(x+a);}, b-a, alpha, h, gauss, k);

    return
        integrate_qaws([f,b,beta](val_t x){return f(x)/pow(b-x,beta);}, a, (a+b)/2, alpha, 0, h, gauss, k)+
        integrate_qaws([f,a,alpha](val_t x){return f(x)/pow(x-a,alpha);}, (a+b)/2, b, 0, beta, h, gauss, k);
}