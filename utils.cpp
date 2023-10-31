#include "utils.h"


val_t my_sin(val_t x)
{
#ifndef MULTIPRECISION
    return std::sin(x);
#else
    x += PI;
    x -= 2*PI * floor(x/PI/2);
    x -= PI;
    if (x > PI/2)  x = PI-x;
    if (x < -PI/2) x = -x-PI;

    val_t tmp(x), prev(0), curr(x);

    for (size_t i = 2; prev - curr != val_t(0); i += 2)
    {
        prev = curr;
        tmp *= -(x*x)/i/(i+1);
        curr += tmp;
    }
    return curr;
#endif
}

val_t my_cos(val_t x)
{
    return my_sin(PI/2-x);
}

val_t my_exp(val_t x)
{
#ifndef MULTIPRECISION
    return std::exp(x);
#else
    val_t prev(0), curr(1), tmp(1);
    for (size_t n = 1; prev != curr; n++)
    {
        prev = curr;
        tmp *= x/n;
        curr += tmp;
    }
    return curr;
#endif
}


val_t my_pow(val_t x, val_t y)
{
#ifndef MULTIPRECISION
    return std::pow(x,y);
#else
    return pow(x,y);

    if (x <= 0 || x < EPS)
        return 0;
    return my_exp(y * my_log(x));
#endif
}


val_t my_log(val_t x)
{
#ifndef MULTIPRECISION
    return std::log(x);
#else
    return log(x);
    size_t k = 0;
    while (x >= val_t(2))
    {
        x /= E;
        k++;
    }
    x -= val_t(1);

    val_t prev(1), curr(0), tmp(x);
    for (size_t i = 1; prev != curr; i++)
    {
        prev = curr;
        curr += tmp/i;
        tmp *= -x;
    }
    return curr+k;
#endif
}


val_t my_abs(val_t x)
{
    return (x > 0 ? x : -x);
}


val_t my_sqrt(val_t x)
{
#ifndef MULTIPRECISION
    return std::sqrt(x);
#else
    return sqrt(x);
    val_t curr = x/2, prev = x;
    while (my_abs(curr-prev) > curr*EPS)
    {
        prev = curr;
        curr = (curr + x/curr)/2;
    }
    return curr;
#endif
}


sturm_sequence::sturm_sequence(Polynomial<val_t> p)
{
    seq.push_back(p);
    seq.push_back(D(p));
    p = -seq[0]%seq[1];

    while (p.deg() > 0)
    {
        seq.push_back(p);
        size_t n = seq.size();
        p = -seq[n-2]%seq[n-1];
    }
    seq.push_back(p);
}



int sturm_sequence::operator()(val_t val) const
{
    int count = 0;

    for (size_t i = 1; i < seq.size(); i++)
    {
        if ((seq[i-1](val) < 0) ^ (seq[i](val) < 0))
            count++;
    }
    return count;
}



std::vector<val_t> roots_quadratic(const Polynomial<val_t>& poly)
{
    val_t a = poly[2], b = poly[1], c = poly[0];
    val_t d = b*b - 4*a*c;

    if (d < 0)
        return {};

    return {(-b + my_sqrt(d))/2/a, (-b - my_sqrt(d))/2/a};
}


std::vector<val_t> roots_cubic(const Polynomial<val_t>& poly)
{
    val_t a = poly[3], b = poly[2], c = poly[1], d = poly[0];
    std::complex<val_t> p = (3*a*c - b*b)/val_t(3)/a/a;
    std::complex<val_t> q = (2*b*b*b - 9*a*b*c + 27*a*a*d)/val_t(27)/a/a/a;

    std::complex<val_t> Q = p*p*p/val_t(27) + q*q/val_t(4);
    
    if (Q.real() >= 0)
        return {};

    std::complex<val_t> rot(val_t(-1)/2, my_sqrt(3)/2);
    std::complex<val_t> alpha = std::pow( std::sqrt(Q) - std::complex<val_t>(q/val_t(2)), val_t(1)/3);
    std::vector<val_t> roots(3);

    val_t shift = b/3/a;
    roots[0] = (alpha + -p/val_t(3)/alpha).real() - shift;
    roots[1] = (alpha*rot - p/val_t(3)/alpha/rot).real() - shift;
    roots[2] = (alpha/rot - p/val_t(3)/alpha*rot).real() - shift;
    return roots;
}


val_t root_binary_search(const Polynomial<val_t>& p, val_t l, val_t r)
{
    bool f = p(l) < 0;
    while (r-l > EPS*my_abs(l))
    //while (r-l > EPS*my_abs(l))
    {
        val_t m = (l+r)/2;
        if ((p(m) < 0) == f)
            l = m;
        else
            r = m;
    }
    return l;
}



void localize_roots(const Polynomial<val_t>& p, const sturm_sequence& seq, val_t l, val_t r, std::vector<val_t>& roots)
{
    int roots_count = seq(l) - seq(r);
    
    if (roots_count == 0)
        return;
    if (roots_count == 1)
        return roots.push_back(root_binary_search(p,l,r));
    
    val_t m = (l+r)/2;
    localize_roots(p, seq, l, m, roots);
    localize_roots(p, seq, m, r, roots);
}



std::vector<val_t> find_roots(const Polynomial<val_t>& p, val_t l, val_t r)
{
    size_t n = p.deg();
    if (n == 0)
        return {};
    if (n == 1)
        return {-p[0]};
    if (n == 2)
        return roots_quadratic(p);
    if (n == 3)
        return roots_cubic(p);
    /*
    std::vector<val_t> roots;
    sturm_sequence seq(p);

    roots.reserve(seq(l)-seq(r));
    localize_roots(p, seq, l, r, roots);
    return roots;
    */
    Polynomial<val_t> dp = D(p,1);

    val_t curr = (l+r)/2;
    val_t prev = curr-1;

    size_t i = 0;

    const int MAX_ITERS = 20000;

    for (; my_abs(curr-prev) > my_abs(curr*EPS) && i < MAX_ITERS; i++)
    {
        prev = curr;
        curr = curr - p(curr)/dp(curr);
    }
    if (i == MAX_ITERS)
        return {};

    std::vector<val_t> res = find_roots(p/Polynomial<val_t>{1,-curr}, l, r);
    res.push_back(curr);
    return res;
} 