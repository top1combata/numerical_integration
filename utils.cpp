#include "utils.h"



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

    return {(-b + sqrt(d))/2/a, (-b - sqrt(d))/2/a};
}


std::vector<val_t> roots_cubic(const Polynomial<val_t>& poly)
{
    val_t a = poly[3], b = poly[2], c = poly[1], d = poly[0];
    std::complex<val_t> p = (3*a*c - b*b)/val_t(3)/a/a;
    std::complex<val_t> q = (2*b*b*b - 9*a*b*c + 27*a*a*d)/val_t(27)/a/a/a;

    std::complex<val_t> Q = p*p*p/val_t(27) + q*q/val_t(4);
    
    if (Q.real() >= 0)
        return {};

    std::complex<val_t> rot(val_t(-1)/2, sqrt(3)/2);
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
    while (r-l > EPS*abs(l))
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

#ifdef MULTIPRECISION
    const int MAX_ITERS = 200*DIGITS;
#else
    const int MAX_ITERS = 2000;
#endif

    for (; abs(curr-prev) > abs(curr*EPS) && i < MAX_ITERS; i++)
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