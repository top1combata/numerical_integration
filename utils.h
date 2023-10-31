#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <iostream>
#include <vector>
#include <cmath>
#include <numbers>
#include <complex>

// #define MULTIPRECISION

#ifdef MULTIPRECISION

#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>

#endif

#include "polynomial.h"

#ifdef MULTIPRECISION
#define DIGITS 100
using val_t = boost::multiprecision::number<boost::multiprecision::backends::gmp_float<DIGITS>>;
#else
using val_t = long double;
#endif

#ifdef MULTIPRECISION
const val_t PI    = boost::math::constants::pi<val_t>();
const val_t E     = boost::math::constants::e<val_t>();
const val_t EPS   = pow(val_t(10),-DIGITS);
#else
const val_t PI    = std::numbers::pi_v<val_t>;
const val_t E     = std::numbers::e_v<val_t>;
const val_t EPS   = 1e-17;
#endif

const val_t A     = val_t(13)/10;
const val_t B     = val_t(11)/5;
const val_t ALPHA = val_t(0);
const val_t BETA  = val_t(5)/6;


val_t my_sin(val_t);
val_t my_cos(val_t);
val_t my_exp(val_t);
val_t my_pow(val_t,val_t);
val_t my_log(val_t);
val_t my_abs(val_t);
val_t my_sqrt(val_t);

std::vector<val_t> roots_cubic(const Polynomial<val_t>& poly);

std::vector<val_t> find_roots(const Polynomial<val_t>&, val_t, val_t);
val_t root_binary_search(const Polynomial<val_t>& p, val_t l, val_t r);

struct sturm_sequence
{
    std::vector<Polynomial<val_t>> seq;
    sturm_sequence(Polynomial<val_t> p);
    int operator()(val_t val) const;
};

#endif 