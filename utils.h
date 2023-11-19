#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <iostream>
#include <vector>
#include <cmath>
#include <numbers>
#include <complex>
#include "polynomial.h"
#include "config.h"


#ifdef MULTIPRECISION
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>
#endif


#ifdef MULTIPRECISION
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

std::vector<val_t> find_roots(const Polynomial<val_t>&, val_t, val_t);
val_t root_binary_search(const Polynomial<val_t>& p, val_t l, val_t r);

struct sturm_sequence
{
    std::vector<Polynomial<val_t>> seq;
    sturm_sequence(Polynomial<val_t> p);
    int operator()(val_t val) const;
};

#endif 