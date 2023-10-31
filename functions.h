#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include <functional>
#include <algorithm>
#include <stdio.h>
#include <cmath>
#include <ctime>

#include "matrix.h"
#include "polynomial.h"
#include "utils.h"

val_t fp(val_t);
val_t f(val_t);
val_t p(val_t);

val_t integrate_riemann(const std::function<val_t(val_t)>& f, val_t l, val_t r, size_t);

struct Solve{val_t res; val_t sum;};
Solve integrate_newton(const std::function<val_t(val_t)>& f, val_t l, val_t r, size_t);

Solve integrate_gauss(const std::function<val_t(val_t)>& f, val_t l, val_t r, size_t);


val_t integrate(const std::function<val_t(val_t)>& f, val_t l, val_t r, val_t h, bool gauss = 1, size_t k = 4);

#endif

