#include "simple_poly.hpp"
#include <gmpxx.h>
#include <iostream>

int main() {
  auto y0 = make_polynomial<mpq_class, 3>({0,-12,0,8});
  
  std:: cout << make_polynomial<mpq_class, 2>({-12,0,8}) + mpq_class(-2) * y0 << '\n';

  std::cout << y0.evaluate(2) << '\n';
}