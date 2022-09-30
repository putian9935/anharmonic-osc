#include "simple_poly.hpp"
#include <gmpxx.h>
#include <iostream>

#define MAX_ORDER 1000
mpq_class gamma_lut[MAX_ORDER] = {mpq_class(1), mpq_class(1)};

void init_gamma_lut() {
  for (int i = 2; i < MAX_ORDER; ++i)
    gamma_lut[i] = mpq_class((i - 1), 2) * gamma_lut[i - 2];
}

template <int deg> mpq_class integrate(Polynomial<mpq_class, deg> p) {
  mpq_class ret = 0;
  for (int i = 0; i <= deg; ++i) {
    ret += gamma_lut[i] * p[i];
  }
  return ret;
}

#include <Eigen/Dense>
#include <Eigen/Eigen>

auto const w = make_polynomial<mpq_class, 4>({0, 0, 0, 0, -1});
int const m = 2;
auto const my_y0 = make_polynomial<mpq_class, m>({-2, 0, 4});

template <typename _Scalar, int deg>
Eigen::Vector<_Scalar, deg + 1>
poly_to_vector(Polynomial<_Scalar, deg> const &p) {
  return Eigen::Vector<_Scalar, deg + 1>(p.coef.data());
}

template <typename _Scalar, int deg>
Eigen::Vector<_Scalar, deg + 3>
poly_to_vector_plus2(Polynomial<_Scalar, deg> const &p) {
  Eigen::Vector<_Scalar, deg + 3> ret;
  for (int i = 0; i < p.coef.size(); ++i)
    ret[i] = p.coef[i];
  return ret;
}

template <typename _Scalar, int deg>
Eigen::Vector<_Scalar, deg + 4>
poly_to_vector_plus3(Polynomial<_Scalar, deg> const &p) {
  Eigen::Vector<_Scalar, deg + 4> ret;
  for (int i = 0; i < p.coef.size(); ++i)
    ret[i] = p.coef[i];
  return ret;
}

template <typename _Scalar, int deg>
Polynomial<_Scalar, deg - 1>
vector_to_poly(Eigen::Vector<_Scalar, deg> const &v) {
  Polynomial<_Scalar, deg - 1> ret;
  for (int i = 0; i < deg; ++i) {
    ret[i] = v[i];
  }
  return ret;
}

template <typename _Scalar, int n, int m>
Polynomial<_Scalar, 4 * n + m> update(Polynomial<_Scalar, 4 * n + 2 * m> b) {
  // Eigen::Matrix<mpq_class, 4 * n + 2 * m + 3, 4 * n + m + 1> b_coef;
  Eigen::Matrix<mpq_class, -1, -1> b_coef;
  b_coef.resize(4 * n + 2 * m + 3, 4 * n + m + 1);
  for (int k = 0; k <= m; ++k) {
    for (int l = 0; l <= 4 * n + m; ++l) {
      b_coef(k + l, l) += -1 * my_y0[k];
    }
    for (int l = 2; l <= 4 * n + m + 2; ++l) {
      b_coef(k + l, l - 2) += my_y0[k];
    }
    for (int l = 1; l <= 4 * n + m; ++l) {
      b_coef(k + l, l) += -2 * l * my_y0[k];
    }
    for (int l = 0; l <= 4 * n + m - 2; ++l) {
      b_coef(k + l, l + 2) += (l + 1) * (l + 2) * my_y0[k];
    }
  }
  for (int k = 0; k <= 4 * n + m; ++k) {
    for (int l = 0; l <= (m); ++l) {
      b_coef(k + l, k) += 1 * my_y0[l];
    }
    for (int l = 2; l <= m + 2; ++l) {
      b_coef(k + l, k) += -my_y0[l - 2];
    }
    for (int l = 1; l <= m; ++l) {
      b_coef(k + l, k) += 2 * l * my_y0[l];
    }
    for (int l = 0; l <= m - 2; ++l) {
      b_coef(k + l, k) += -(l + 1) * (l + 2) * my_y0[l + 2];
    }
  }

  Eigen::Vector<mpq_class, 4 *n + m + 1> ret_vec =
      b_coef.fullPivLu().solve(poly_to_vector_plus2(b));

  return vector_to_poly(ret_vec);
}

int main() {
  init_gamma_lut();

  mpq_class e1 = integrate(my_y0 * w * my_y0) / integrate(my_y0 * my_y0);
  std::cout << e1 << '\n';

  auto y1 =
      update<mpq_class, 1, m>(my_y0 * (w * my_y0 + mpq_class(-e1) * my_y0));
  mpq_class e2 = integrate(my_y0 * (w * y1 + mpq_class(-e1) * y1)) /
                 integrate(my_y0 * my_y0);
  std::cout << e2 << '\n';

  auto y2 = update<mpq_class, 2, m>(
      my_y0 * (w * y1 + mpq_class(-e1) * y1 + mpq_class(-e2) * my_y0));
  mpq_class e3 =
      integrate(my_y0 * (w * y2 + mpq_class(-e1) * y2 + mpq_class(-e2) * y1)) /
      integrate(my_y0 * my_y0);
  std::cout << e3 << '\n';

    auto y3 = update<mpq_class, 3, m>(my_y0 * (w * y2 + mpq_class(-e1) * y2+mpq_class(-e2) * y1+ mpq_class(-e3) * my_y0));
    mpq_class e4 = integrate(my_y0 * (w * y3 + mpq_class(-e1) * y3+mpq_class(-e2) * y2+mpq_class(-e3) * y1)) / integrate(my_y0 * my_y0);
    std::cout << e4 << '\n';
    
  
    auto y4 = update<mpq_class, 4, m>(my_y0 * (w * y3 + mpq_class(-e1) * y3+mpq_class(-e2) * y2+mpq_class(-e3) * y1+ mpq_class(-e4) * my_y0));
    mpq_class e5 = integrate(my_y0 * (w * y4 + mpq_class(-e1) * y4+mpq_class(-e2) * y3+mpq_class(-e3) * y2+mpq_class(-e4) * y1)) / integrate(my_y0 * my_y0);
    std::cout << e5 << '\n';
    
  
    auto y5 = update<mpq_class, 5, m>(my_y0 * (w * y4 + mpq_class(-e1) * y4+mpq_class(-e2) * y3+mpq_class(-e3) * y2+mpq_class(-e4) * y1+ mpq_class(-e5) * my_y0));
    mpq_class e6 = integrate(my_y0 * (w * y5 + mpq_class(-e1) * y5+mpq_class(-e2) * y4+mpq_class(-e3) * y3+mpq_class(-e4) * y2+mpq_class(-e5) * y1)) / integrate(my_y0 * my_y0);
    std::cout << e6 << '\n';
    
  
    auto y6 = update<mpq_class, 6, m>(my_y0 * (w * y5 + mpq_class(-e1) * y5+mpq_class(-e2) * y4+mpq_class(-e3) * y3+mpq_class(-e4) * y2+mpq_class(-e5) * y1+ mpq_class(-e6) * my_y0));
    mpq_class e7 = integrate(my_y0 * (w * y6 + mpq_class(-e1) * y6+mpq_class(-e2) * y5+mpq_class(-e3) * y4+mpq_class(-e4) * y3+mpq_class(-e5) * y2+mpq_class(-e6) * y1)) / integrate(my_y0 * my_y0);
    std::cout << e7 << '\n';
    
  
    auto y7 = update<mpq_class, 7, m>(my_y0 * (w * y6 + mpq_class(-e1) * y6+mpq_class(-e2) * y5+mpq_class(-e3) * y4+mpq_class(-e4) * y3+mpq_class(-e5) * y2+mpq_class(-e6) * y1+ mpq_class(-e7) * my_y0));
    mpq_class e8 = integrate(my_y0 * (w * y7 + mpq_class(-e1) * y7+mpq_class(-e2) * y6+mpq_class(-e3) * y5+mpq_class(-e4) * y4+mpq_class(-e5) * y3+mpq_class(-e6) * y2+mpq_class(-e7) * y1)) / integrate(my_y0 * my_y0);
    std::cout << e8 << '\n';
    
  
    auto y8 = update<mpq_class, 8, m>(my_y0 * (w * y7 + mpq_class(-e1) * y7+mpq_class(-e2) * y6+mpq_class(-e3) * y5+mpq_class(-e4) * y4+mpq_class(-e5) * y3+mpq_class(-e6) * y2+mpq_class(-e7) * y1+ mpq_class(-e8) * my_y0));
    mpq_class e9 = integrate(my_y0 * (w * y8 + mpq_class(-e1) * y8+mpq_class(-e2) * y7+mpq_class(-e3) * y6+mpq_class(-e4) * y5+mpq_class(-e5) * y4+mpq_class(-e6) * y3+mpq_class(-e7) * y2+mpq_class(-e8) * y1)) / integrate(my_y0 * my_y0);
    std::cout << e9 << '\n';
    
  
    auto y9 = update<mpq_class, 9, m>(my_y0 * (w * y8 + mpq_class(-e1) * y8+mpq_class(-e2) * y7+mpq_class(-e3) * y6+mpq_class(-e4) * y5+mpq_class(-e5) * y4+mpq_class(-e6) * y3+mpq_class(-e7) * y2+mpq_class(-e8) * y1+ mpq_class(-e9) * my_y0));
    mpq_class e10 = integrate(my_y0 * (w * y9 + mpq_class(-e1) * y9+mpq_class(-e2) * y8+mpq_class(-e3) * y7+mpq_class(-e4) * y6+mpq_class(-e5) * y5+mpq_class(-e6) * y4+mpq_class(-e7) * y3+mpq_class(-e8) * y2+mpq_class(-e9) * y1)) / integrate(my_y0 * my_y0);
    std::cout << e10 << '\n';
}