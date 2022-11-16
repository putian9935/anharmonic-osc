#include <gmpxx.h>
#include <iostream>

#include "simple_poly.hpp"
#include <Eigen/Dense>
#include <iomanip>

#define N 50
#define M 50
#define COEFF_NUM N + M + 1

#include "pade_coef.hpp"

template <typename _Scalar, int deg>
Polynomial<_Scalar, deg - 1>
vector_to_poly(Eigen::Vector<_Scalar, deg> const &v) {
  Polynomial<_Scalar, deg - 1> ret;
  for (int i = 0; i < deg; ++i) {
    ret[i] = v[i];
  }
  return ret;
}

int main() {
  Eigen::Matrix<mpq_class, M, M> a_coef;
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < M; ++j) {
      a_coef(i, j) = coef[N + i - j];
    }
  }

  Eigen::Vector<mpq_class, M> rhs;
  for (int i = 0; i < M; ++i) {
    rhs(i) = -coef[N + i + 1];
  }
  Eigen::Vector<mpq_class, M> b = a_coef.fullPivLu().solve(rhs);

  Eigen::Vector<mpq_class, N + 1> a;
  for (int i = 0; i <= N; ++i) {
    a(i) = coef[i];
    for (int j = 0; j < i; ++j) {
      a(i) += coef[i - j - 1] * b(j);
    }
  }

  std::cout << std::setprecision(15); 
  
  auto poly_b = vector_to_poly(b);
  for (int num = -100; num <= 100; ++num) {

    mpq_class const eps(num, 10);
    mpq_class const res =
        vector_to_poly(a).evaluate(eps) /
        (poly_b * make_polynomial<mpq_class, 1>({mpq_class(0), 1}) +
         make_polynomial<mpq_class, 0>({mpq_class(1)}))
            .evaluate(eps);
    std::cout << res.get_d() << '\n';
  }
}