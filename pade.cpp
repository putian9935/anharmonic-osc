#include <gmpxx.h>
#include <iostream>

#include <Eigen/Dense>

#define N 4
#define M 4
#define COEFF_NUM N + M + 1

// Pade(4,4) approximant to e^x
mpq_class coef[COEFF_NUM] = {
    mpq_class("1"),    mpq_class("1"),        mpq_class(1, 2),
    mpq_class(1, 6),   mpq_class(1, 24),      mpq_class(1, 120),
    mpq_class(1, 720), mpq_class(1, 720 * 7), mpq_class(1, 720 * 7 * 8)};

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
  std::cout << "Denominator: " << 1 << ' ' << b.transpose() << '\n';
  Eigen::Vector<mpq_class, N + 1> a;
  for (int i = 0; i <= N; ++i) {
    a(i) = coef[i];
    for (int j = 0; j < i; ++j) {
      a(i) += coef[i - j - 1] * b(j);
    }
  }
  std::cout << "Numerator: " << a.transpose() << '\n';
}