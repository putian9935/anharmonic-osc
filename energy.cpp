#include "simple_poly.hpp"
#include <gmpxx.h>
#include <iostream>

mpq_class const gamma_lut[] = {
    1,
    1,
    mpq_class(1, 2),
    1,
    mpq_class(3, 4),
    2,
    mpq_class(15, 8),
    6,
    mpq_class(105, 16),
    24,
    mpq_class(945, 32),
    120,
    mpq_class(10395, 64),
    720,
};

template <int deg> mpq_class integrate(Polynomial<mpq_class, deg> p) {
  mpq_class ret = 0;
  for (int i = 0; i <= deg; ++i) {
    ret += gamma_lut[i] * p[i];
  }
  return ret;
}

#include <Eigen/Dense>
#include <Eigen/Eigen>

template <typename _Scalar, int deg>
Eigen::Vector<_Scalar, deg + 1>
poly_to_vector(Polynomial<_Scalar, deg> const &p) {
  return Eigen::Vector<_Scalar, deg + 1>(p.coef.data());
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

// template <typename _Scalar, int n> 
// Polynomial<_Scalar, 4 *n> 
// update() {

//   auto p1 = w * y0 + mpq_class(-e1) * y0;
//   Eigen::Matrix<mpq_class, 4 * n + 2 * m - 1, 4 * n + 2 * m - 2> q1_coef;
//   for (int i = 0; i < 4 * n + 2 * m - 3; ++i)
//     q1_coef(i, i + 1) = i + 1;
//   for (int i = 1; i < 4 * n + 2 * m - 2; ++i)
//     q1_coef(i, i - 1) = -2;
 

//   auto b1 = y0 * p1;
//   Eigen::Vector<mpq_class, 4 *n + 2 *m - 2> q1_vec =
//       q1_coef.fullPivLu().solve(poly_to_vector(b1));
//   auto q1 = vector_to_poly(q1_vec);

//   Eigen::Matrix<mpq_class, 4 * n + 2 * m - 2, 4 * n> r1_coef;
//   for (int i = 0; i < 4 * n; ++i) {
//     for (int l = 0; l < std::min(2 * m + 1, 4 * n + 2 * m - 2 - i); ++l) {
//       r1_coef(l + i, i) = ysq[l] * (i + 1);
//     }
//   }


//   Eigen::Vector<mpq_class, 4 * n> r1_vec =
//       r1_coef.fullPivLu().solve(poly_to_vector(q1));
//   auto r1 = vector_to_poly(r1_vec);

//   auto y1 = y0 * r1 * make_polynomial<mpq_class, 1>({0, 1});
//   mpq_class e1 = integrate(y0 * (w * y1 + mpq_class(-e1) * y1)) /
//                    integrate(y0 * y0);
//   std::cout << e1 << '\n';
// }

int main() {
  int const m = 2, n = 1;
  auto w = make_polynomial<mpq_class, 4>({0, 0, 0, 0, -1});
  auto y0 = make_polynomial<mpq_class, 1>({0, 2});
  auto ysq = y0 * y0;

  mpq_class e1 = integrate(y0 * w * y0) / integrate(y0 * y0);
  std::cout << e1 << '\n';

  auto p1 = w * y0 + mpq_class(-e1) * y0;
  Eigen::Matrix<mpq_class, 4 * n + 2 * m - 1, 4 * n + 2 * m - 2> q1_coef;
  for (int i = 0; i < 4 * n + 2 * m - 3; ++i)
    q1_coef(i, i + 1) = i + 1;
  for (int i = 1; i < 4 * n + 2 * m - 2; ++i)
    q1_coef(i, i - 1) = -2;
 

  auto b1 = y0 * p1;
  Eigen::Vector<mpq_class, 4 *n + 2 *m - 2> q1_vec =
      q1_coef.fullPivLu().solve(poly_to_vector(b1));
  auto q1 = vector_to_poly(q1_vec);

  Eigen::Matrix<mpq_class, 4 * n + 2 * m - 2, 4 * n> r1_coef;
  for (int i = 0; i < 4 * n; ++i) {
    for (int l = 0; l < std::min(2 * m + 1, 4 * n + 2 * m - 2 - i); ++l) {
      r1_coef(l + i, i) = ysq[l] * (i + 1);
    }
  }


  std::cout << q1 << "???\n";

  // Eigen::Vector<mpq_class, 4 * n> r1_vec =
  //     r1_coef.fullPivLu().solve(poly_to_vector(q1));
  // auto r1 = vector_to_poly(r1_vec);


  // std::cout << "??\n";
  // auto y1 = y0 * r1 * make_polynomial<mpq_class, 1>({0, 1});
  // mpq_class e2 = integrate(y0 * (w * y1 + mpq_class(-e1) * y1)) /
  //                  integrate(y0 * y0);
  // std::cout << e2 << '\n';
/*==========*/
  // int const nn = 2; 
  // auto p2 = w * y1 + mpq_class(-e2) * y0 + mpq_class(-e1) * y1; 

  // Eigen::Matrix<mpq_class, 4 * nn + 2 * m - 1, 4 * nn + 2 * m - 2> q2_coef;
  // for (int i = 0; i < 4 * nn + 2 * m - 3; ++i)
  //   q2_coef(i, i + 1) = i + 1;
  // for (int i = 1; i < 4 * nn + 2 * m - 2; ++i)
  //   q2_coef(i, i - 1) = -2;

  // auto b2 = y0 * p2;
  // Eigen::Vector<mpq_class, 4 * nn + 2 *m - 2> q2_vec =
  //     q2_coef.fullPivLu().solve(poly_to_vector(b2));
  // auto q2 = vector_to_poly(q2_vec);

  // Eigen::Matrix<mpq_class, 4 * nn + 2 * m - 2, 4 * nn> r2_coef;
  // for (int i = 0; i < 4 * nn; ++i) {
  //   for (int l = 0; l < std::min(2 * m + 1, 4 * nn + 2 * m - 2 - i); ++l) {
  //     r2_coef(l + i, i) = ysq[l] * (i + 1);
  //   }
  // }


  // Eigen::Vector<mpq_class, 4 * nn> r2_vec =
  //     r2_coef.fullPivLu().solve(poly_to_vector(q2));
  // auto r2 = vector_to_poly(r2_vec);
  // std::cout << r2 << '\n';
  // auto y2 = y0 * r2 * make_polynomial<mpq_class, 1>({0, 1});
  // auto e3 = integrate(y0 * (w * y2 + mpq_class(-e1) * y2 + mpq_class(-e2) * y1)) /
  //                  integrate(y0 * y0);
  // std::cout << e3 << '\n';


}