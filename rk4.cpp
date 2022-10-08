#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>

using vec = Eigen::Vector4d;
using mat = Eigen::Matrix4<double>;

mat c, r2, r4;
double const step_size = 5e-5;
double const energy_res = 5e-7;

vec f(double const r, vec const &y) {
  return (c + r * r * r2 + r * r * r * r * r4) * y;
}

typedef enum { SYM, ASYM } bvc_t;

template <bvc_t bvc = SYM> double evaluate_scalar(double const energy) {
  c(0, 2) = c(1, 3) = 1;
  c(2, 0) = -1. / 2 * energy;
  c(2, 1) = -std::sqrt(3) / 2 * energy;
  c(3, 0) = std::sqrt(3) / 2 * energy;
  c(3, 1) = -1. / 2 * energy;

  vec y;
  if constexpr (bvc == SYM)
    y << 1, 0, 0, 0;
  else
    y << 0, 0, 1, 0;
  for (double r = 0; r < 10; r += step_size) {
    vec k1 = step_size * f(r, y);
    vec k2 = step_size * f(r + step_size / 3., y + k1 / 3.);
    vec k3 = step_size * f(r + 2. * step_size / 3., y - k1 / 3. + k2);
    vec k4 = step_size * f(r + step_size, y + k1 - k2 + k3);
    y += (k1 + 3 * k2 + 3 * k3 + k4) / 8.;
  }
  return y[0];
}

template <bvc_t bvc = SYM>
double find_energy_bisect(double energy_left, double energy_right) {
  double left = evaluate_scalar<bvc>(energy_left),
         right = evaluate_scalar<bvc>(energy_right);
  if (left * right > 0) {
    std::cout << "Wrong boundary value. Aborted.\n";
    exit(-1);
  }

  double mid;
  double energy_mid;
  while (energy_right > energy_left + energy_res) {
    energy_mid = (energy_left + energy_right) / 2.;
    mid = evaluate_scalar<bvc>(energy_mid);
    if (left * mid < 0) {
      right = mid;
      energy_right = energy_mid;
    } else if (right * mid < 0) {
      left = mid;
      energy_left = energy_mid;
    } else {
      std::cout << "Failed around " << energy_mid << " , " << mid << '\n';
      exit(-1);
    }
  }
  return energy_mid;
}

int main() {
  double const eps = -0.1;

  r2(2, 0) = -1. / 2;
  r2(2, 1) = +std::sqrt(3) / 2;
  r2(3, 0) = -std::sqrt(3) / 2;
  r2(3, 1) = -1. / 2;

  r4(2, 0) = eps;
  r4(3, 1) = eps;

  double energy_left = 3, energy_right = 3.4;

  std::cout << std::setprecision(16);
  std::cout << find_energy_bisect<bvc_t::ASYM>(energy_left, energy_right)
            << '\n';
  return 0;
}