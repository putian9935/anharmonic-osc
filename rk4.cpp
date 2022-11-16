#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>

#ifdef BVC_ASYM
#define BVC bvc_t::ASYM
#elif defined(BVC_SYM)
#define BVC bvc_t::SYM 
#endif 

using vec = Eigen::Vector4d;
using mat = Eigen::Matrix4<double>;

mat c, r2, r4;
double const step_size = 2e-4; // for negative
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
  for (double r = 0; r < 40; r += step_size) { // negative
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
  if (!(std::signbit(left) ^ std::signbit(right))) {
    std::cout << "Wrong boundary value. Aborted.\n";
    std::cout << "Left: " << left << '\n';
    std::cout << "Right: " << right << '\n';
    exit(-1);
  }

  double mid;
  double energy_mid;
  while (energy_right > energy_left + energy_res) {
    energy_mid = (energy_left + energy_right) / 2.;
    mid = evaluate_scalar<bvc>(energy_mid);
    if (std::signbit(left) ^ std::signbit(mid)) {
      right = mid;
      energy_right = energy_mid;
    } else if (std::signbit(right) ^ std::signbit(mid)) {
      left = mid;
      energy_left = energy_mid;
    } else {
      std::cout << "Failed around " << energy_mid << " , " << mid << '\n';
      exit(-1);
    }
  }
  return energy_mid;
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    std::cout << "Expect 3 arguments. Aborted.\n";
    exit(-1);
  }
  auto const eps = std::atof(argv[1]), energy_left = std::atof(argv[2]),
             energy_right = std::atof(argv[3]);

  r2(2, 0) = -1. / 2;
  r2(2, 1) = +std::sqrt(3) / 2;
  r2(3, 0) = -std::sqrt(3) / 2;
  r2(3, 1) = -1. / 2;

  r4(2, 0) = eps;
  r4(3, 1) = eps;

  std::cout << std::setprecision(16);
  std::cout << find_energy_bisect<BVC>(energy_left, energy_right)
            << '\n';
  return 0;
}