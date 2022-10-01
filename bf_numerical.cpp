#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cmath>

#define N_BASIS 50

int main(int argc, char *argv[]) {
  using mat = Eigen::Matrix<double, N_BASIS, N_BASIS>;

  if(argc != 2) {
    std::cout << "Expect only one input. Aborted.\n";
    exit(0); 
  }
  auto const eps = std::atof(argv[1]);
  
  Eigen::Vector<double, N_BASIS> diagonal;
  for (int i = 0; i < N_BASIS; ++i) {
    diagonal[i] = 2 * i + 1;
  }

  mat h0 = diagonal.asDiagonal();

  mat perturbation;
  for (int i = 0; i < N_BASIS; ++i) {
    perturbation(i, i) = (6. * i * i + 6. * i + 3.) / 4.;
    if (i + 2 < N_BASIS)
      perturbation(i, i + 2) = perturbation(i + 2, i) =
          (i + 3. / 2) * std::sqrt((i + 1) * (i + 2));
    if (i + 4 < N_BASIS)
      perturbation(i, i + 4) = perturbation(i + 4, i) =
          std::sqrt((i + 1) * (i + 2) * (i + 3) * (i + 4)) / 4;
  }
 
  mat full_h = h0 - perturbation * eps;
  Eigen::SelfAdjointEigenSolver<mat> es;
  es.compute(full_h, Eigen::EigenvaluesOnly);
  std::vector<double> ev(N_BASIS, 0);

  for (int i = 0; i < N_BASIS; ++i) {
    ev[i] = full_h.eigenvalues()[i].real();
    if (full_h.eigenvalues()[i].imag() > 1e-10) {
        std::cout << "Bad eigenvalue" << full_h.eigenvalues()[i] << '\n';
        exit(0);
    }
  }
  std::sort(ev.begin(), ev.end());
  std::cout << std::setprecision(15);
  for (int i = 0; i < 10; ++i) {
    std::cout << ev[i] << '\n';
  }
}