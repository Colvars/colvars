#include "colvars_system.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>

// Exhaust the original triangular lattice inside the sphere defined by the
// unshifted displacement. This oracle uses no reduced-basis information.
double minimum_distance(cvm::rvector delta, cvm::rvector a, cvm::rvector b, cvm::rvector c)
{
  double best2 = delta.norm2();
  const double radius = std::sqrt(best2);
  const int first_z = static_cast<int>(std::ceil((delta.z - radius) / c.z));
  const int last_z = static_cast<int>(std::floor((delta.z + radius) / c.z));
  for (int z = first_z; z <= last_z; ++z) {
    const auto plane = delta - z*c;
    const double remaining2 = best2 - plane.z*plane.z;
    if (remaining2 <= 0) continue;
    const double dy = std::sqrt(remaining2);
    const int first_y = static_cast<int>(std::ceil((plane.y - dy) / b.y));
    const int last_y = static_cast<int>(std::floor((plane.y + dy) / b.y));
    for (int y = first_y; y <= last_y; ++y) {
      auto image = plane - y*b;
      image.x -= a.x * std::round(image.x / a.x);
      best2 = std::min(best2, image.norm2());
    }
  }
  return std::sqrt(best2);
}

int main()
{
  std::mt19937_64 rng(1763);
  std::uniform_real_distribution<double> uniform(-0.5, 0.5);
  double maximum_error = 0.0;
  for (int t = 0; t < 5000; ++t) {
    const cvm::rvector a(std::exp(7 * uniform(rng)), 0, 0);
    const cvm::rvector b(a.x * uniform(rng), std::exp(7 * uniform(rng)), 0);
    const cvm::rvector c(a.x * uniform(rng), b.y * uniform(rng), std::exp(7 * uniform(rng)));
    const auto delta = uniform(rng)*a + uniform(rng)*b + uniform(rng)*c;
    cvm::system_boundary_conditions cell;
    if (!cell.set_triclinic_minimum_image(a, b, c)) return 1;
    const double expected = minimum_distance(delta, a, b, c);
    const double actual = std::sqrt(cell.position_distance({}, delta).norm2());
    const double error = std::abs(actual - expected);
    if (!(error <= 1e-12)) {
      std::cerr << "Trial " << t << ", expected " << expected << ", actual " << actual
                << ", absolute error " << error << '\n';
      return 1;
    }
    maximum_error = std::max(maximum_error, error);
  }
  std::cout << "5000 triclinic displacements, maximum distance error " << maximum_error << " nm\n";
}
