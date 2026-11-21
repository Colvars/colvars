#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "colvarproxy.h"
#include "colvarproxy_stub.h"
#include "colvar_UIestimator.h"

namespace {
struct record {
  std::vector<double> position;
  int count;
  std::vector<double> gradient;
};
int failures = 0;

// Files are created in CTest's binary working directory, never in the sources.
void check(std::string const &name, int dim, std::vector<record> const &rows,
           std::vector<int> const &counts, std::vector<double> const &gradients,
           std::string const &bad_count = "", std::string const &bad_grad = "",
           bool duplicate_file = false) {
  std::string prefix = "ui-test-" + name;
  {
    std::ofstream count(prefix + ".UI.count"), grad(prefix + ".UI.grad");
    count.precision(17);
    grad.precision(17);
    count << "# " << dim << '\n';
    grad << "# " << dim << '\n';
    for (int d = 0; d < dim; ++d) {
      // The loader takes the number of records from the count header.
      int n = d == 0 ? static_cast<int>(rows.size()) + (!bad_count.empty() || !bad_grad.empty()) : 1;
      count << "# 0 1 " << n << " 0\n";
      grad << "# 0 1 " << n << " 0\n";
    }
    for (auto const &row : rows) {
      for (double p : row.position) { count << p << ' '; grad << p << ' '; }
      count << row.count << '\n';
      for (double g : row.gradient) grad << g << ' ';
      grad << '\n';
    }
    count << bad_count;
    grad << bad_grad;
  }
  colvarproxy_stub proxy;
  std::vector<double> lower(dim, 0), upper(dim, 2), width(dim, 1), force(dim, 1);
  std::vector<std::string> inputs{prefix};
  if (duplicate_file) inputs.push_back(prefix);
  UIestimator::UIestimator estimator(proxy.cvmodule, lower, upper, width, force,
                                   prefix + "-out", true, inputs, 300);
  // Zero new samples: output exposes exactly the counts/gradients loaded.
  estimator.write_files();
  std::ifstream c(prefix + "-out.UI.count"), g(prefix + "-out.UI.grad");
  std::string line;
  auto skip_header = [&](std::ifstream &file) {
    for (int i = 0; i <= dim; ++i) std::getline(file, line);
  };
  skip_header(c); skip_header(g);
  bool ok = true;
  for (size_t bin = 0; bin < counts.size(); ++bin) {
    double ignored, value;
    int count = -1;
    for (int d = 0; d < dim; ++d) { c >> ignored; g >> ignored; }
    c >> count;
    ok = ok && bool(c) && count == counts[bin];
    for (int d = 0; d < dim; ++d) {
      g >> value;
      ok = ok && bool(g) && std::fabs(value - gradients[bin * dim + d]) < 1.e-5;
    }
  }
  std::cout << (ok ? "PASS " : "FAIL ") << name << '\n';
  failures += !ok;
}
}
int main() {
  check("valid-1d", 1, {{{0.5}, 2, {3}}, {{1.5}, 4, {7}}}, {2,4}, {3,7});
  check("valid-2d", 2, {{{0.5,0.5},2,{3,4}}, {{1.5,1.5},5,{6,7}}},
        {2,0,0,5}, {3,4,0,0,0,0,6,7});
  check("all-dimensions", 2, {{{-1,0.5},9,{99,99}}, {{0.5,3},9,{99,99}},
        {{3,0.5},9,{99,99}}, {{0.5,-1},9,{99,99}}, {{0.5,0.5},2,{3,4}}},
        {2,0,0,0}, {3,4,0,0,0,0,0,0});
  check("boundaries", 1, {{{0},2,{3}}, {{-0.5e-6},2,{5}}, {{-2e-6},9,{99}},
        {{2},9,{99}}, {{2-0.5e-6},9,{99}}, {{2-1.e-6},9,{99}},
        {{2-1.0001e-6},9,{99}}, {{2-2e-6},4,{7}}}, {4,4}, {4,7});
  check("upper-edge-2d", 2, {{{0.5,2-1.e-6},9,{99,99}},
        {{2-1.e-6,0.5},9,{99,99}}, {{0.5,0.5},2,{3,4}}},
        {2,0,0,0}, {3,4,0,0,0,0,0,0});
  check("zero-count", 1, {{{0.5},0,{99}}, {{1.5},2,{7}}}, {0,2}, {0,7});
  check("weighted-duplicates", 1, {{{0.5},2,{3}}, {{0.5},6,{7}}}, {8,0}, {6,0});
  // A second input prefix in 1D made position_temp[i] read past its end.
  check("multiple-files", 1, {{{0.5},2,{3}}}, {4,0}, {3,0}, "", "", true);
  check("truncated-count", 1, {{{0.5},2,{3}}}, {2,0}, {3,0}, "1.5", "1.5 99\n");
  check("truncated-grad", 1, {{{0.5},2,{3}}}, {2,0}, {3,0}, "1.5 8\n", "1.5");
  check("malformed-count", 1, {{{0.5},2,{3}}}, {2,0}, {3,0}, "oops 8\n", "1.5 99\n");
  check("malformed-grad", 1, {{{0.5},2,{3}}}, {2,0}, {3,0}, "1.5 8\n", "1.5 oops\n");
  check("truncated-2d", 2, {{{0.5,0.5},2,{3,4}}}, {2,0,0,0},
        {3,4,0,0,0,0,0,0}, "1.5 1.5 8\n", "1.5 1.5 9");
  return failures ? 1 : 0;
}
