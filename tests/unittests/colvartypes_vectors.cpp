// -*- c++ -*-

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "colvarmodule.h"
#include "colvartypes.h"


// Simple test assertion helpers
static int test_failures = 0;

#define CHECK(cond, msg)                                                       \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::cerr << "FAIL: " << msg << std::endl;                              \
      test_failures++;                                                         \
    }                                                                          \
  } while (0)

#define CHECK_NEAR(a, b, eps, msg)                                             \
  CHECK(std::fabs((a) - (b)) < (eps), msg << ": " << (a) << " != " << (b))


// Test cvm::rvector construction and access
static void test_rvector_construction()
{
  // Default constructor resets to zero
  cvm::rvector v0;
  CHECK(v0.x == 0.0 && v0.y == 0.0 && v0.z == 0.0,
        "rvector default constructor should zero all components");

  // Parameterized constructor
  cvm::rvector v1(1.0, 2.0, 3.0);
  CHECK(v1.x == 1.0 && v1.y == 2.0 && v1.z == 3.0,
        "rvector parameterized constructor");

  // Scalar constructor
  cvm::rvector v2(5.0);
  CHECK(v2.x == 5.0 && v2.y == 5.0 && v2.z == 5.0,
        "rvector scalar constructor");

  // Index access
  cvm::rvector v3(4.0, 5.0, 6.0);
  CHECK(v3[0] == 4.0 && v3[1] == 5.0 && v3[2] == 6.0,
        "rvector index access");

  // set() method
  cvm::rvector v4;
  v4.set(7.0, 8.0, 9.0);
  CHECK(v4.x == 7.0 && v4.y == 8.0 && v4.z == 9.0,
        "rvector set(x,y,z)");

  v4.set(3.0);
  CHECK(v4.x == 3.0 && v4.y == 3.0 && v4.z == 3.0,
        "rvector set(scalar)");

  // reset()
  cvm::rvector v5(1.0, 2.0, 3.0);
  v5.reset();
  CHECK(v5.x == 0.0 && v5.y == 0.0 && v5.z == 0.0,
        "rvector reset()");
}


// Test cvm::rvector arithmetic
static void test_rvector_arithmetic()
{
  cvm::rvector a(1.0, 2.0, 3.0);
  cvm::rvector b(4.0, 5.0, 6.0);

  // Addition
  cvm::rvector sum = a + b;
  CHECK(sum.x == 5.0 && sum.y == 7.0 && sum.z == 9.0,
        "rvector addition");

  // Subtraction
  cvm::rvector diff = b - a;
  CHECK(diff.x == 3.0 && diff.y == 3.0 && diff.z == 3.0,
        "rvector subtraction");

  // Negation
  cvm::rvector neg = -a;
  CHECK(neg.x == -1.0 && neg.y == -2.0 && neg.z == -3.0,
        "rvector negation");

  // Scalar multiplication
  cvm::rvector scaled = 2.0 * a;
  CHECK(scaled.x == 2.0 && scaled.y == 4.0 && scaled.z == 6.0,
        "rvector scalar multiplication (scalar*vec)");

  cvm::rvector scaled2 = a * 3.0;
  CHECK(scaled2.x == 3.0 && scaled2.y == 6.0 && scaled2.z == 9.0,
        "rvector scalar multiplication (vec*scalar)");

  // Scalar division
  cvm::rvector divided = a / 2.0;
  CHECK_NEAR(divided.x, 0.5, 1e-14, "rvector division x");
  CHECK_NEAR(divided.y, 1.0, 1e-14, "rvector division y");
  CHECK_NEAR(divided.z, 1.5, 1e-14, "rvector division z");

  // In-place operators
  cvm::rvector c = a;
  c += b;
  CHECK(c.x == 5.0 && c.y == 7.0 && c.z == 9.0,
        "rvector operator+=");

  cvm::rvector d = b;
  d -= a;
  CHECK(d.x == 3.0 && d.y == 3.0 && d.z == 3.0,
        "rvector operator-=");

  cvm::rvector e = a;
  e *= 2.0;
  CHECK(e.x == 2.0 && e.y == 4.0 && e.z == 6.0,
        "rvector operator*=");

  cvm::rvector f = a;
  f /= 2.0;
  CHECK_NEAR(f.x, 0.5, 1e-14, "rvector operator/= x");
  CHECK_NEAR(f.y, 1.0, 1e-14, "rvector operator/= y");
  CHECK_NEAR(f.z, 1.5, 1e-14, "rvector operator/= z");
}


// Test cvm::rvector dot product and cross product
static void test_rvector_products()
{
  cvm::rvector a(1.0, 0.0, 0.0);
  cvm::rvector b(0.0, 1.0, 0.0);
  cvm::rvector c(0.0, 0.0, 1.0);

  // Dot product: unit vectors
  CHECK_NEAR(a * a, 1.0, 1e-14, "rvector dot product: (1,0,0)*(1,0,0)");
  CHECK_NEAR(a * b, 0.0, 1e-14, "rvector dot product: (1,0,0)*(0,1,0)");
  CHECK_NEAR(a * c, 0.0, 1e-14, "rvector dot product: (1,0,0)*(0,0,1)");

  cvm::rvector v1(3.0, 4.0, 0.0);
  cvm::rvector v2(4.0, 3.0, 0.0);
  CHECK_NEAR(v1 * v2, 24.0, 1e-14, "rvector dot product (3,4,0)*(4,3,0)");

  // Cross product: basic cross products of unit vectors
  cvm::rvector axb = cvm::rvector::outer(a, b);
  CHECK_NEAR(axb.x, 0.0, 1e-14, "rvector cross product x*y = z: x");
  CHECK_NEAR(axb.y, 0.0, 1e-14, "rvector cross product x*y = z: y");
  CHECK_NEAR(axb.z, 1.0, 1e-14, "rvector cross product x*y = z: z");

  cvm::rvector bxa = cvm::rvector::outer(b, a);
  CHECK_NEAR(bxa.z, -1.0, 1e-14, "rvector cross product y*x = -z");

  // Cross product: general case
  cvm::rvector p(1.0, 2.0, 3.0);
  cvm::rvector q(4.0, 5.0, 6.0);
  cvm::rvector pxq = cvm::rvector::outer(p, q);
  // (2*6-5*3, -(1*6-4*3), 1*5-4*2) = (-3, 6, -3)
  CHECK_NEAR(pxq.x, -3.0, 1e-14, "rvector cross product general: x");
  CHECK_NEAR(pxq.y, 6.0, 1e-14, "rvector cross product general: y");
  CHECK_NEAR(pxq.z, -3.0, 1e-14, "rvector cross product general: z");
}


// Test cvm::rvector norm and unit vector
static void test_rvector_norm()
{
  cvm::rvector v1(3.0, 4.0, 0.0);
  CHECK_NEAR(v1.norm2(), 25.0, 1e-14, "rvector norm2 (3,4,0)");
  CHECK_NEAR(v1.norm(), 5.0, 1e-14, "rvector norm (3,4,0)");

  cvm::rvector unit_v1 = v1.unit();
  CHECK_NEAR(unit_v1.x, 0.6, 1e-14, "rvector unit x");
  CHECK_NEAR(unit_v1.y, 0.8, 1e-14, "rvector unit y");
  CHECK_NEAR(unit_v1.z, 0.0, 1e-14, "rvector unit z");
  CHECK_NEAR(unit_v1.norm(), 1.0, 1e-14, "rvector unit norm should be 1");

  // Zero vector - unit() returns (1,0,0) by convention
  cvm::rvector v0;
  cvm::rvector unit_v0 = v0.unit();
  CHECK_NEAR(unit_v0.x, 1.0, 1e-14, "rvector unit of zero vec: x");
  CHECK_NEAR(unit_v0.y, 0.0, 1e-14, "rvector unit of zero vec: y");
  CHECK_NEAR(unit_v0.z, 0.0, 1e-14, "rvector unit of zero vec: z");

  // General 3D case
  cvm::rvector v2(1.0, 1.0, 1.0);
  CHECK_NEAR(v2.norm2(), 3.0, 1e-14, "rvector norm2 (1,1,1)");
  CHECK_NEAR(v2.norm(), std::sqrt(3.0), 1e-14, "rvector norm (1,1,1)");
}


// Test cvm::vector1d<double> construction and basic operations
static void test_vector1d_construction()
{
  // Default constructor
  cvm::vector1d<double> v0(5);
  CHECK(v0.size() == 5, "vector1d default constructor size");
  for (size_t i = 0; i < v0.size(); i++) {
    CHECK(v0[i] == 0.0, "vector1d default constructor zeros");
  }

  // Constructor from C array
  double arr[] = {1.0, 2.0, 3.0, 4.0};
  cvm::vector1d<double> v1(4, arr);
  CHECK(v1.size() == 4, "vector1d from C-array size");
  CHECK(v1[0] == 1.0 && v1[1] == 2.0 && v1[2] == 3.0 && v1[3] == 4.0,
        "vector1d from C-array values");

  // Copy constructor
  cvm::vector1d<double> v2(v1);
  CHECK(v2.size() == 4, "vector1d copy constructor size");
  CHECK(v2[0] == 1.0 && v2[1] == 2.0, "vector1d copy constructor values");

  // reset()
  cvm::vector1d<double> v3(v1);
  v3.reset();
  for (size_t i = 0; i < v3.size(); i++) {
    CHECK(v3[i] == 0.0, "vector1d reset()");
  }

  // resize()
  cvm::vector1d<double> v4(3);
  v4.resize(6);
  CHECK(v4.size() == 6, "vector1d resize()");

  // clear()
  cvm::vector1d<double> v5(5);
  v5.clear();
  CHECK(v5.size() == 0, "vector1d clear()");
}


// Test cvm::vector1d arithmetic
static void test_vector1d_arithmetic()
{
  double arr1[] = {1.0, 2.0, 3.0, 4.0};
  double arr2[] = {5.0, 6.0, 7.0, 8.0};
  cvm::vector1d<double> a(4, arr1);
  cvm::vector1d<double> b(4, arr2);

  // Addition
  cvm::vector1d<double> sum = a + b;
  CHECK(sum[0] == 6.0 && sum[1] == 8.0 && sum[2] == 10.0 && sum[3] == 12.0,
        "vector1d addition");

  // Subtraction
  cvm::vector1d<double> diff = b - a;
  CHECK(diff[0] == 4.0 && diff[1] == 4.0 && diff[2] == 4.0 && diff[3] == 4.0,
        "vector1d subtraction");

  // Scalar multiplication
  cvm::vector1d<double> scaled = a * 3.0;
  CHECK(scaled[0] == 3.0 && scaled[1] == 6.0 && scaled[2] == 9.0 && scaled[3] == 12.0,
        "vector1d scalar multiplication");

  cvm::vector1d<double> scaled2 = 2.0 * a;
  CHECK(scaled2[0] == 2.0 && scaled2[1] == 4.0,
        "vector1d scalar multiplication (reverse)");

  // Scalar division
  cvm::vector1d<double> divided = b / 2.0;
  CHECK_NEAR(divided[0], 2.5, 1e-14, "vector1d division [0]");
  CHECK_NEAR(divided[1], 3.0, 1e-14, "vector1d division [1]");

  // In-place operators
  cvm::vector1d<double> c(4, arr1);
  c += b;
  CHECK(c[0] == 6.0 && c[1] == 8.0, "vector1d operator+=");

  cvm::vector1d<double> d(4, arr2);
  d -= a;
  CHECK(d[0] == 4.0 && d[1] == 4.0, "vector1d operator-=");

  cvm::vector1d<double> e(4, arr1);
  e *= 2.0;
  CHECK(e[0] == 2.0 && e[1] == 4.0, "vector1d operator*=");

  cvm::vector1d<double> f(4, arr2);
  f /= 2.0;
  CHECK_NEAR(f[0], 2.5, 1e-14, "vector1d operator/=");
}


// Test cvm::vector1d inner product, norm and sum
static void test_vector1d_math()
{
  double arr[] = {3.0, 4.0, 0.0, 0.0};
  cvm::vector1d<double> v(4, arr);

  // Inner product
  CHECK_NEAR(v * v, 25.0, 1e-14, "vector1d inner product with itself");

  double arr2[] = {1.0, 2.0, 2.0, 0.0};
  cvm::vector1d<double> v2(4, arr2);
  CHECK_NEAR(v * v2, 11.0, 1e-14, "vector1d inner product");

  // Norm2 and norm
  CHECK_NEAR(v.norm2(), 25.0, 1e-14, "vector1d norm2");
  CHECK_NEAR(v.norm(), 5.0, 1e-14, "vector1d norm");

  // Sum
  double arr3[] = {1.0, 2.0, 3.0, 4.0};
  cvm::vector1d<double> v3(4, arr3);
  CHECK_NEAR(v3.sum(), 10.0, 1e-14, "vector1d sum");
}


// Test cvm::vector1d slice and sliceassign
static void test_vector1d_slice()
{
  double arr[] = {1.0, 2.0, 3.0, 4.0, 5.0};
  cvm::vector1d<double> v(5, arr);

  // Slice [1, 3) -> {2.0, 3.0}
  cvm::vector1d<double> sl = v.slice(1, 3);
  CHECK(sl.size() == 2, "vector1d slice size");
  CHECK(sl[0] == 2.0 && sl[1] == 3.0, "vector1d slice values");

  // Sliceassign: replace elements [1,3) with {10.0, 20.0}
  double rep[] = {10.0, 20.0};
  cvm::vector1d<double> repl(2, rep);
  v.sliceassign(1, 3, repl);
  CHECK(v[0] == 1.0 && v[1] == 10.0 && v[2] == 20.0 && v[3] == 4.0 && v[4] == 5.0,
        "vector1d sliceassign");
}


// Test cvm::vector1d serialization (to_simple_string / from_simple_string)
static void test_vector1d_serialization()
{
  double arr[] = {1.5, 2.5, 3.5};
  cvm::vector1d<double> v(3, arr);

  std::string s = v.to_simple_string();
  CHECK(s.size() > 0, "vector1d to_simple_string not empty");

  cvm::vector1d<double> v2(3);
  int err = v2.from_simple_string(s);
  CHECK(err == COLVARS_OK, "vector1d from_simple_string return OK");
  CHECK_NEAR(v2[0], 1.5, 1e-10, "vector1d round-trip [0]");
  CHECK_NEAR(v2[1], 2.5, 1e-10, "vector1d round-trip [1]");
  CHECK_NEAR(v2[2], 3.5, 1e-10, "vector1d round-trip [2]");
}


int main(int argc, char *argv[])
{
  test_rvector_construction();
  test_rvector_arithmetic();
  test_rvector_products();
  test_rvector_norm();
  test_vector1d_construction();
  test_vector1d_arithmetic();
  test_vector1d_math();
  test_vector1d_slice();
  test_vector1d_serialization();

  if (test_failures > 0) {
    std::cerr << test_failures << " test(s) failed." << std::endl;
    return 1;
  }

  std::cout << "All colvartypes vector tests passed." << std::endl;
  return 0;
}
