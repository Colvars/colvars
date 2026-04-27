// -*- c++ -*-

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvarvalue.h"


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


// -------------------------------------------------------
// Scalar tests
// -------------------------------------------------------

static void test_scalar_construction()
{
  // Default: type_scalar, value 0
  colvarvalue cv;
  CHECK(cv.type() == colvarvalue::type_scalar, "default colvarvalue is scalar");
  CHECK_NEAR(cv.real_value, 0.0, 1e-14, "default scalar value is 0");

  // Construct from real
  colvarvalue cv2(3.14);
  CHECK(cv2.type() == colvarvalue::type_scalar, "colvarvalue from real is scalar");
  CHECK_NEAR(cv2.real_value, 3.14, 1e-14, "scalar value from real");

  // Construct from type flag
  colvarvalue cv3(colvarvalue::type_scalar);
  CHECK(cv3.type() == colvarvalue::type_scalar, "colvarvalue type_scalar construction");
  CHECK_NEAR(cv3.real_value, 0.0, 1e-14, "scalar value is 0 after type construction");

  // Copy constructor
  colvarvalue cv4(cv2);
  CHECK(cv4.type() == colvarvalue::type_scalar, "copy constructor preserves type");
  CHECK_NEAR(cv4.real_value, 3.14, 1e-14, "copy constructor preserves value");
}


static void test_scalar_arithmetic()
{
  colvarvalue a(2.0);
  colvarvalue b(3.0);

  // Addition
  colvarvalue sum = a + b;
  CHECK(sum.type() == colvarvalue::type_scalar, "scalar sum type");
  CHECK_NEAR(sum.real_value, 5.0, 1e-14, "scalar sum value");

  // Subtraction
  colvarvalue diff = b - a;
  CHECK_NEAR(diff.real_value, 1.0, 1e-14, "scalar subtraction");

  // Multiplication by real
  colvarvalue scaled = 2.0 * a;
  CHECK_NEAR(scaled.real_value, 4.0, 1e-14, "scalar * real");

  colvarvalue scaled2 = a * 3.0;
  CHECK_NEAR(scaled2.real_value, 6.0, 1e-14, "scalar * real (reversed)");

  // Division by real
  colvarvalue divided = b / 3.0;
  CHECK_NEAR(divided.real_value, 1.0, 1e-14, "scalar / real");

  // In-place operators
  colvarvalue c(2.0);
  c += b;
  CHECK_NEAR(c.real_value, 5.0, 1e-14, "scalar +=");

  colvarvalue d(5.0);
  d -= a;
  CHECK_NEAR(d.real_value, 3.0, 1e-14, "scalar -=");

  colvarvalue e(2.0);
  e *= 4.0;
  CHECK_NEAR(e.real_value, 8.0, 1e-14, "scalar *=");

  colvarvalue f(6.0);
  f /= 3.0;
  CHECK_NEAR(f.real_value, 2.0, 1e-14, "scalar /=");

  // Inner product
  cvm::real ip = a * b;
  CHECK_NEAR(ip, 6.0, 1e-14, "scalar inner product");
}


static void test_scalar_norm_sum()
{
  colvarvalue a(3.0);
  CHECK_NEAR(a.norm2(), 9.0, 1e-14, "scalar norm2");
  CHECK_NEAR(a.norm(), 3.0, 1e-14, "scalar norm");
  CHECK_NEAR(a.sum(), 3.0, 1e-14, "scalar sum");

  colvarvalue neg(-4.0);
  CHECK_NEAR(neg.norm2(), 16.0, 1e-14, "negative scalar norm2");
}


static void test_scalar_dist2()
{
  colvarvalue a(2.0);
  colvarvalue b(5.0);

  CHECK_NEAR(a.dist2(b), 9.0, 1e-14, "scalar dist2");
  CHECK_NEAR(a.dist2(a), 0.0, 1e-14, "scalar dist2 with itself");

  colvarvalue grad = a.dist2_grad(b);
  CHECK(grad.type() == colvarvalue::type_scalar, "scalar dist2_grad type");
  CHECK_NEAR(grad.real_value, -6.0, 1e-14, "scalar dist2_grad: 2*(a-b)");
}


static void test_scalar_interpolate()
{
  colvarvalue a(0.0);
  colvarvalue b(10.0);

  colvarvalue mid = colvarvalue::interpolate(a, b, 0.5);
  CHECK(mid.type() == colvarvalue::type_scalar, "scalar interpolate type");
  CHECK_NEAR(mid.real_value, 5.0, 1e-14, "scalar interpolate midpoint");

  colvarvalue quarter = colvarvalue::interpolate(a, b, 0.25);
  CHECK_NEAR(quarter.real_value, 2.5, 1e-14, "scalar interpolate quarter point");
}


static void test_scalar_ones_reset()
{
  colvarvalue a(3.0);

  // set_ones() sets all components to 1
  colvarvalue ones(colvarvalue::type_scalar);
  ones.set_ones();
  CHECK(ones.type() == colvarvalue::type_scalar, "scalar set_ones() type");
  CHECK_NEAR(ones.real_value, 1.0, 1e-14, "scalar set_ones() value");

  a.reset();
  CHECK_NEAR(a.real_value, 0.0, 1e-14, "scalar reset()");
}


// -------------------------------------------------------
// 3-vector tests
// -------------------------------------------------------

static void test_3vector_construction()
{
  // Construct from rvector with type_3vector (default)
  cvm::rvector v(1.0, 2.0, 3.0);
  colvarvalue cv(v);
  CHECK(cv.type() == colvarvalue::type_3vector, "3vector type");
  CHECK(cv.rvector_value.x == 1.0 && cv.rvector_value.y == 2.0 &&
        cv.rvector_value.z == 3.0, "3vector values");

  // From type flag
  colvarvalue cv2(colvarvalue::type_3vector);
  CHECK(cv2.type() == colvarvalue::type_3vector, "3vector from type flag");
  CHECK(cv2.rvector_value.x == 0.0, "3vector zero-initialized");

  // Copy constructor
  colvarvalue cv3(cv);
  CHECK(cv3.type() == colvarvalue::type_3vector, "3vector copy constructor type");
  CHECK(cv3.rvector_value.x == 1.0, "3vector copy constructor value");
}


static void test_3vector_arithmetic()
{
  cvm::rvector va(1.0, 2.0, 3.0);
  cvm::rvector vb(4.0, 5.0, 6.0);
  colvarvalue a(va);
  colvarvalue b(vb);

  // Addition
  colvarvalue sum = a + b;
  CHECK(sum.type() == colvarvalue::type_3vector, "3vector sum type");
  CHECK_NEAR(sum.rvector_value.x, 5.0, 1e-14, "3vector sum x");
  CHECK_NEAR(sum.rvector_value.y, 7.0, 1e-14, "3vector sum y");
  CHECK_NEAR(sum.rvector_value.z, 9.0, 1e-14, "3vector sum z");

  // Subtraction
  colvarvalue diff = b - a;
  CHECK_NEAR(diff.rvector_value.x, 3.0, 1e-14, "3vector subtraction x");
  CHECK_NEAR(diff.rvector_value.y, 3.0, 1e-14, "3vector subtraction y");
  CHECK_NEAR(diff.rvector_value.z, 3.0, 1e-14, "3vector subtraction z");

  // Scalar multiplication
  colvarvalue scaled = 2.0 * a;
  CHECK_NEAR(scaled.rvector_value.x, 2.0, 1e-14, "3vector scaled x");

  // Scalar division
  colvarvalue divided = a / 2.0;
  CHECK_NEAR(divided.rvector_value.x, 0.5, 1e-14, "3vector divided x");

  // Inner product
  cvm::real ip = a * b;
  CHECK_NEAR(ip, 32.0, 1e-14, "3vector inner product");
}


static void test_3vector_norm_sum()
{
  cvm::rvector v(3.0, 4.0, 0.0);
  colvarvalue cv(v);
  CHECK_NEAR(cv.norm2(), 25.0, 1e-14, "3vector norm2");
  CHECK_NEAR(cv.norm(), 5.0, 1e-14, "3vector norm");
  CHECK_NEAR(cv.sum(), 7.0, 1e-14, "3vector sum");
}


static void test_3vector_dist2()
{
  cvm::rvector va(1.0, 0.0, 0.0);
  cvm::rvector vb(0.0, 1.0, 0.0);
  colvarvalue a(va);
  colvarvalue b(vb);

  CHECK_NEAR(a.dist2(a), 0.0, 1e-14, "3vector dist2 with itself");
  CHECK_NEAR(a.dist2(b), 2.0, 1e-14, "3vector dist2");

  // dist2_grad: 2*(a - b)
  colvarvalue grad = a.dist2_grad(b);
  CHECK(grad.type() == colvarvalue::type_3vector, "3vector dist2_grad type");
  CHECK_NEAR(grad.rvector_value.x, 2.0, 1e-14, "3vector dist2_grad x");
  CHECK_NEAR(grad.rvector_value.y, -2.0, 1e-14, "3vector dist2_grad y");
  CHECK_NEAR(grad.rvector_value.z, 0.0, 1e-14, "3vector dist2_grad z");
}


static void test_3vector_interpolate()
{
  cvm::rvector va(0.0, 0.0, 0.0);
  cvm::rvector vb(2.0, 4.0, 6.0);
  colvarvalue a(va);
  colvarvalue b(vb);

  colvarvalue mid = colvarvalue::interpolate(a, b, 0.5);
  CHECK(mid.type() == colvarvalue::type_3vector, "3vector interpolate type");
  CHECK_NEAR(mid.rvector_value.x, 1.0, 1e-14, "3vector interpolate x");
  CHECK_NEAR(mid.rvector_value.y, 2.0, 1e-14, "3vector interpolate y");
  CHECK_NEAR(mid.rvector_value.z, 3.0, 1e-14, "3vector interpolate z");
}


static void test_3vector_ones_reset()
{
  cvm::rvector v(3.0, 4.0, 5.0);
  colvarvalue a(v);

  colvarvalue ones(colvarvalue::type_3vector);
  ones.set_ones();
  CHECK(ones.type() == colvarvalue::type_3vector, "3vector set_ones() type");
  CHECK_NEAR(ones.rvector_value.x, 1.0, 1e-14, "3vector set_ones() x");
  CHECK_NEAR(ones.rvector_value.y, 1.0, 1e-14, "3vector set_ones() y");
  CHECK_NEAR(ones.rvector_value.z, 1.0, 1e-14, "3vector set_ones() z");

  a.reset();
  CHECK_NEAR(a.rvector_value.x, 0.0, 1e-14, "3vector reset() x");
  CHECK_NEAR(a.rvector_value.y, 0.0, 1e-14, "3vector reset() y");
  CHECK_NEAR(a.rvector_value.z, 0.0, 1e-14, "3vector reset() z");
}


// -------------------------------------------------------
// Quaternion colvarvalue tests
// -------------------------------------------------------

static void test_quaternion_cv_construction()
{
  cvm::quaternion q(0.5, 0.5, 0.5, 0.5);
  colvarvalue cv(q);
  CHECK(cv.type() == colvarvalue::type_quaternion, "quaternion cv type");
  CHECK_NEAR(cv.quaternion_value.q0, 0.5, 1e-14, "quaternion cv q0");

  // From type flag
  colvarvalue cv2(colvarvalue::type_quaternion);
  CHECK(cv2.type() == colvarvalue::type_quaternion, "quaternion cv from type flag");
}


static void test_quaternion_cv_arithmetic()
{
  cvm::quaternion qa(1.0, 0.0, 0.0, 0.0);
  cvm::quaternion qb(0.0, 1.0, 0.0, 0.0);
  colvarvalue a(qa);
  colvarvalue b(qb);

  colvarvalue sum = a + b;
  CHECK(sum.type() == colvarvalue::type_quaternion, "quaternion cv sum type");
  CHECK_NEAR(sum.quaternion_value.q0, 1.0, 1e-14, "quaternion cv sum q0");
  CHECK_NEAR(sum.quaternion_value.q1, 1.0, 1e-14, "quaternion cv sum q1");

  colvarvalue scaled = 2.0 * a;
  CHECK_NEAR(scaled.quaternion_value.q0, 2.0, 1e-14, "quaternion cv scaled q0");
}


static void test_quaternion_cv_norm()
{
  cvm::quaternion q(0.5, 0.5, 0.5, 0.5);
  colvarvalue cv(q);
  CHECK_NEAR(cv.norm2(), 1.0, 1e-14, "quaternion cv norm2");
  CHECK_NEAR(cv.norm(), 1.0, 1e-14, "quaternion cv norm");
  CHECK_NEAR(cv.sum(), 2.0, 1e-14, "quaternion cv sum");
}


static void test_quaternion_cv_dist2()
{
  cvm::quaternion q1(1.0, 0.0, 0.0, 0.0);
  cvm::real const angle = PI / 2.0;
  cvm::quaternion q2(std::cos(angle / 2.0), 0.0, 0.0, std::sin(angle / 2.0));
  colvarvalue a(q1);
  colvarvalue b(q2);

  CHECK_NEAR(a.dist2(a), 0.0, 1e-14, "quaternion cv dist2 with itself");
  // Inner product = cos(pi/4), geodesic distance = pi/4, dist2 = (pi/4)^2
  CHECK_NEAR(a.dist2(b), (PI / 4.0) * (PI / 4.0), 1e-10,
             "quaternion cv dist2 for 90-degree rotation");
}


static void test_quaternion_cv_apply_constraints()
{
  // apply_constraints() on a quaternion should normalize it
  cvm::quaternion q(2.0, 0.0, 0.0, 0.0);
  colvarvalue cv(q);
  cv.apply_constraints();
  CHECK_NEAR(cv.quaternion_value.norm(), 1.0, 1e-14,
             "quaternion cv apply_constraints normalizes");
}


// -------------------------------------------------------
// Generic vector (type_vector) tests
// -------------------------------------------------------

static void test_vector_cv_construction()
{
  double vdata[] = {1.0, 2.0, 3.0, 4.0, 5.0};
  cvm::vector1d<cvm::real> v(5, vdata);
  colvarvalue cv(v, colvarvalue::type_vector);
  CHECK(cv.type() == colvarvalue::type_vector, "vector cv type");
  CHECK(cv.vector1d_value.size() == 5, "vector cv size");
  CHECK_NEAR(cv.vector1d_value[0], 1.0, 1e-14, "vector cv [0]");
  CHECK_NEAR(cv.vector1d_value[4], 5.0, 1e-14, "vector cv [4]");
}


static void test_vector_cv_arithmetic()
{
  double va[] = {1.0, 2.0, 3.0};
  double vb[] = {4.0, 5.0, 6.0};
  cvm::vector1d<cvm::real> a(3, va);
  cvm::vector1d<cvm::real> b(3, vb);
  colvarvalue cva(a, colvarvalue::type_vector);
  colvarvalue cvb(b, colvarvalue::type_vector);

  colvarvalue sum = cva + cvb;
  CHECK(sum.type() == colvarvalue::type_vector, "vector cv sum type");
  CHECK_NEAR(sum.vector1d_value[0], 5.0, 1e-14, "vector cv sum [0]");
  CHECK_NEAR(sum.vector1d_value[1], 7.0, 1e-14, "vector cv sum [1]");
  CHECK_NEAR(sum.vector1d_value[2], 9.0, 1e-14, "vector cv sum [2]");

  colvarvalue diff = cvb - cva;
  CHECK_NEAR(diff.vector1d_value[0], 3.0, 1e-14, "vector cv diff [0]");

  colvarvalue scaled = 2.0 * cva;
  CHECK_NEAR(scaled.vector1d_value[0], 2.0, 1e-14, "vector cv scaled [0]");

  // Inner product
  cvm::real ip = cva * cvb;
  CHECK_NEAR(ip, 32.0, 1e-14, "vector cv inner product");
}


static void test_vector_cv_norm_sum()
{
  double vdata[] = {3.0, 4.0, 0.0};
  cvm::vector1d<cvm::real> v(3, vdata);
  colvarvalue cv(v, colvarvalue::type_vector);

  CHECK_NEAR(cv.norm2(), 25.0, 1e-14, "vector cv norm2");
  CHECK_NEAR(cv.norm(), 5.0, 1e-14, "vector cv norm");
  CHECK_NEAR(cv.sum(), 7.0, 1e-14, "vector cv sum");
}


static void test_vector_cv_dist2()
{
  double va[] = {1.0, 0.0, 0.0};
  double vb[] = {0.0, 1.0, 0.0};
  cvm::vector1d<cvm::real> a(3, va);
  cvm::vector1d<cvm::real> b(3, vb);
  colvarvalue cva(a, colvarvalue::type_vector);
  colvarvalue cvb(b, colvarvalue::type_vector);

  CHECK_NEAR(cva.dist2(cva), 0.0, 1e-14, "vector cv dist2 with itself");
  CHECK_NEAR(cva.dist2(cvb), 2.0, 1e-14, "vector cv dist2");

  colvarvalue grad = cva.dist2_grad(cvb);
  CHECK(grad.type() == colvarvalue::type_vector, "vector cv dist2_grad type");
  CHECK_NEAR(grad.vector1d_value[0], 2.0, 1e-14, "vector cv dist2_grad [0]");
  CHECK_NEAR(grad.vector1d_value[1], -2.0, 1e-14, "vector cv dist2_grad [1]");
  CHECK_NEAR(grad.vector1d_value[2], 0.0, 1e-14, "vector cv dist2_grad [2]");
}


static void test_vector_cv_interpolate()
{
  double va[] = {0.0, 0.0};
  double vb[] = {4.0, 8.0};
  cvm::vector1d<cvm::real> a(2, va);
  cvm::vector1d<cvm::real> b(2, vb);
  colvarvalue cva(a, colvarvalue::type_vector);
  colvarvalue cvb(b, colvarvalue::type_vector);

  colvarvalue mid = colvarvalue::interpolate(cva, cvb, 0.5);
  CHECK(mid.type() == colvarvalue::type_vector, "vector cv interpolate type");
  CHECK_NEAR(mid.vector1d_value[0], 2.0, 1e-14, "vector cv interpolate [0]");
  CHECK_NEAR(mid.vector1d_value[1], 4.0, 1e-14, "vector cv interpolate [1]");
}


static void test_vector_cv_ones_reset()
{
  double vdata[] = {5.0, 6.0, 7.0};
  cvm::vector1d<cvm::real> v(3, vdata);
  colvarvalue cv(v, colvarvalue::type_vector);

  colvarvalue ones(v, colvarvalue::type_vector);
  ones.set_ones();
  CHECK(ones.type() == colvarvalue::type_vector, "vector cv set_ones() type");
  CHECK_NEAR(ones.vector1d_value[0], 1.0, 1e-14, "vector cv set_ones() [0]");
  CHECK_NEAR(ones.vector1d_value[1], 1.0, 1e-14, "vector cv set_ones() [1]");
  CHECK_NEAR(ones.vector1d_value[2], 1.0, 1e-14, "vector cv set_ones() [2]");

  cv.reset();
  for (size_t i = 0; i < cv.vector1d_value.size(); i++) {
    CHECK_NEAR(cv.vector1d_value[i], 0.0, 1e-14, "vector cv reset()");
  }
}


// -------------------------------------------------------
// is_derivative() tests
// -------------------------------------------------------

static void test_is_derivative()
{
  colvarvalue cv_u3v(colvarvalue::type_unit3vector);
  cv_u3v.is_derivative();
  CHECK(cv_u3v.type() == colvarvalue::type_unit3vectorderiv,
        "unit3vector is_derivative -> unit3vectorderiv");

  colvarvalue cv_q(colvarvalue::type_quaternion);
  cv_q.is_derivative();
  CHECK(cv_q.type() == colvarvalue::type_quaternionderiv,
        "quaternion is_derivative -> quaternionderiv");
}


// -------------------------------------------------------
// type_desc / type_keyword / num_dimensions tests
// -------------------------------------------------------

static void test_type_metadata()
{
  CHECK(colvarvalue::type_desc(colvarvalue::type_scalar) == "scalar number",
        "type_desc scalar");
  CHECK(colvarvalue::type_desc(colvarvalue::type_3vector) == "3-dimensional vector",
        "type_desc 3vector");
  CHECK(colvarvalue::type_desc(colvarvalue::type_unit3vector) == "3-dimensional unit vector",
        "type_desc unit3vector");
  CHECK(colvarvalue::type_desc(colvarvalue::type_quaternion) == "4-dimensional unit quaternion",
        "type_desc quaternion");
  CHECK(colvarvalue::type_desc(colvarvalue::type_vector) == "n-dimensional vector",
        "type_desc vector");

  CHECK(colvarvalue::type_keyword(colvarvalue::type_scalar) == "scalar",
        "type_keyword scalar");
  CHECK(colvarvalue::type_keyword(colvarvalue::type_3vector) == "vector3",
        "type_keyword 3vector");
  CHECK(colvarvalue::type_keyword(colvarvalue::type_quaternion) == "unit_quaternion",
        "type_keyword quaternion");
  CHECK(colvarvalue::type_keyword(colvarvalue::type_vector) == "vector",
        "type_keyword vector");

  CHECK(colvarvalue::num_dimensions(colvarvalue::type_scalar) == 1,
        "num_dimensions scalar");
  CHECK(colvarvalue::num_dimensions(colvarvalue::type_3vector) == 3,
        "num_dimensions 3vector");
  CHECK(colvarvalue::num_dimensions(colvarvalue::type_quaternion) == 4,
        "num_dimensions quaternion");
  // vector has no fixed dimensions
  CHECK(colvarvalue::num_dimensions(colvarvalue::type_vector) == 0,
        "num_dimensions vector (dynamic)");
}


// -------------------------------------------------------
// size() tests
// -------------------------------------------------------

static void test_size()
{
  colvarvalue cv_s(1.0);
  CHECK(cv_s.size() == 1, "scalar size() == 1");

  colvarvalue cv_v(cvm::rvector(1.0, 2.0, 3.0));
  CHECK(cv_v.size() == 3, "3vector size() == 3");

  colvarvalue cv_q(cvm::quaternion(1.0, 0.0, 0.0, 0.0));
  CHECK(cv_q.size() == 4, "quaternion size() == 4");

  double vdata[] = {1.0, 2.0, 3.0, 4.0, 5.0};
  cvm::vector1d<cvm::real> v(5, vdata);
  colvarvalue cv_vec(v, colvarvalue::type_vector);
  CHECK(cv_vec.size() == 5, "generic vector size() == 5");
}


int main(int argc, char *argv[])
{
  test_scalar_construction();
  test_scalar_arithmetic();
  test_scalar_norm_sum();
  test_scalar_dist2();
  test_scalar_interpolate();
  test_scalar_ones_reset();

  test_3vector_construction();
  test_3vector_arithmetic();
  test_3vector_norm_sum();
  test_3vector_dist2();
  test_3vector_interpolate();
  test_3vector_ones_reset();

  test_quaternion_cv_construction();
  test_quaternion_cv_arithmetic();
  test_quaternion_cv_norm();
  test_quaternion_cv_dist2();
  test_quaternion_cv_apply_constraints();

  test_vector_cv_construction();
  test_vector_cv_arithmetic();
  test_vector_cv_norm_sum();
  test_vector_cv_dist2();
  test_vector_cv_interpolate();
  test_vector_cv_ones_reset();

  test_is_derivative();
  test_type_metadata();
  test_size();

  if (test_failures > 0) {
    std::cerr << test_failures << " test(s) failed." << std::endl;
    return 1;
  }

  std::cout << "All colvarvalue tests passed." << std::endl;
  return 0;
}
