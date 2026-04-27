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


// Test quaternion construction
static void test_quaternion_construction()
{
  // Default constructor: null quaternion
  cvm::quaternion q0;
  CHECK(q0.q0 == 0.0 && q0.q1 == 0.0 && q0.q2 == 0.0 && q0.q3 == 0.0,
        "quaternion default constructor should zero all components");

  // Parameterized constructor
  cvm::quaternion q1(1.0, 0.0, 0.0, 0.0);
  CHECK(q1.q0 == 1.0 && q1.q1 == 0.0 && q1.q2 == 0.0 && q1.q3 == 0.0,
        "quaternion parameterized constructor");

  // From C-array
  cvm::real arr[] = {0.5, 0.5, 0.5, 0.5};
  cvm::quaternion q2(arr);
  CHECK_NEAR(q2.q0, 0.5, 1e-14, "quaternion from array q0");
  CHECK_NEAR(q2.q1, 0.5, 1e-14, "quaternion from array q1");
  CHECK_NEAR(q2.q2, 0.5, 1e-14, "quaternion from array q2");
  CHECK_NEAR(q2.q3, 0.5, 1e-14, "quaternion from array q3");

  // From vector1d
  double vdata[] = {0.5, 0.5, 0.5, 0.5};
  cvm::vector1d<cvm::real> v(4, vdata);
  cvm::quaternion q3(v);
  CHECK_NEAR(q3.q0, 0.5, 1e-14, "quaternion from vector1d");

  // Index access
  cvm::quaternion q4(1.0, 2.0, 3.0, 4.0);
  CHECK(q4[0] == 1.0 && q4[1] == 2.0 && q4[2] == 3.0 && q4[3] == 4.0,
        "quaternion index access");

  // reset()
  cvm::quaternion q5(1.0, 2.0, 3.0, 4.0);
  q5.reset();
  CHECK(q5.q0 == 0.0 && q5.q1 == 0.0 && q5.q2 == 0.0 && q5.q3 == 0.0,
        "quaternion reset()");
}


// Test quaternion arithmetic
static void test_quaternion_arithmetic()
{
  cvm::quaternion p(1.0, 2.0, 3.0, 4.0);
  cvm::quaternion q(5.0, 6.0, 7.0, 8.0);

  // Addition
  cvm::quaternion sum = p + q;
  CHECK(sum.q0 == 6.0 && sum.q1 == 8.0 && sum.q2 == 10.0 && sum.q3 == 12.0,
        "quaternion addition");

  // Subtraction
  cvm::quaternion diff = q - p;
  CHECK(diff.q0 == 4.0 && diff.q1 == 4.0 && diff.q2 == 4.0 && diff.q3 == 4.0,
        "quaternion subtraction");

  // Scalar multiplication
  cvm::quaternion scaled = 2.0 * p;
  CHECK(scaled.q0 == 2.0 && scaled.q1 == 4.0 && scaled.q2 == 6.0 && scaled.q3 == 8.0,
        "quaternion scalar multiplication (scalar*q)");

  cvm::quaternion scaled2 = p * 3.0;
  CHECK(scaled2.q0 == 3.0 && scaled2.q1 == 6.0 && scaled2.q2 == 9.0 && scaled2.q3 == 12.0,
        "quaternion scalar multiplication (q*scalar)");

  // Scalar division
  cvm::quaternion divided = p / 2.0;
  CHECK_NEAR(divided.q0, 0.5, 1e-14, "quaternion division q0");
  CHECK_NEAR(divided.q1, 1.0, 1e-14, "quaternion division q1");
  CHECK_NEAR(divided.q2, 1.5, 1e-14, "quaternion division q2");
  CHECK_NEAR(divided.q3, 2.0, 1e-14, "quaternion division q3");

  // In-place operators
  cvm::quaternion a = p;
  a += q;
  CHECK(a.q0 == 6.0 && a.q1 == 8.0, "quaternion operator+=");

  cvm::quaternion b = q;
  b -= p;
  CHECK(b.q0 == 4.0 && b.q1 == 4.0, "quaternion operator-=");

  cvm::quaternion c = p;
  c *= 2.0;
  CHECK(c.q0 == 2.0 && c.q1 == 4.0, "quaternion operator*=");

  cvm::quaternion d = p;
  d /= 2.0;
  CHECK_NEAR(d.q0, 0.5, 1e-14, "quaternion operator/=");
}


// Test quaternion product (Hamilton product)
static void test_quaternion_product()
{
  // Identity: (1,0,0,0) * q = q
  cvm::quaternion ident(1.0, 0.0, 0.0, 0.0);
  cvm::quaternion q(0.5, 0.5, 0.5, 0.5);
  cvm::quaternion iq = ident * q;
  CHECK_NEAR(iq.q0, q.q0, 1e-14, "quaternion product: identity*q q0");
  CHECK_NEAR(iq.q1, q.q1, 1e-14, "quaternion product: identity*q q1");
  CHECK_NEAR(iq.q2, q.q2, 1e-14, "quaternion product: identity*q q2");
  CHECK_NEAR(iq.q3, q.q3, 1e-14, "quaternion product: identity*q q3");

  // q * conjugate(q) = (|q|^2, 0, 0, 0)
  cvm::quaternion conj_q = q.conjugate();
  cvm::quaternion qcq = q * conj_q;
  CHECK_NEAR(qcq.q0, q.norm2(), 1e-14, "q * conj(q) = |q|^2 (q0)");
  CHECK_NEAR(qcq.q1, 0.0, 1e-14, "q * conj(q) = 0 (q1)");
  CHECK_NEAR(qcq.q2, 0.0, 1e-14, "q * conj(q) = 0 (q2)");
  CHECK_NEAR(qcq.q3, 0.0, 1e-14, "q * conj(q) = 0 (q3)");

  // Non-commutativity: i*j = k, j*i = -k
  cvm::quaternion qi(0.0, 1.0, 0.0, 0.0);  // unit i
  cvm::quaternion qj(0.0, 0.0, 1.0, 0.0);  // unit j
  cvm::quaternion qk(0.0, 0.0, 0.0, 1.0);  // unit k
  cvm::quaternion ij = qi * qj;
  CHECK_NEAR(ij.q0, 0.0, 1e-14, "i*j = k: q0");
  CHECK_NEAR(ij.q1, 0.0, 1e-14, "i*j = k: q1");
  CHECK_NEAR(ij.q2, 0.0, 1e-14, "i*j = k: q2");
  CHECK_NEAR(ij.q3, 1.0, 1e-14, "i*j = k: q3");

  cvm::quaternion ji = qj * qi;
  CHECK_NEAR(ji.q3, -1.0, 1e-14, "j*i = -k: q3");
}


// Test quaternion norm, conjugate, inner product
static void test_quaternion_norm()
{
  cvm::quaternion q(0.5, 0.5, 0.5, 0.5);

  // Norm2 of (0.5, 0.5, 0.5, 0.5) = 1.0
  CHECK_NEAR(q.norm2(), 1.0, 1e-14, "quaternion norm2");
  CHECK_NEAR(q.norm(), 1.0, 1e-14, "quaternion norm");

  cvm::quaternion q2(1.0, 2.0, 3.0, 4.0);
  CHECK_NEAR(q2.norm2(), 30.0, 1e-14, "quaternion norm2 general");
  CHECK_NEAR(q2.norm(), std::sqrt(30.0), 1e-14, "quaternion norm general");

  // Conjugate
  cvm::quaternion conj = q.conjugate();
  CHECK_NEAR(conj.q0, 0.5, 1e-14, "quaternion conjugate q0");
  CHECK_NEAR(conj.q1, -0.5, 1e-14, "quaternion conjugate q1");
  CHECK_NEAR(conj.q2, -0.5, 1e-14, "quaternion conjugate q2");
  CHECK_NEAR(conj.q3, -0.5, 1e-14, "quaternion conjugate q3");

  // Inner product
  cvm::quaternion qa(1.0, 0.0, 0.0, 0.0);
  cvm::quaternion qb(0.0, 1.0, 0.0, 0.0);
  CHECK_NEAR(qa.inner(qb), 0.0, 1e-14, "quaternion inner product (orthogonal)");
  CHECK_NEAR(qa.inner(qa), 1.0, 1e-14, "quaternion inner product (self)");

  cvm::quaternion qc(0.5, 0.5, 0.5, 0.5);
  cvm::quaternion qd(0.5, 0.5, 0.5, 0.5);
  CHECK_NEAR(qc.inner(qd), 1.0, 1e-14, "quaternion inner product (equal unit quaternions)");
}


// Test quaternion set_positive()
static void test_quaternion_set_positive()
{
  cvm::quaternion q1(0.5, 0.5, 0.5, 0.5);
  q1.set_positive();
  CHECK(q1.q0 > 0.0, "set_positive: already positive q0 unchanged");

  cvm::quaternion q2(-0.5, 0.5, 0.5, 0.5);
  q2.set_positive();
  CHECK(q2.q0 > 0.0, "set_positive: flipped negative q0");
  CHECK_NEAR(q2.q0, 0.5, 1e-14, "set_positive: correct q0 after flip");
  CHECK_NEAR(q2.q1, -0.5, 1e-14, "set_positive: correct q1 after flip");
}


// Test quaternion rotation (90° around z-axis)
static void test_quaternion_rotation()
{
  // Rotation of 90 degrees around z-axis:
  // q = (cos(pi/4), 0, 0, sin(pi/4))
  cvm::real const angle = PI / 2.0;
  cvm::quaternion q_rot(std::cos(angle / 2.0), 0.0, 0.0, std::sin(angle / 2.0));

  // Rotate x-axis -> should give y-axis
  cvm::rvector x_axis(1.0, 0.0, 0.0);
  cvm::rvector rotated = q_rot.rotate(x_axis);
  CHECK_NEAR(rotated.x, 0.0, 1e-14, "rotate x by 90 deg around z: x-component");
  CHECK_NEAR(rotated.y, 1.0, 1e-14, "rotate x by 90 deg around z: y-component");
  CHECK_NEAR(rotated.z, 0.0, 1e-14, "rotate x by 90 deg around z: z-component");

  // Rotate y-axis -> should give -x-axis
  cvm::rvector y_axis(0.0, 1.0, 0.0);
  cvm::rvector rotated_y = q_rot.rotate(y_axis);
  CHECK_NEAR(rotated_y.x, -1.0, 1e-14, "rotate y by 90 deg around z: x-component");
  CHECK_NEAR(rotated_y.y, 0.0, 1e-14, "rotate y by 90 deg around z: y-component");
  CHECK_NEAR(rotated_y.z, 0.0, 1e-14, "rotate y by 90 deg around z: z-component");

  // z-axis should be unchanged
  cvm::rvector z_axis(0.0, 0.0, 1.0);
  cvm::rvector rotated_z = q_rot.rotate(z_axis);
  CHECK_NEAR(rotated_z.x, 0.0, 1e-14, "rotate z by 90 deg around z: x-component");
  CHECK_NEAR(rotated_z.y, 0.0, 1e-14, "rotate z by 90 deg around z: y-component");
  CHECK_NEAR(rotated_z.z, 1.0, 1e-14, "rotate z by 90 deg around z: z-component");
}


// Test quaternion rotation_matrix()
static void test_quaternion_rotation_matrix()
{
  // Identity rotation: q = (1, 0, 0, 0) -> should give identity matrix
  cvm::quaternion q_id(1.0, 0.0, 0.0, 0.0);
  cvm::rmatrix R_id = q_id.rotation_matrix();
  CHECK_NEAR(R_id.xx, 1.0, 1e-14, "identity rotation matrix: xx");
  CHECK_NEAR(R_id.yy, 1.0, 1e-14, "identity rotation matrix: yy");
  CHECK_NEAR(R_id.zz, 1.0, 1e-14, "identity rotation matrix: zz");
  CHECK_NEAR(R_id.xy, 0.0, 1e-14, "identity rotation matrix: xy");
  CHECK_NEAR(R_id.xz, 0.0, 1e-14, "identity rotation matrix: xz");
  CHECK_NEAR(R_id.yx, 0.0, 1e-14, "identity rotation matrix: yx");
  CHECK_NEAR(R_id.yz, 0.0, 1e-14, "identity rotation matrix: yz");
  CHECK_NEAR(R_id.zx, 0.0, 1e-14, "identity rotation matrix: zx");
  CHECK_NEAR(R_id.zy, 0.0, 1e-14, "identity rotation matrix: zy");

  // 90 degrees around z-axis: q = (cos(pi/4), 0, 0, sin(pi/4))
  cvm::real const angle = PI / 2.0;
  cvm::quaternion q_rot(std::cos(angle / 2.0), 0.0, 0.0, std::sin(angle / 2.0));
  cvm::rmatrix R = q_rot.rotation_matrix();

  // Rotation matrix for 90 deg around z:
  // [ 0 -1  0 ]
  // [ 1  0  0 ]
  // [ 0  0  1 ]
  CHECK_NEAR(R.xx, 0.0, 1e-14, "90deg-z rotation matrix: xx");
  CHECK_NEAR(R.xy, -1.0, 1e-14, "90deg-z rotation matrix: xy");
  CHECK_NEAR(R.xz, 0.0, 1e-14, "90deg-z rotation matrix: xz");
  CHECK_NEAR(R.yx, 1.0, 1e-14, "90deg-z rotation matrix: yx");
  CHECK_NEAR(R.yy, 0.0, 1e-14, "90deg-z rotation matrix: yy");
  CHECK_NEAR(R.yz, 0.0, 1e-14, "90deg-z rotation matrix: yz");
  CHECK_NEAR(R.zx, 0.0, 1e-14, "90deg-z rotation matrix: zx");
  CHECK_NEAR(R.zy, 0.0, 1e-14, "90deg-z rotation matrix: zy");
  CHECK_NEAR(R.zz, 1.0, 1e-14, "90deg-z rotation matrix: zz");

  // Apply rotation matrix to x-axis and check the result
  cvm::rvector x_axis(1.0, 0.0, 0.0);
  cvm::rvector result = R * x_axis;
  CHECK_NEAR(result.x, 0.0, 1e-14, "rotation matrix * x-axis: x");
  CHECK_NEAR(result.y, 1.0, 1e-14, "rotation matrix * x-axis: y");
  CHECK_NEAR(result.z, 0.0, 1e-14, "rotation matrix * x-axis: z");
}


// Test quaternion dist2 and dist2_grad
static void test_quaternion_dist2()
{
  // dist2 of identical quaternions should be 0
  cvm::quaternion q(1.0, 0.0, 0.0, 0.0);
  CHECK_NEAR(q.dist2(q), 0.0, 1e-14, "quaternion dist2(q,q)");

  // dist2 of q and -q should also be 0 (they represent the same rotation)
  cvm::quaternion qneg(-1.0, 0.0, 0.0, 0.0);
  CHECK_NEAR(q.dist2(qneg), 0.0, 1e-14, "quaternion dist2(q, -q)");

  // 90-degree rotation around z: q1=(1,0,0,0), q2=(cos(pi/4),0,0,sin(pi/4))
  // Inner product = cos(pi/4), omega = acos(cos(pi/4)) = pi/4
  // dist2 = (pi/4)^2 (geodesic distance on S^3 is half the rotation angle)
  cvm::real const angle = PI / 2.0;
  cvm::quaternion q2(std::cos(angle / 2.0), 0.0, 0.0, std::sin(angle / 2.0));
  cvm::real d2 = q.dist2(q2);
  CHECK_NEAR(d2, (PI / 4.0) * (PI / 4.0), 1e-12,
             "quaternion dist2 for 90-degree rotation");

  // dist2_grad at identical quaternions: should be zero
  cvm::quaternion grad = q.dist2_grad(q);
  CHECK_NEAR(grad.norm2(), 0.0, 1e-14, "quaternion dist2_grad at identical quaternions");
}


// Test quaternion cosine
static void test_quaternion_cosine()
{
  // cosine of identical quaternions should be 1.0
  cvm::quaternion q(1.0, 0.0, 0.0, 0.0);
  CHECK_NEAR(q.cosine(q), 1.0, 1e-14, "quaternion cosine(q, q)");

  // cos(omega) = 2*|inner|^2 - 1
  // For orthogonal 4D unit vectors, inner=0, cos=-1
  cvm::quaternion q2(0.0, 1.0, 0.0, 0.0);
  CHECK_NEAR(q.cosine(q2), -1.0, 1e-14, "quaternion cosine(q, orthogonal q)");
}


// Test quaternion as_vector() and get_vector()
static void test_quaternion_conversions()
{
  cvm::quaternion q(1.0, 2.0, 3.0, 4.0);

  // as_vector()
  cvm::vector1d<cvm::real> v = q.as_vector();
  CHECK(v.size() == 4, "quaternion as_vector size");
  CHECK_NEAR(v[0], 1.0, 1e-14, "quaternion as_vector [0]");
  CHECK_NEAR(v[1], 2.0, 1e-14, "quaternion as_vector [1]");
  CHECK_NEAR(v[2], 3.0, 1e-14, "quaternion as_vector [2]");
  CHECK_NEAR(v[3], 4.0, 1e-14, "quaternion as_vector [3]");

  // get_vector()
  cvm::rvector vec = q.get_vector();
  CHECK_NEAR(vec.x, 2.0, 1e-14, "quaternion get_vector x");
  CHECK_NEAR(vec.y, 3.0, 1e-14, "quaternion get_vector y");
  CHECK_NEAR(vec.z, 4.0, 1e-14, "quaternion get_vector z");
}


int main(int argc, char *argv[])
{
  test_quaternion_construction();
  test_quaternion_arithmetic();
  test_quaternion_product();
  test_quaternion_norm();
  test_quaternion_set_positive();
  test_quaternion_rotation();
  test_quaternion_rotation_matrix();
  test_quaternion_dist2();
  test_quaternion_cosine();
  test_quaternion_conversions();

  if (test_failures > 0) {
    std::cerr << test_failures << " test(s) failed." << std::endl;
    return 1;
  }

  std::cout << "All quaternion tests passed." << std::endl;
  return 0;
}
