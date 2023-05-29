/*
 * Copyright 2012-2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

// check memory allocation in some method
#define EIGEN_RUNTIME_NO_MALLOC

// includes
// std
#include <iostream>

// boost
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE LogDiff test
#include <boost/test/unit_test.hpp>

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

// Exponential of matrices
#include <unsupported/Eigen/MatrixFunctions>

using namespace sva;
using namespace Eigen;

BOOST_AUTO_TEST_CASE(Cbrt)
{
  auto testd = [](double x, double eps = std::numeric_limits<double>::epsilon())
  {
    double y = details::cbrt(x);
    BOOST_CHECK_LE((y * y * y - x), 2 * std::abs(x) * eps);
  };

  auto testf = [](float x, float eps = std::numeric_limits<float>::epsilon())
  {
    float y = details::cbrt(x);
    BOOST_CHECK_LE((y * y * y - x), 2 * std::abs(x) * eps);
  };

  for(int i = -12; i <= 12; ++i) testd(std::pow(10., i));
  for(int i = -12; i <= 12; ++i) testd(-std::pow(10., i));
  testd(0);
  for(int i = -6; i <= 6; ++i) testf(std::pow(10.f, static_cast<float>(i)));
  for(int i = -6; i <= 6; ++i) testf(-std::pow(10.f, static_cast<float>(i)));
  testf(0);
}

BOOST_AUTO_TEST_CASE(SO3JacF2)
{
  auto testd = [](double x, double res, double eps = std::numeric_limits<double>::epsilon())
  {
    double f2 = details::SO3JacF2(x);
    BOOST_CHECK_SMALL(std::abs(f2 - res) / res, eps);
    f2 = details::SO3JacF2(-x);
    BOOST_CHECK_SMALL(std::abs(f2 - res) / res, eps);
  };

  auto testf = [](float x, float res, float eps = 4 * std::numeric_limits<float>::epsilon())
  {
    float f2 = details::SO3JacF2(x);
    BOOST_CHECK_SMALL(std::abs(f2 - res) / res, eps);
    f2 = details::SO3JacF2(-x);
    BOOST_CHECK_SMALL(std::abs(f2 - res) / res, eps);
  };

  // Ground truth numbers were computed with matlab vpa function and 50 digits precision.
  testd(2, 0.0894768460164173242483950);
  // The two next computations do not use the taylor approximation yet, but are subtracting numbers quite close, hence a
  // small loss of accuracy
  testd(1., 0.0847561391437740403659902, 5e-16);
  testd(1e-1, 0.0833472255299274574976738, 1e-13);
  testd(1e-2, 0.0833334722225529108796317);
  testd(1e-3, 0.0833333347222222552910061);
  testd(1e-4, 0.0833333333472222222255291);
  testd(1e-5, 0.0833333333334722222222225);
  testd(1e-6, 0.0833333333333347222222222);
  testd(1e-7, 0.0833333333333333472222222);
  testd(1e-8, 0.0833333333333333334722222);
  testd(1e-9, 0.0833333333333333333347222);
  testd(1e-10, 0.0833333333333333333333472);

  testf(2.f, 0.0894768460164173242483950f);
  testf(1.f, 0.0847561391437740403659902f);
  testf(1e-1f, 0.0833472255299274574976738f);
  testf(1e-2f, 0.0833334722225529108796317f);
  testf(1e-3f, 0.0833333347222222552910061f);
  testf(1e-4f, 0.0833333333472222222255291f);
  testf(1e-5f, 0.0833333333334722222222225f);
  testf(1e-6f, 0.0833333333333347222222222f);
  testf(1e-7f, 0.0833333333333333472222222f);
  testf(1e-8f, 0.0833333333333333334722222f);
  testf(1e-9f, 0.0833333333333333333347222f);
  testf(1e-10f, 0.0833333333333333333333472f);
}

BOOST_AUTO_TEST_CASE(dSO3JacF2)
{
  std::cout.precision(17);
  auto testd = [](double x, double res, double eps = std::numeric_limits<double>::epsilon())
  {
    double f2 = details::dSO3JacF2(x);
    BOOST_CHECK_SMALL(std::abs(f2 - res) / res, eps);
    f2 = -details::dSO3JacF2(-x);
    BOOST_CHECK_SMALL(std::abs(f2 - res) / res, eps);
  };

  auto testf = [](float x, float res, float eps = std::numeric_limits<float>::epsilon())
  {
    float f2 = details::dSO3JacF2(x);
    BOOST_CHECK_SMALL(std::abs(f2 - res) / res, eps);
    f2 = -details::dSO3JacF2(-x);
    BOOST_CHECK_SMALL(std::abs(f2 - res) / res, eps);
  };

  // Ground truth numbers were computed with matlab vpa function and 100 digits precision (50 is not enough for the
  // smallest entries).

  testd(2., 0.00679694292146532720196937, 1e-14);
  testd(1., 0.0029151856912366650224031, 1e-13);
  // This is about the worst scenario: 0.267 ~ taylor_12_bound, but still using the general formula that subtract two
  // close numbers
  testd(0.267, 0.000744191160073626215859056, 3e-10);
  testd(1e-1, 0.000277910102529934204614325);
  testd(0.058, 0.0001611369228328263168718797); // 0.058 ~ taylor_8_bound
  testd(1e-2, 0.0000277779100534060863262305);
  testd(1e-3, 0.00000277777791005291501322768);
  testd(1e-4, 0.000000277777777910052910102513);
  testd(1e-5, 0.0000000277777777779100529100534, 3e-16);
  testd(1e-6, 0.00000000277777777777791005291005);
  testd(1e-7, 0.000000000277777777777777910052910);
  testd(1e-8, 0.0000000000277777777777777779100529);
  testd(1e-9, 0.00000000000277777777777777777791005);
  testd(1e-10, 0.000000000000277777777777777777777910);

  testf(2.f, 0.00679694292146532720196937f, 3e-6f);
  testf(1.f, 0.0029151856912366650224031f);
  testf(1e-1f, 0.000277910102529934204614325f);
  testf(1e-2f, 0.0000277779100534060863262305f);
  testf(1e-3f, 0.00000277777791005291501322768f);
  testf(1e-4f, 0.000000277777777910052910102513f);
  testf(1e-5f, 0.0000000277777777779100529100534f);
  testf(1e-6f, 0.00000000277777777777791005291005f);
  testf(1e-7f, 0.000000000277777777777777910052910f);
  testf(1e-8f, 0.0000000000277777777777777779100529f);
  testf(1e-9f, 0.00000000000277777777777777777791005f);
  testf(1e-10f, 0.000000000000277777777777777777777910f);
}

BOOST_AUTO_TEST_CASE(SO3RightJacInvTest)
{
  for(int i = 0; i < 1000; ++i)
  {
#if EIGEN_VERSION_AT_LEAST(3, 3, 0)
    Matrix3d R1 = Quaterniond::UnitRandom().toRotationMatrix();
    Matrix3d R2 = Quaterniond::UnitRandom().toRotationMatrix();
#else
    Matrix3d R1 = Matrix3d::Random().householderQr().householderQ();
    Matrix3d R2 = Matrix3d::Random().householderQr().householderQ();
    assert(std::abs(R1.determinant() - 1) < 1e-14); // We use the fact that the QR implementation gives a determinant
    assert(std::abs(R2.determinant() - 1) < 1e-14); // of -(-1)^n with n the size of the matrix. However this is
                                                    // implementation-dependent, so we keep the check.
#endif

    // Computations are getting less precise as we approach the singularity of the log for a rotation angle of pi.
    // We skip those cases so that we can keep a tight bound on the precision of the remaining tests.
    // This loss of precision is not a problem in practice because (i) in the target applications of sva, rotation
    // angles close to pi is no something we should write/encounter and (ii) even if the precision deteriorates, the
    // derivative is still not bad, and provides a correct descent direction for gradient-based computation schemes.
    if((R2 * R1.transpose()).trace() < -0.95)
    {
      --i;
      continue;
    }
    Matrix3d I = Matrix3d::Identity();
    double h = 1e-8;

    Vector3d u0 = rotationError(R1, R2);
    auto invJr = SO3RightJacInv(u0);

    // Finite differences for f(X,Y) = rotationError(X, Y)
    // Jd1 is Df/DX, Jd2 is Df/DY
    Matrix3d J1d, J2d;
    for(int i = 0; i < 3; ++i)
    {
      // R1 "+" h e_i
      Matrix3d R1p = R1 * vector3ToCrossMatrix(Vector3d(h * I.col(i))).exp();
      // R2 "+" h e_i
      Matrix3d R2p = R2 * vector3ToCrossMatrix(Vector3d(h * I.col(i))).exp();
      J1d.col(i) = (rotationError(R1p, R2) - u0) / h; // ( f(R1 "+" h e_i, R2) "-" f(R1,R2) ) / h
      J2d.col(i) = (rotationError(R1, R2p) - u0) / h; // ( f(R1, R2 "+" h e_i) "-" f(R1,R2) ) / h
    }

    BOOST_CHECK_SMALL((J1d - invJr).norm(), 3e-6); // Df/DX = Jr^-1
    BOOST_CHECK_SMALL((J2d + invJr.transpose()).norm(), 3e-6); // Df/DY = -Jl^-1 = -Jr^-T
  }
}

BOOST_AUTO_TEST_CASE(SO3RightJacInvDotTest)
{
  for(int i = 0; i < 10; ++i)
  {
    Vector3d u0 = Vector3d::Random();
    Vector3d du = Vector3d::Random();

    double dt = 1e-8;
    Matrix3d M0 = SO3RightJacInv(u0);
    Matrix3d Mt = SO3RightJacInv(Vector3d(u0 + dt * du));
    Matrix3d dM = (Mt - M0) / dt;
    Matrix3d J0 = SO3RightJacInvDot(u0, du);

    BOOST_CHECK_SMALL((dM - J0).norm(), 1e-6);
  }
}
