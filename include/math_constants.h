// Copyright 2026 Alexander Lapin
// Mathematical constants used for numerical optimization.

#ifndef INTEGRAL_OPTIMIZATION_MATH_CONSTANTS_H_
#define INTEGRAL_OPTIMIZATION_MATH_CONSTANTS_H_

namespace integral_optimization {

// Accuracy threshold for Newton's method convergence.
constexpr double kEps = 1e-8;

// Small delta used for numerical differentiation (Jacobi matrix calculation).
constexpr double kDelta = 1e-5;

// Mathematical constant Pi.
constexpr double kPi = 3.1415926535897932;

// Maximum number of iterations allowed for the Newton solver.
constexpr int kMaxSteps = 1000;

}  // namespace integral_optimization

#endif  // INTEGRAL_OPTIMIZATION_MATH_CONSTANTS_H_
