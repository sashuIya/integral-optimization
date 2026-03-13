// Copyright 2026 Alexander Lapin
// Implementation of the Shooting Method for Boundary Value Problems.

#ifndef INTEGRAL_OPTIMIZATION_SHOOTING_METHOD_H_
#define INTEGRAL_OPTIMIZATION_SHOOTING_METHOD_H_

#include "ode_solver.h"
#include "system_state.h"

namespace integral_optimization {

// Finds the missing initial values for p1(0) and p2(0) to satisfy the
// boundary conditions at the end of the interval (T = PI/2).
//
// uses a modified Newton's method with Fedorenko normalization for convergence.
// final_error: Output parameter containing the final residual norm.
SystemState SolveBoundaryValueProblem(SystemState initial_guess, double tau, int n, double alpha,
                                      const OdeSystem& system, double& final_error);

}  // namespace integral_optimization

#endif  // INTEGRAL_OPTIMIZATION_SHOOTING_METHOD_H_
