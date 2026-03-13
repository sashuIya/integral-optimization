#pragma once

#include "ode_solver.h"
#include "system_state.h"

namespace integral_optimization {

// Solves the two-point boundary value problem using a modified Newton's method
// (Shooting Method) with Fedorenko normalization.
//
// initial_guess: Initial values for the missing boundary conditions (p1(0) and p2(0))
// tau: Runge-Kutta step size
// n: Number of Runge-Kutta steps
// alpha: Parameter alpha in the ODE system
// system: The system of ODEs to integrate
// final_error: Output parameter to store the final error norm
SystemState SolveBoundaryValueProblem(SystemState initial_guess, double tau, int n, double alpha,
                                      const OdeSystem& system, double& final_error);

}  // namespace integral_optimization
