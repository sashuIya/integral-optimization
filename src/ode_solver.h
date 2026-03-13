// Methods for solving ODE systems and integrating functionals.

#ifndef INTEGRAL_OPTIMIZATION_ODE_SOLVER_H_
#define INTEGRAL_OPTIMIZATION_ODE_SOLVER_H_

#include <functional>

#include "system_state.h"

namespace integral_optimization {

// Defines a system of ODEs as a function object.
// Returns the derivatives (SystemState) based on current state and parameter.
using OdeSystem = std::function<SystemState(const SystemState&, double)>;

// Solves an ODE system over a specified number of steps using RK4.
// n: Number of integration steps.
// flag_print: If true, writes results (p2 over time) to trajectory.csv.
// alpha: The parameter alpha used by the ODE system.
SystemState ApplyRungeKutta(const SystemState& starting_condition, double tau, int n,
                            bool flag_print, double alpha, const OdeSystem& system);

// Computes the value of the integrand (u^2 = p2^2) for the objective function.
double Integrand(const SystemState& state);

// Approximates the integral over a single sub-interval [t, t + 2*tau].
double SimpsonRule(const SystemState& state1, const SystemState& state2, const SystemState& state3,
                   double tau);

// Integrates the functional across the entire domain using Simpson's Rule.
// t: The final time (usually PI/2).
// count_steps: Number of integration intervals.
double IntegrateFunctional(const SystemState& starting_condition, double t, int count_steps,
                           double alpha, const OdeSystem& system);

}  // namespace integral_optimization

#endif  // INTEGRAL_OPTIMIZATION_ODE_SOLVER_H_
