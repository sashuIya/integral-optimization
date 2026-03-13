#pragma once

#include <functional>

#include "system_state.h"

namespace integral_optimization {

// Defines a system of ODEs: returns the derivatives of the state given the current state and a
// parameter.
using OdeSystem = std::function<SystemState(const SystemState&, double)>;

// Solves an ODE system using the 4th-order Runge-Kutta method over n steps.
SystemState ApplyRungeKutta(const SystemState& starting_condition, double tau, int n,
                            bool flag_print, double alpha, const OdeSystem& system);

// Computes the objective integrand for a given state.
double Integrand(const SystemState& state);

// Approximates the integral using Simpson's rule over a single sub-interval.
double SimpsonRule(const SystemState& state1, const SystemState& state2, const SystemState& state3,
                   double tau);

// Calculates the integral of the functional over the full time interval.
double IntegrateFunctional(const SystemState& starting_condition, double t, int count_steps,
                           double alpha, const OdeSystem& system);

}  // namespace integral_optimization
