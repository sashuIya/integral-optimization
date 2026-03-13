#include "ode_solver.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace integral_optimization;

// A simple exponential decay/growth ODE for testing Runge-Kutta
// dx1/dt = alpha * x1
SystemState ExponentialOde(const SystemState& state, double alpha) {
  SystemState derivative;
  derivative.x1 = alpha * state.x1;
  return derivative;
}

TEST(OdeSolverTest, RungeKuttaExponentialGrowth) {
  SystemState initial_state;
  initial_state.x1 = 1.0;  // y(0) = 1

  double alpha = 2.0;
  double tau = 0.01;
  int steps = 100;  // t = 1.0

  // Runge-Kutta should approximate y(1) = e^2
  SystemState final_state =
      ApplyRungeKutta(initial_state, tau, steps, false, alpha, ExponentialOde);

  double expected = std::exp(2.0);
  EXPECT_NEAR(final_state.x1, expected, 1e-4);
}

TEST(OdeSolverTest, SimpsonRuleIntegration) {
  // Test Simpson's rule integrating a constant integrand
  SystemState s1;
  s1.p2 = 2.0;  // Integrand = 4.0
  SystemState s2;
  s2.p2 = 2.0;  // Integrand = 4.0
  SystemState s3;
  s3.p2 = 2.0;  // Integrand = 4.0

  double tau = 0.5;
  // Int(4.0) over interval [0, 1] width tau*2 (tau is half step in Simpson) -> wait, standard
  // simpson formula used here: (tau/3)*(f(0) + 4f(1) + f(2)) The integration interval length is
  // 2*tau. Area = 4.0 * 2*tau = 4.0 * 1.0 = 4.0
  double result = SimpsonRule(s1, s2, s3, tau);

  EXPECT_DOUBLE_EQ(result, 4.0);
}
