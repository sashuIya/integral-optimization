#include "shooting_method.h"

#include <gtest/gtest.h>

using namespace integral_optimization;

// A simple mock ODE system:
// dx1/dt = p1
// dx2/dt = p2
// dp1/dt = 0
// dp2/dt = 0
// This means x1(t) = x1(0) + p1(0)*t, x2(t) = x2(0) + p2(0)*t
// Our boundary target is x1(T) = 0, x2(T) = 1.
// If x1(0)=0, x2(0)=0, and T=1, then the required p1(0)=0, p2(0)=1.
SystemState MockOdeSystem(const SystemState& state, double /* alpha */) {
  SystemState derivative;
  derivative.x1 = state.p1;
  derivative.x2 = state.p2;
  derivative.p1 = 0.0;
  derivative.p2 = 0.0;
  return derivative;
}

TEST(ShootingMethodTest, SimpleConvergence) {
  SystemState initial_guess;
  initial_guess.x1 = 0.0;
  initial_guess.x2 = 0.0;
  initial_guess.p1 = -0.5;  // bad guess
  initial_guess.p2 = 0.5;   // bad guess

  double tau = 0.1;
  int n = 10;          // T = tau * n = 1.0
  double alpha = 0.0;  // unused
  double final_error = 0.0;

  SystemState result =
      SolveBoundaryValueProblem(initial_guess, tau, n, alpha, MockOdeSystem, final_error);

  EXPECT_NEAR(result.p1, 0.0, 1e-4);
  EXPECT_NEAR(result.p2, 1.0, 1e-4);
  EXPECT_NEAR(final_error, 0.0, 1e-4);
}
