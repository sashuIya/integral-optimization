#include "shooting_method.h"

#include <cmath>

#include "math_constants.h"

namespace integral_optimization {

SystemState SolveBoundaryValueProblem(SystemState initial_guess, double tau, int n, double alpha,
                                      const OdeSystem& system, double& final_error) {
  double gamma, det, a1, a2, next_a1, next_a2, x11, x12, x21, x22, k1, k2;
  double next_error;
  int counter = 0;
  SystemState result_vector1, result_vector2, starting_condition;

  a1 = initial_guess.p1;
  a2 = initial_guess.p2;

  starting_condition.x2 = initial_guess.x2;

  // Find Jacobi matrix elements
  starting_condition.p1 = a1;
  starting_condition.p2 = a2;
  starting_condition.x1 = a1;

  result_vector1 = ApplyRungeKutta(starting_condition, tau, n, false, alpha, system);

  starting_condition.p1 = a1 + kDelta;
  starting_condition.p2 = a2;
  starting_condition.x1 = a1 + kDelta;

  result_vector2 = ApplyRungeKutta(starting_condition, tau, n, false, alpha, system);

  x11 = (result_vector2.p1 - result_vector1.p1) / kDelta;
  x21 = (result_vector2.p2 - result_vector1.p2) / kDelta;

  starting_condition.p1 = a1;
  starting_condition.p2 = a2 + kDelta;
  starting_condition.x1 = a1;

  result_vector2 = ApplyRungeKutta(starting_condition, tau, n, false, alpha, system);

  x12 = (result_vector2.p1 - result_vector1.p1) / kDelta;
  x22 = (result_vector2.p2 - result_vector1.p2) / kDelta;

  // Calculate residual with Fedorenko normalization
  k1 = x11 * x11 + x12 * x12;
  k2 = x21 * x21 + x22 * x22;

  final_error = std::sqrt(result_vector1.x1 * result_vector1.x1 / k1 +
                          (result_vector1.x2 - 1.0) * (result_vector1.x2 - 1.0) / k2);

  while (counter < kMaxSteps) {
    gamma = 0.1;

    if (final_error < kEps) {
      starting_condition.p1 = a1;
      starting_condition.p2 = a2;
      starting_condition.x1 = a1;
      return starting_condition;
    }

    next_error = 1e+300;
    next_a1 = a1;
    next_a2 = a2;

    det = x11 * x22 - x12 * x21;
    double add_a1 = (result_vector1.x1 * x22 - (result_vector1.x2 - 1.0) * x12) / det;
    double add_a2 = ((result_vector1.x2 - 1.0) * x11 - result_vector1.x1 * x21) / det;

    while (next_error > final_error) {
      next_a1 = a1 - gamma * add_a1;
      next_a2 = a2 - gamma * add_a2;

      gamma *= 0.5;

      starting_condition.p1 = next_a1;
      starting_condition.p2 = next_a2;
      starting_condition.x1 = next_a1;

      result_vector1 = ApplyRungeKutta(starting_condition, tau, n, false, alpha, system);

      starting_condition.p1 = next_a1 + kDelta;
      starting_condition.p2 = next_a2;
      starting_condition.x1 = next_a1 + kDelta;

      result_vector2 = ApplyRungeKutta(starting_condition, tau, n, false, alpha, system);

      x11 = (result_vector2.p1 - result_vector1.p1) / kDelta;
      x21 = (result_vector2.p2 - result_vector1.p2) / kDelta;

      starting_condition.p1 = next_a1;
      starting_condition.p2 = next_a2 + kDelta;
      starting_condition.x1 = next_a1;

      result_vector2 = ApplyRungeKutta(starting_condition, tau, n, false, alpha, system);

      x12 = (result_vector2.p1 - result_vector1.p1) / kDelta;
      x22 = (result_vector2.p2 - result_vector1.p2) / kDelta;

      k1 = x11 * x11 + x12 * x12;
      k2 = x21 * x21 + x22 * x22;

      next_error = std::sqrt(result_vector1.x1 * result_vector1.x1 / k1 +
                             (result_vector1.x2 - 1.0) * (result_vector1.x2 - 1.0) / k2);
    }
    a1 = next_a1;
    a2 = next_a2;

    final_error = next_error;
    ++counter;
  }

  starting_condition.p1 = a1;
  starting_condition.p2 = a2;
  starting_condition.x1 = a1;

  return starting_condition;
}

}  // namespace integral_optimization
