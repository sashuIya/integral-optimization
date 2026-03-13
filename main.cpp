#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "math_constants.h"
#include "ode_solver.h"
#include "system_state.h"

using integral_optimization::ApplyRungeKutta;
using integral_optimization::IntegrateFunctional;
using integral_optimization::kDelta;
using integral_optimization::kEps;
using integral_optimization::kMaxSteps;
using integral_optimization::kPi;
using integral_optimization::SystemState;

SystemState get_function_vector(const SystemState& current_vector, double alpha) {
  SystemState result_vector;

  result_vector.x1 = current_vector.x2;
  result_vector.x2 =
      current_vector.p2 - current_vector.x1 * std::exp(-1.0 * alpha * current_vector.x1);
  result_vector.p1 = std::exp(-1.0 * alpha * current_vector.x1) * (1 - alpha * current_vector.x1) *
                     current_vector.p2;
  result_vector.p2 = -1.0 * current_vector.p1;

  return result_vector;
}

double norm(double x, double y) { return std::sqrt(x * x + y * y); }

SystemState newton_solver(SystemState boundary_condition, double tau, int n, double alpha,
                          double& error) {
  double gamma, det, a1, a2, next_a1, next_a2, x11, x12, x21, x22, k1, k2;
  double next_error;
  int counter = 0;
  SystemState result_vector1, result_vector2, starting_condition;

  a1 = boundary_condition.p1;
  a2 = boundary_condition.p2;

  starting_condition.x2 = boundary_condition.x2;

  // Find Jacobi matrix elements
  starting_condition.p1 = a1;
  starting_condition.p2 = a2;
  starting_condition.x1 = a1;

  result_vector1 = ApplyRungeKutta(starting_condition, tau, n, false, alpha, get_function_vector);

  starting_condition.p1 = a1 + kDelta;
  starting_condition.p2 = a2;
  starting_condition.x1 = a1 + kDelta;

  result_vector2 = ApplyRungeKutta(starting_condition, tau, n, false, alpha, get_function_vector);

  x11 = (result_vector2.p1 - result_vector1.p1) / kDelta;
  x21 = (result_vector2.p2 - result_vector1.p2) / kDelta;

  starting_condition.p1 = a1;
  starting_condition.p2 = a2 + kDelta;
  starting_condition.x1 = a1;

  result_vector2 = ApplyRungeKutta(starting_condition, tau, n, false, alpha, get_function_vector);

  x12 = (result_vector2.p1 - result_vector1.p1) / kDelta;
  x22 = (result_vector2.p2 - result_vector1.p2) / kDelta;

  // Calculate residual with Fedorenko normalization
  k1 = x11 * x11 + x12 * x12;
  k2 = x21 * x21 + x22 * x22;

  error = std::sqrt(result_vector1.x1 * result_vector1.x1 / k1 +
                    (result_vector1.x2 - 1.0) * (result_vector1.x2 - 1.0) / k2);

  while (counter < kMaxSteps) {
    gamma = 0.1;

    if (error < kEps) {
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

    while (next_error > error) {
      next_a1 = a1 - gamma * add_a1;
      next_a2 = a2 - gamma * add_a2;

      gamma *= 0.5;

      starting_condition.p1 = next_a1;
      starting_condition.p2 = next_a2;
      starting_condition.x1 = next_a1;

      result_vector1 =
          ApplyRungeKutta(starting_condition, tau, n, false, alpha, get_function_vector);

      starting_condition.p1 = next_a1 + kDelta;
      starting_condition.p2 = next_a2;
      starting_condition.x1 = next_a1 + kDelta;

      result_vector2 =
          ApplyRungeKutta(starting_condition, tau, n, false, alpha, get_function_vector);

      x11 = (result_vector2.p1 - result_vector1.p1) / kDelta;
      x21 = (result_vector2.p2 - result_vector1.p2) / kDelta;

      starting_condition.p1 = next_a1;
      starting_condition.p2 = next_a2 + kDelta;
      starting_condition.x1 = next_a1;

      result_vector2 =
          ApplyRungeKutta(starting_condition, tau, n, false, alpha, get_function_vector);

      x12 = (result_vector2.p1 - result_vector1.p1) / kDelta;
      x22 = (result_vector2.p2 - result_vector1.p2) / kDelta;

      k1 = x11 * x11 + x12 * x12;
      k2 = x21 * x21 + x22 * x22;

      next_error = std::sqrt(result_vector1.x1 * result_vector1.x1 / k1 +
                             (result_vector1.x2 - 1.0) * (result_vector1.x2 - 1.0) / k2);
    }
    a1 = next_a1;
    a2 = next_a2;

    error = next_error;
    ++counter;
  }

  starting_condition.p1 = a1;
  starting_condition.p2 = a2;
  starting_condition.x1 = a1;

  return starting_condition;
}

int main() {
  double alpha, tau, error;
  int n;

  std::cout << "alpha:\n";
  std::cin >> alpha;

  std::cout << "Runge-Kutta steps count:\n";
  std::cin >> n;

  tau = (kPi / 2.0) / n;

  SystemState boundary_condition;
  boundary_condition.x2 = 0.0;
  double min_error = 1e+300;
  double best_p1 = 0.0, best_p2 = 0.0;

  double pp1 = 0.0, pp2 = 0.0;
  int alpha_num = static_cast<int>(std::ceil(alpha));
  switch (alpha_num) {
    case 5:
      pp1 = -0.05508700;
      pp2 = -1.22477297;
      break;
    case 10:
      pp1 = 0.147;
      pp2 = -0.734;
      break;
    case 15:
      pp1 = 0.256;
      pp2 = -0.606;
      break;
    case 20:
      pp1 = 0.41193;
      pp2 = -0.44501;
      break;
    case 25:
      pp1 = 0.661;
      pp2 = -0.636;
      break;
    default:
      pp1 = 0.0;
      pp2 = 0.0;
      break;
  }

  best_p1 = pp1;
  best_p2 = pp2;

  // Grid search
  int perform_grid_search = 1;
  if (perform_grid_search) {
    for (double p1 = pp1 - 1; p1 < pp1 + 1; p1 += 0.01) {
      for (double p2 = pp2 - 1; p2 < pp2 + 1; p2 += 0.01) {
        boundary_condition.p1 = p1;
        boundary_condition.p2 = p2;
        boundary_condition.x1 = boundary_condition.p1;

        SystemState result_vector =
            ApplyRungeKutta(boundary_condition, tau, n, false, alpha, get_function_vector);

        double error_vector[2];
        error_vector[0] = result_vector.x1;
        error_vector[1] = result_vector.x2 - 1.0;

        double residual = norm(error_vector[0], error_vector[1]);
        if (residual < min_error) {
          min_error = residual;
          best_p1 = p1;
          best_p2 = p2;
        }
      }
    }
  }

  boundary_condition.p1 = best_p1;
  boundary_condition.p2 = best_p2;
  boundary_condition.x1 = boundary_condition.p1;

  SystemState result_vector =
      ApplyRungeKutta(boundary_condition, tau, n, false, alpha, get_function_vector);

  result_vector = newton_solver(boundary_condition, tau, n, alpha, error);
  ApplyRungeKutta(result_vector, tau, n, true, alpha, get_function_vector);

  double result = IntegrateFunctional(result_vector, 0.5 * kPi, n, alpha, get_function_vector);
  printf("\n------------------------------RESULTS-------------------------------\n");
  printf("Runge-kutta step (tau)                       =    %.16lf\n", tau);
  printf("Integration step (h_sim)                     =    %.16lf\n", tau);
  printf("First unknown boundary condition (p1 (0))    =    %.16lf\n", result_vector.p1);
  printf("Second unknown boundary condition (p2 (0))   =    %.16lf\n", result_vector.p2);
  printf("Infinum (result)                             =    %.16lf\n\n\n", result);
  printf("Error                                        =    %.16lf\n\n\n", error);

  return 0;
}
