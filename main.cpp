#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "math_constants.h"
#include "ode_solver.h"
#include "shooting_method.h"
#include "system_state.h"

using integral_optimization::ApplyRungeKutta;
using integral_optimization::IntegrateFunctional;
using integral_optimization::kPi;
using integral_optimization::SolveBoundaryValueProblem;
using integral_optimization::SystemState;

SystemState GetFunctionVector(const SystemState& current_vector, double alpha) {
  SystemState result_vector;

  result_vector.x1 = current_vector.x2;
  result_vector.x2 =
      current_vector.p2 - current_vector.x1 * std::exp(-1.0 * alpha * current_vector.x1);
  result_vector.p1 = std::exp(-1.0 * alpha * current_vector.x1) * (1 - alpha * current_vector.x1) *
                     current_vector.p2;
  result_vector.p2 = -1.0 * current_vector.p1;

  return result_vector;
}

double CalculateNorm(double x, double y) { return std::sqrt(x * x + y * y); }

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

  // Grid search for better initial guesses
  int perform_grid_search = 1;
  if (perform_grid_search) {
    for (double p1 = pp1 - 1; p1 < pp1 + 1; p1 += 0.01) {
      for (double p2 = pp2 - 1; p2 < pp2 + 1; p2 += 0.01) {
        boundary_condition.p1 = p1;
        boundary_condition.p2 = p2;
        boundary_condition.x1 = boundary_condition.p1;

        SystemState result_vector =
            ApplyRungeKutta(boundary_condition, tau, n, false, alpha, GetFunctionVector);

        double error_vector[2];
        error_vector[0] = result_vector.x1;
        error_vector[1] = result_vector.x2 - 1.0;

        double residual = CalculateNorm(error_vector[0], error_vector[1]);
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
      SolveBoundaryValueProblem(boundary_condition, tau, n, alpha, GetFunctionVector, error);

  // Apply Runge Kutta one last time to write the final path to u.gnuplot
  ApplyRungeKutta(result_vector, tau, n, true, alpha, GetFunctionVector);

  double result = IntegrateFunctional(result_vector, 0.5 * kPi, n, alpha, GetFunctionVector);

  printf("\n------------------------------RESULTS-------------------------------\n");
  printf("Runge-kutta step (tau)                       =    %.16lf\n", tau);
  printf("Integration step (h_sim)                     =    %.16lf\n", tau);
  printf("First unknown boundary condition (p1 (0))    =    %.16lf\n", result_vector.p1);
  printf("Second unknown boundary condition (p2 (0))   =    %.16lf\n", result_vector.p2);
  printf("Infinum (result)                             =    %.16lf\n\n\n", result);
  printf("Error                                        =    %.16lf\n\n\n", error);

  return 0;
}
