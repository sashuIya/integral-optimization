// Copyright 2026 Alexander Lapin
// Main entry point for the Integral Optimization Solver.
// This application solves a two-point boundary value problem using the
// Shooting Method (Runge-Kutta 4 + Newton's method with Fedorenko
// normalization) to find the infimum of a given functional.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "math_constants.h"
#include "ode_solver.h"
#include "shooting_method.h"
#include "system_state.h"

using integral_optimization::ApplyRungeKutta;
using integral_optimization::IntegrateFunctional;
using integral_optimization::kPi;
using integral_optimization::SolveBoundaryValueProblem;
using integral_optimization::SystemState;

// Defines the system of differential equations derived from the
// Pontryagin Maximum Principle for this specific control problem.
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

// L2 Norm calculation for error vectors.
double CalculateNorm(double x, double y) { return std::sqrt(x * x + y * y); }

struct Config {
  double alpha = 5.0;
  int steps = 157;
  double p1 = 0.0;
  double p2 = 0.0;
  bool p_provided = false;
  bool use_grid_search = true;
};

Config ParseArgs(int argc, char** argv) {
  Config config;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--alpha" && i + 1 < argc) {
      config.alpha = std::stod(argv[++i]);
    } else if (arg == "--steps" && i + 1 < argc) {
      config.steps = std::stoi(argv[++i]);
    } else if (arg == "--p1" && i + 1 < argc) {
      config.p1 = std::stod(argv[++i]);
      config.p_provided = true;
    } else if (arg == "--p2" && i + 1 < argc) {
      config.p2 = std::stod(argv[++i]);
      config.p_provided = true;
    } else if (arg == "--no-grid-search") {
      config.use_grid_search = false;
    } else if (arg == "--help" || arg == "-h") {
      std::cout << "Usage: integral_optimization [options]\n"
                << "Options:\n"
                << "  --alpha <value>       Parameter alpha (default: 5.0)\n"
                << "  --steps <value>       Number of Runge-Kutta steps (default: 157)\n"
                << "  --p1 <value>          Initial guess for p1(0)\n"
                << "  --p2 <value>          Initial guess for p2(0)\n"
                << "  --no-grid-search      Disable grid search for initial guess\n";
      std::exit(0);
    }
  }
  return config;
}

int main(int argc, char** argv) {
  Config config = ParseArgs(argc, argv);

  double tau = (kPi / 2.0) / config.steps;
  double error = 0.0;

  double pp1 = 0.0, pp2 = 0.0;

  if (config.p_provided) {
    pp1 = config.p1;
    pp2 = config.p2;
  } else {
    // Fallback guesses based on alpha
    int alpha_num = static_cast<int>(std::ceil(config.alpha));
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
  }

  double best_p1 = pp1;
  double best_p2 = pp2;
  double min_error = 1e+300;

  SystemState boundary_condition;
  boundary_condition.x2 = 0.0;

  // Grid search for better initial guesses
  if (config.use_grid_search) {
    for (double p1 = pp1 - 1; p1 < pp1 + 1; p1 += 0.01) {
      for (double p2 = pp2 - 1; p2 < pp2 + 1; p2 += 0.01) {
        boundary_condition.p1 = p1;
        boundary_condition.p2 = p2;
        boundary_condition.x1 = boundary_condition.p1;

        SystemState result_vector = ApplyRungeKutta(boundary_condition, tau, config.steps, false,
                                                    config.alpha, GetFunctionVector);

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

  SystemState result_vector = SolveBoundaryValueProblem(boundary_condition, tau, config.steps,
                                                        config.alpha, GetFunctionVector, error);

  // Apply Runge Kutta one last time to write the final path to u.gnuplot
  ApplyRungeKutta(result_vector, tau, config.steps, true, config.alpha, GetFunctionVector);

  double result =
      IntegrateFunctional(result_vector, 0.5 * kPi, config.steps, config.alpha, GetFunctionVector);

  printf("\n------------------------------RESULTS-------------------------------\n");
  printf("Runge-kutta step (tau)                       =    %.16lf\n", tau);
  printf("Integration step (h_sim)                     =    %.16lf\n", tau);
  printf("First unknown boundary condition (p1 (0))    =    %.16lf\n", result_vector.p1);
  printf("Second unknown boundary condition (p2 (0))   =    %.16lf\n", result_vector.p2);
  printf("Infinum (result)                             =    %.16lf\n\n\n", result);
  printf("Error                                        =    %.16lf\n\n\n", error);

  return 0;
}
