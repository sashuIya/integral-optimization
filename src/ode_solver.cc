#include "ode_solver.h"

#include <cstdio>

namespace integral_optimization {

SystemState ApplyRungeKutta(const SystemState& starting_condition, double tau, int n,
                            bool flag_print, double alpha, const OdeSystem& system) {
  SystemState k1, k2, k3, k4, answer;
  FILE* out = nullptr;
  double t = 0.0;

  if (flag_print) {
    out = std::fopen("./trajectory.csv", "w");
    if (out) {
      std::fprintf(out, "t,p2\n");
    }
  }

  answer = starting_condition;

  for (int index = 0; index < n; ++index) {
    if (flag_print && out) {
      std::fprintf(out, "%lf,%lf\n", t, answer.p2);
    }
    k1 = system(answer, alpha) * tau;
    k2 = system(answer + k1 * 0.5, alpha) * tau;
    k3 = system(answer + k2 * 0.5, alpha) * tau;
    k4 = system(answer + k3, alpha) * tau;
    answer = answer + (k1 + k2 * 2.0 + k3 * 2.0 + k4) / 6.0;
    t += tau;
  }

  if (flag_print && out) {
    std::fprintf(out, "%lf,%lf\n", t, answer.p2);
    std::fclose(out);
  }

  return answer;
}

double Integrand(const SystemState& state) { return state.p2 * state.p2; }

double SimpsonRule(const SystemState& state1, const SystemState& state2, const SystemState& state3,
                   double tau) {
  return (tau / 3.0) * (Integrand(state1) + 4.0 * Integrand(state2) + Integrand(state3));
}

double IntegrateFunctional(const SystemState& starting_condition, double t, int count_steps,
                           double alpha, const OdeSystem& system) {
  SystemState result_vector1, result_vector2, start;
  double tau = t / (2.0 * count_steps);
  double sum = 0.0;

  start = starting_condition;

  for (int index = 0; index < count_steps; ++index) {
    result_vector1 = ApplyRungeKutta(start, tau, 1, false, alpha, system);
    result_vector2 = ApplyRungeKutta(result_vector1, tau, 1, false, alpha, system);
    sum += SimpsonRule(start, result_vector1, result_vector2, tau);
    start = result_vector2;
  }

  return sum;
}

}  // namespace integral_optimization
