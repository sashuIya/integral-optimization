// Defines the state of the system for both state and adjoint variables.

#ifndef INTEGRAL_OPTIMIZATION_SYSTEM_STATE_H_
#define INTEGRAL_OPTIMIZATION_SYSTEM_STATE_H_

namespace integral_optimization {

// Represents the four-dimensional state of the system at any given point in time.
// x1, x2 are state variables; p1, p2 are adjoint variables from the
// Pontryagin Maximum Principle.
struct SystemState {
  double x1 = 0.0;
  double x2 = 0.0;
  double p1 = 0.0;
  double p2 = 0.0;

  // Standard arithmetic overloads to support RK4 vector operations.
  SystemState operator+(const SystemState& other) const;
  SystemState operator*(double scalar) const;
  SystemState operator/(double scalar) const;
};

// Allows scalar multiplication on the left side (e.g., 0.5 * state).
SystemState operator*(double scalar, const SystemState& state);

}  // namespace integral_optimization

#endif  // INTEGRAL_OPTIMIZATION_SYSTEM_STATE_H_
