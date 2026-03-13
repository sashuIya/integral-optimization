#pragma once

namespace integral_optimization {

struct SystemState {
  double x1 = 0.0;
  double x2 = 0.0;
  double p1 = 0.0;
  double p2 = 0.0;

  SystemState operator+(const SystemState& other) const;
  SystemState operator*(double scalar) const;
  SystemState operator/(double scalar) const;
};

SystemState operator*(double scalar, const SystemState& state);

}  // namespace integral_optimization
