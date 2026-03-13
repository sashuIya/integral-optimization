#include "system_state.h"

namespace integral_optimization {

SystemState SystemState::operator+(const SystemState& other) const {
  return {x1 + other.x1, x2 + other.x2, p1 + other.p1, p2 + other.p2};
}

SystemState SystemState::operator*(double scalar) const {
  return {x1 * scalar, x2 * scalar, p1 * scalar, p2 * scalar};
}

SystemState SystemState::operator/(double scalar) const {
  return {x1 / scalar, x2 / scalar, p1 / scalar, p2 / scalar};
}

SystemState operator*(double scalar, const SystemState& state) { return state * scalar; }

}  // namespace integral_optimization
