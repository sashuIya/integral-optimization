#include "system_state.h"

#include <gtest/gtest.h>

using integral_optimization::SystemState;

TEST(SystemStateTest, DefaultConstructor) {
  SystemState state;
  EXPECT_DOUBLE_EQ(state.x1, 0.0);
  EXPECT_DOUBLE_EQ(state.x2, 0.0);
  EXPECT_DOUBLE_EQ(state.p1, 0.0);
  EXPECT_DOUBLE_EQ(state.p2, 0.0);
}

TEST(SystemStateTest, Addition) {
  SystemState a{1.0, 2.0, -3.0, 4.0};
  SystemState b{5.0, -6.0, 7.0, 8.0};
  SystemState c = a + b;
  EXPECT_DOUBLE_EQ(c.x1, 6.0);
  EXPECT_DOUBLE_EQ(c.x2, -4.0);
  EXPECT_DOUBLE_EQ(c.p1, 4.0);
  EXPECT_DOUBLE_EQ(c.p2, 12.0);
}

TEST(SystemStateTest, ScalarMultiplication) {
  SystemState a{1.0, 2.0, 3.0, 4.0};
  SystemState c = a * 2.0;
  SystemState d = 3.0 * a;

  EXPECT_DOUBLE_EQ(c.x1, 2.0);
  EXPECT_DOUBLE_EQ(c.x2, 4.0);
  EXPECT_DOUBLE_EQ(c.p1, 6.0);
  EXPECT_DOUBLE_EQ(c.p2, 8.0);

  EXPECT_DOUBLE_EQ(d.x1, 3.0);
  EXPECT_DOUBLE_EQ(d.x2, 6.0);
  EXPECT_DOUBLE_EQ(d.p1, 9.0);
  EXPECT_DOUBLE_EQ(d.p2, 12.0);
}

TEST(SystemStateTest, ScalarDivision) {
  SystemState a{2.0, 4.0, 6.0, 8.0};
  SystemState c = a / 2.0;

  EXPECT_DOUBLE_EQ(c.x1, 1.0);
  EXPECT_DOUBLE_EQ(c.x2, 2.0);
  EXPECT_DOUBLE_EQ(c.p1, 3.0);
  EXPECT_DOUBLE_EQ(c.p2, 4.0);
}
