#include <gtest/gtest.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <regex>

// Helper to run the executable and capture its output
std::string RunOptimizationAndGetOutput(double alpha, int n) {
  std::string command = "echo \"" + std::to_string(alpha) + " " + std::to_string(n) + "\" | ./integral_optimization";
  
  char buffer[128];
  std::string result = "";
  FILE* pipe = popen(command.c_str(), "r");
  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer, sizeof buffer, pipe) != nullptr) {
    result += buffer;
  }
  pclose(pipe);
  return result;
}

// Helper to extract the infimum result using regex
double ExtractInfimum(const std::string& output) {
  std::regex inf_regex("Infinum \\(result\\)\\s*=\\s*([0-9\\.]+)");
  std::smatch match;
  if (std::regex_search(output, match, inf_regex)) {
    return std::stod(match[1]);
  }
  return -1.0;
}

TEST(BaselineTest, Alpha5) {
  // tau = (PI / 2) / n. For tau ~ 1e-2, n ~ 1.57 / 0.01 = 157
  std::string output = RunOptimizationAndGetOutput(5.0, 157);
  double infimum = ExtractInfimum(output);
  
  // Expected: 1.051751810209
  EXPECT_NEAR(infimum, 1.05175, 1e-4);
}

TEST(BaselineTest, Alpha20) {
  std::string output = RunOptimizationAndGetOutput(20.0, 157);
  double infimum = ExtractInfimum(output);
  
  // Expected: 0.723388841
  EXPECT_NEAR(infimum, 0.72338, 1e-4);
}
