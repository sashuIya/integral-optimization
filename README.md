# Integral Optimization Solver

## Problem Statement

The goal is to find the numerical solution to the following optimal control problem:

$$
\int_0^{\frac{\pi}{2}} u^2 \, dt + x^2(0) \rightarrow \inf
$$

Subject to the constraints:

$$
\dot{x}(0) = x(\pi/2) = 0
$$

$$
\dot{x}(\pi/2) = 1
$$

$$
\ddot{x} + x \exp(-\alpha x) = u
$$

$$
\alpha \in [0.0, 25.0]
$$

## Algorithm

### Reduction to ODE System with Boundary Conditions

By applying the Pontryagin Maximum Principle, the problem is reduced to a system of Ordinary Differential Equations (ODEs) with boundary conditions:

$$
\dot{x}_1 = x_2
$$

$$
\dot{x}_2 = u - x_1 e^{-\alpha x_1}
$$

$$
\dot{p}_1 = e^{-\alpha x_1} (1 - \alpha x_1) p_2
$$

$$
\dot{p}_2 = -p_1
$$

$$
p_1(0) = x_1(0)
$$

$$
p_2 = u
$$

$$
x_2(0) = 0
$$

$$
x_1(\pi/2) = 0
$$

$$
x_2(\pi/2) = 1
$$

### Numerical Solution

Since we have 4 unknown functions but only two boundary conditions at each end, we need to find the missing boundary conditions at the start of the interval ($t=0$). Let:

$$
p_1(0) = a_1
$$

$$
p_2(0) = a_2
$$

We find the parameters $a_1$ and $a_2$ using a modified **Newton's Method** (Shooting Method) with Fedorenko normalization. The underlying ODEs are solved iteratively using the **4th-order Runge-Kutta method**.

The objective functional (the integral) is then calculated using **Simpson's rule**.

## Expected Numerical Results

Using an error tolerance $\epsilon < 10^{-7}$ and step length $\tau = 10^{-2}$, the expected results are:

| $\alpha$ | $\tau$ | $h_{sim}$ | $p_1(0)$ | $p_2(0)$ | result (infimum) |
|----------|--------|-----------|----------|----------|------------------|
| 0        | 1e-2   | 1e-2      | -0.6816222 | -0.4339341 | 0.21701335789 |
| 5        | 1e-2   | 1e-2      | -0.0550866 | -1.224772669 | 1.051751810209 |
| 10       | 1e-2   | 1e-2      | 0.14338431 | -0.727341818 | 0.7684414831 |
| 15       | 1e-2   | 1e-2      | 0.2333231563 | -0.5793537491 | 0.7335585168 |
| 20       | 1e-2   | 1e-2      | 0.2932859841 | -0.5102745366 | 0.723388841 |
| 25       | 1e-2   | 1e-2      | 0.334583262 | -0.461841896 | 0.71827762311 |

## Building and Running

### Prerequisites
- CMake 3.14+
- C++17 compiler (GCC 9+, Clang 10+, or MSVC 2019+)

### Build Steps
```bash
mkdir -p build && cd build
cmake ..
make -j
```

### Running the Program
```bash
./integral_optimization --alpha 20.0 --steps 157
```

### Running Tests
```bash
ctest
```
