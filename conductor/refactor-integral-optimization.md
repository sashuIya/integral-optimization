# Plan: Modernizing Integral Optimization Solver

## Objective
Refactor the existing `main.cpp` into a modular, modern C++ application following the Google C++ Style Guide. The goal is to improve maintainability, performance, and readability while preserving mathematical correctness. Work will be split into small, meaningful commits, with unit tests written for each refactoring step.

## Key Constraints & Standards
1. **File Extensions**: Use `.h` for headers and `.cc` for implementation files.
2. **Language**: Translate all Russian comments, output, and variable names to English.
3. **Style**: Follow the Google C++ Style Guide (enforced via `.clang-format`).
4. **Commits**: Make small, logical, self-contained commits.
5. **Testing**: Write unit tests before or during each refactoring step to ensure correctness (using a framework like Google Test).
6. **Documentation**: Create a comprehensive `README.md` translating the mathematical problem and results from `chm_task2.tex`.

## Phase 1: Infrastructure, Documentation, and Baseline
1. **Setup Environment**: 
   - Add `CMakeLists.txt` with C++17 support, Google Test integration, and strict compilation flags.
   - Add `.clang-format` configured for Google Style.
2. **Documentation**:
   - Create `README.md` based on `chm_task2.tex` (Problem statement, Algorithm, Expected Results).
3. **Baseline Tests**:
   - Write a high-level integration/baseline test that compiles and runs the existing `main.cpp` logic to verify output matches the expected results for $\alpha \in \{0, 5, 10, 15, 20, 25\}$.
4. *Commit: "chore: add build system, formatting, documentation, and baseline tests"*

## Phase 2: Modularization - System State & Constants
1. **Math Constants**: Extract `PI`, `EPS`, `DELTA` into `include/math_constants.h`.
2. **System State**: Extract `FunctionVector` into `include/system_state.h` and `src/system_state.cc`. Rename to `SystemState` (Google Style: `CamelCase` for types). Translate comments.
3. **Unit Tests**: Add tests for `SystemState` operator overloads.
4. **Integration**: Update `main.cpp` to use these new components.
5. *Commit: "refactor: extract SystemState and math constants with unit tests"*

## Phase 3: Modularization - ODE Solver
1. **ODE Solver Module**: Create `include/ode_solver.h` and `src/ode_solver.cc` to encapsulate the Runge-Kutta 4 step (`apply_runge_kutta`) and Simpson's rule (`integration_functions`).
2. **Translation**: Ensure all variables and comments are in English.
3. **Unit Tests**: Write tests for the `OdeSolver` (e.g., integrating a simple known function).
4. **Integration**: Update `main.cpp` to use `OdeSolver`.
5. *Commit: "refactor: extract OdeSolver module with unit tests"*

## Phase 4: Modularization - Shooting Method
1. **Shooting Method Module**: Extract Newton's method (`newton_solver`) logic into `include/shooting_method.h` and `src/shooting_method.cc`.
2. **Unit Tests**: Add tests to ensure the shooting method converges for a known simple problem.
3. **Integration**: Update `main.cpp` to use `ShootingMethod`.
4. *Commit: "refactor: extract ShootingMethod module with unit tests"*

## Phase 5: Modernization & Enhancements
1. **CLI Arguments**: Replace hardcoded `std::cin` inputs in `main.cpp` with command-line argument parsing (e.g., using `getopt` or simple `argc`/`argv` parsing) for $\alpha$, `n`, and initial guesses.
2. **Remove Hardcoded Values**: Eliminate the `switch(alpha)` statement and allow passing initial guesses via CLI or configuration.
3. **Output format**: Output results in a clean format (e.g., CSV) instead of a hardcoded `u.gnuplot` file.
4. *Commit: "feat: add CLI argument parsing and improve output formatting"*

## Verification
- Run the full test suite (baseline + unit tests) to ensure everything passes and matches the original numerical results.
