#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
using namespace std;


const double EPS = 1e-8;
const int MAX_STEPS = 1000;
const double DELTA = 1e-5;
const double PI    = 3.1415926535897932;

struct FunctionVector
{
  double x1, x2;
  double p1, p2;
};

FunctionVector operator + (FunctionVector first_vector, FunctionVector second_vector)
{
  FunctionVector result_vector;

  result_vector.x1 = first_vector.x1 + second_vector.x1;
  result_vector.x2 = first_vector.x2 + second_vector.x2;
  result_vector.p1 = first_vector.p1 + second_vector.p1;
  result_vector.p2 = first_vector.p2 + second_vector.p2;

  return result_vector;
}

FunctionVector operator * (double number, FunctionVector vector)
{
  FunctionVector result_vector;

  result_vector.x1 = number * vector.x1;
  result_vector.x2 = number * vector.x2;
  result_vector.p1 = number * vector.p1;
  result_vector.p2 = number * vector.p2;

  return result_vector;
}

FunctionVector operator * (FunctionVector vector, double number)
{
  FunctionVector result_vector;

  result_vector.x1 = number * vector.x1;
  result_vector.x2 = number * vector.x2;
  result_vector.p1 = number * vector.p1;
  result_vector.p2 = number * vector.p2;

  return result_vector;
}

FunctionVector operator / (FunctionVector vector, double number)
{
  FunctionVector result_vector;

  result_vector.x1 = vector.x1 / number;
  result_vector.x2 = vector.x2 / number;
  result_vector.p1 = vector.p1 / number;
  result_vector.p2 = vector.p2 / number;

  return result_vector;
}

FunctionVector get_function_vector (
    FunctionVector current_vector, 
    double alpha
  )
{
  FunctionVector result_vector;
  result_vector.x1 = current_vector.x2;
  result_vector.x2 = 1e+13;
  result_vector.p1 = 1e+13;
  result_vector.p2 = -current_vector.p1;

  result_vector.x1 = current_vector.x2;
  result_vector.x2 = current_vector.p2 - current_vector.x1 * exp((-1.0) * alpha * current_vector.x1);
  result_vector.p1 = exp((-1.0) * alpha * current_vector.x1) * (1 - alpha * current_vector.x1) * current_vector.p2;
  result_vector.p2 = (-1.0) * current_vector.p1;

  return result_vector;
}

// Алгоритм Рунге-Кутта
FunctionVector apply_runge_kutta (
    FunctionVector starting_condition, 
    double tau, 
    int n, 
    bool flag_print, 
    double alpha
  )
{
  int index;
  FunctionVector k1, k2, k3, k4, answer;
  FILE *out = 0;
  double t = 0.0;


  if (flag_print)
    out = fopen ("./u.gnuplot", "w");

  answer = starting_condition;

  for (index = 0; index < n; ++index)
    {
      if (flag_print)
        fprintf (out, "%lf %lf\n", t, answer.p2);
      k1 = get_function_vector (answer, alpha) * tau;
      k2 = get_function_vector (answer + 0.5 * k1, alpha) * tau;
      k3 = get_function_vector (answer + 0.5 * k2, alpha) * tau;
      k4 = get_function_vector (answer + k3, alpha) * tau;
      answer = answer + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
      t += tau;
    }
  if (flag_print)
    {
       fprintf (out, "%lf %lf\n", t, answer.p2);
       fclose  (out);
     }  

  return answer;
}

double integrand (FunctionVector functionVector)
{
  return functionVector.p2 * functionVector.p2;
}

double simpson_rule (FunctionVector FunctionVector1, FunctionVector FunctionVector2, FunctionVector FunctionVector3, double tau)
{
  return (tau / 3) * (integrand(FunctionVector1) + 4 * integrand(FunctionVector2) + integrand(FunctionVector3));
}

double integration_functions (FunctionVector starting_condition, double t, int count_steps, double alpha)
{
  FunctionVector result_vector1, result_vector2, start;
  double tau = t / (2 * count_steps);
  double sum = 0.0;

  start = starting_condition;

  for (int index = 0; index < count_steps; ++index)
    {
      result_vector1 = apply_runge_kutta(start, tau, 1, false, alpha);
      result_vector2 = apply_runge_kutta(result_vector1, tau, 1, false, alpha);
      sum += simpson_rule(start, result_vector1, result_vector2, tau);
      start = result_vector2;
    }

//  sum += starting_condition.x1 * starting_condition.x1;
  return sum;
}

double norm (double x, double y)
{
  return sqrt(x * x + y * y);
}

FunctionVector newton_solver (FunctionVector boundary_condition, double tau, int n, double alpha, double &error)
{
  double gamma, det, a1, a2, next_a1, next_a2, x11, x12, x21, x22, k1, k2;
  double next_error;
  int counter = 0;
  FunctionVector result_vector1, result_vector2, starting_condition;

  a1 = boundary_condition.p1;
  a2 = boundary_condition.p2;

  starting_condition.x2 = boundary_condition.x2;

  // Найдём элементы матрицы Якоби
  starting_condition.p1 = a1;
  starting_condition.p2 = a2;
  starting_condition.x1 = a1;

  result_vector1 = apply_runge_kutta (starting_condition, tau, n, false, alpha);

  starting_condition.p1 = a1 + DELTA;
  starting_condition.p2 = a2;
  starting_condition.x1 = a1 + DELTA;

  result_vector2 = apply_runge_kutta (starting_condition, tau, n, false, alpha);

  x11 = (result_vector2.p1 - result_vector1.p1) / DELTA;
  x21 = (result_vector2.p2 - result_vector1.p2) / DELTA;

  starting_condition.p1 = a1;
  starting_condition.p2 = a2 + DELTA;
  starting_condition.x1 = a1;

  result_vector2 = apply_runge_kutta (starting_condition, tau, n, false, alpha);

  x12 = (result_vector2.p1 - result_vector1.p1) / DELTA;
  x22 = (result_vector2.p2 - result_vector1.p2) / DELTA;

  cout << x11 << " " << x12 << " " << endl << x21 << " " << x22 << endl;
  // Считаем невязку с нормировкой Федоренко
  k1 = x11 * x11 + x12 * x12;
  k2 = x21 * x21 + x22 * x22;

  error = sqrt (result_vector1.x1 * result_vector1.x1 / k1 + (result_vector1.x2 - 1.0) * (result_vector1.x2 - 1.0) / k2);
    
  while (counter < MAX_STEPS)
    {
      // Инициализация итерации
      gamma = 0.1;

      printf ("Error = %.14lf\n", error);
      // Если норма достаточно маленькая, заканчиваем алгоритм
      if(error < EPS)
        {
          starting_condition.p1 = a1;
          starting_condition.p2 = a2;
          starting_condition.x1 = a1;
          return starting_condition;
        }

      next_error = 1e+300;
      next_a1 = a1;
      next_a2 = a2;
      // Если нет, ищем начальные условия так, чтобы норма уменьшилась относительно предыдущего шага
      
      det = x11 * x22 - x12 * x21;
//      cout << "det: " << det << endl;
      double add_a1 = (result_vector1.x1 * x22 - (result_vector1.x2 - 1.0) * x12) / det;
      double add_a2 = ((result_vector1.x2 - 1.0) * x11 - result_vector1.x1 * x21) / det;
      while (next_error > error)
        {
          next_a1 = a1 - gamma * add_a1;
          next_a2 = a2 - gamma * add_a2;

          gamma *=  0.5;

          starting_condition.p1 = next_a1;
          starting_condition.p2 = next_a2;
          starting_condition.x1 = next_a1;
 //         cout << starting_condition.p2 << endl;

          result_vector1 = apply_runge_kutta (starting_condition, tau, n, false, alpha);

  //        cout << result_vector1.x1 << endl;

          starting_condition.p1 = next_a1 + DELTA;
          starting_condition.p2 = next_a2;
          starting_condition.x1 = next_a1 + DELTA;

          result_vector2 = apply_runge_kutta (starting_condition, tau, n, false, alpha);

          x11 = (result_vector2.p1 - result_vector1.p1) / DELTA;
          x21 = (result_vector2.p2 - result_vector1.p2) / DELTA;

          starting_condition.p1 = next_a1;
          starting_condition.p2 = next_a2 + DELTA;
          starting_condition.x1 = next_a1;

          result_vector2 = apply_runge_kutta (starting_condition, tau, n, false, alpha);

          x12 = (result_vector2.p1 - result_vector1.p1) / DELTA;
          x22 = (result_vector2.p2 - result_vector1.p2) / DELTA;

          // Считаем невязку с нормировкой Федоренко
          k1 = x11 * x11 + x12 * x12;
          k2 = x21 * x21 + x22 * x22;

          next_error = sqrt (result_vector1.x1 * result_vector1.x1 / k1 + (result_vector1.x2 - 1.0) * (result_vector1.x2 - 1.0) / k2);
//          cout << "next_error: " << next_error << endl;
        }
      // Заканчиваем итерацию
      a1 = next_a1;
      a2 = next_a2;

      error = next_error;
      ++counter;
    }

  starting_condition.p1 = a1;
  starting_condition.p2 = a2;
  starting_condition.x1 = a1;

  return starting_condition;
}

int main()
{
  // Ввод данных
  double alpha, tau, n, error;
  // Параметр
  std::cout << "alpha:";
  std::cin >> alpha;

  std::cout << "Runge-Kutta steps count:";
  std::cin >> n;

  // step
  tau = (PI / 2.0) / n;

  // Численное решение
  // Нахождение недостающих начальных условий

  // Алгоритмом Рунге-Кутты решаем задачу
  FunctionVector boundary_condition;
  

  boundary_condition.x2 = 0.0;
  double min_error = 1e+300;
  double best_p1 = 0, best_p2 = 0;

  double pp1, pp2;
  int alpha_num = (int)ceil(alpha);
  switch (alpha_num)
    {
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
        pp1 = 0;
        pp2 = 0;
        break;
    }

  best_p1 = pp1;
  best_p2 = pp2;

  int PER = 1;
  if (PER)
    for (double p1 = pp1 - 1; p1 < pp1 + 1; p1 += 0.01)
      for (double p2 = pp2 - 1; p2 < pp2 + 1; p2 += 0.01)
        {
          boundary_condition.p1 = p1;//-0.398;
          boundary_condition.p2 = p2; //-1.089;
          boundary_condition.x1 = boundary_condition.p1;

          FunctionVector result_vector = apply_runge_kutta (boundary_condition, tau, n, false, alpha);

          double error_vector[2];
          error_vector[0] = result_vector.x1;
          error_vector[1] = result_vector.x2 - 1.0;

          double residual = norm (error_vector[0], error_vector[1]);
          if (residual < min_error)
            {
              min_error = residual;
              best_p1 = p1;
              best_p2 = p2;

              cout << "MINIMUM ERROR: " << min_error << endl;
              cout << "best_p1:       " << best_p1 << endl;
              cout << "best_p2:       " << best_p2 << endl;
            }
        }

  cout << "-------------------------------------" << endl;
  cout << "MINIMUM ERROR: " << min_error << endl;
  cout << "best_p1:       " << best_p1 << endl;
  cout << "best_p2:       " << best_p2 << endl;
  cout << "-------------------------------------" << endl;

  boundary_condition.p1 = best_p1;//-0.398;
  boundary_condition.p2 = best_p2; //-1.089;
  boundary_condition.x1 = boundary_condition.p1;

  FunctionVector result_vector = apply_runge_kutta (boundary_condition, tau, n, false, alpha);

  // Считаем вектор ошибки и невязку
  double error_vector[2];
  error_vector[0] = result_vector.x1;
  error_vector[1] = result_vector.x2 - 1.0;

  double residual = norm (error_vector[0], error_vector[1]);

  // Проверяем, маленькая невязка или нет
  bool newton_flag;
  if (residual < EPS)
    newton_flag = false;
  else
    newton_flag = true;

  // Если нет применяем алгоритм Ньютона
//  if(newton_flag)
  result_vector = newton_solver (boundary_condition, tau, n, alpha, error);
  cout << "p1: " << result_vector.p1 << "; p2: " << result_vector.p2 << endl;

  apply_runge_kutta (result_vector, tau, n, true, alpha);

  // Методом Симпсона считаем значение функционала
  double result = integration_functions (result_vector, 0.5 * PI, n, alpha);
  printf ("\n------------------------------RESULTS-------------------------------\n");
  printf ("Runge-kutta step (tau)                       =    %.16lf\n" , tau);
  printf ("Integration step (h_sim)                     =    %.16lf\n" , tau);
  printf ("First unknown boundary condition (p1 (0))    =    %.16lf\n", result_vector.p1);
  printf ("Second unknown boundary condition (p2 (0))   =    %.16lf\n", result_vector.p2);
  printf ("Infinum (result)                             =    %.16lf\n\n\n", result);
  printf ("Error                                        =    %.16lf\n\n\n", error);

  // Нахождение значения функционала
  // Выдача ответа
  return 0;
}

