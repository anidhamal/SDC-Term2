#include "PID.h"
#include <iostream>

using namespace std;

/*
* TODO: Complete the PID class.
*/

PID::PID() {
  Kp = 0.0;
  Ki = 0.0;
  Kd = 0.0;

  p_error = 0.0;
  i_error = 0.0;
  d_error = 0.0;
}

PID::~PID() {}

void PID::Init(double kp, double ki, double kd) {
    Kp = kp;
    Ki = ki;
    Kd = kd;
    std::cout << Kp << std::endl;
    std::cout << Ki << std::endl;
    std::cout << Kd << std::endl;
    sum_cte = 0.0;
}

void PID::UpdateError(double cte) {
  p_error = cte;
  d_error = cte - prev_cte;
  prev_cte = cte;
  sum_cte += cte;
  i_error = sum_cte;
}

double PID::TotalError() {
  return Kp * p_error + Kd * d_error + Ki * i_error;
}

