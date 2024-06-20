#include "xm_motor.h"
#include <Arduino.h>
#include <math.h>
#include <BasicLinearAlgebra.h>

#define PI 3.1415926535897931

#define l1 0.09264
#define l2 0.150
#define l3 0.150
#define l4 0.050

#define Sgn11 -1
#define Sgn12 -1
#define Sgn13 1
#define Bias13 0.366519

#define Sgn21 1
#define Sgn22 1
#define Sgn23 -1
#define Bias23 0.366519

#define Sgn31 1
#define Sgn32 1
#define Sgn33 -1
#define Bias33 0.366519

#define Sgn41 -1
#define Sgn42 -1
#define Sgn43 1
#define Bias43 0.366519

#define SPD_LIMIT 15
#define CUR_LIMIT 5
#define KP 0.2
#define KI 0.0158
#define DelayMS 150

constexpr double M_SQRT1_3 = 1 / sqrt(3);
constexpr float DEG2RAD = PI / 180;
constexpr float RAD2DEG = 180 / PI;
constexpr double af = PI / 4;
constexpr double bt = 3 * PI / 4;

const double saf = sin(af);
const double caf = cos(af);
const double sbt = sin(bt);
const double cbt = cos(bt);

int ErrorFlag = 0;

using namespace BLA;

BLA::Matrix<3, 1> w = {M_SQRT1_3, M_SQRT1_3, M_SQRT1_3};
BLA::Matrix<3, 3> ww;

BLA::Matrix<6, 1> S1 = {1, 0, 0, 0, 0, 0};
BLA::Matrix<6, 1> S2 = {0, 1, 0, 0, 0, 0};
BLA::Matrix<6, 1> S3 = {sbt, cbt, 0, l3 *cbt, -l3 *cbt, -l1 *sbt};

const int CAN_INT_PIN = 2;
void setup()
{
  while (!Serial)
    ;
  Serial.begin(115200);

  // 初始化can设置
  xm_can_start();

  // 位置置0

  //   motor_pos_zero(11);
  //   delay(50);
  //   motor_pos_zero(12);
  //   delay(50);
  //   motor_pos_zero(13);
  //   delay(50);
  //  motor_pos_zero(21);
  // delay(50);
  // motor_pos_zero(22);
  // delay(50);
  // motor_pos_zero(23);k
  // delay(50);
  motor_pos_zero(31);
  delay(50);
  motor_pos_zero(32);
  delay(50);
  motor_pos_zero(33);
  delay(50);
  motor_pos_zero(41);
  delay(50);
  motor_pos_zero(42);
  delay(50);
  motor_pos_zero(43);
  delay(50);

  //  使能id电机
  // motor_enable(11);
  // delay(50);
  // motor_enable(12);
  // delay(50);
  // motor_enable(13);
  // delay(50);
  // motor_enable(21);
  // delay(50);
  // motor_enable(22);
  // delay(50);
  // motor_enable(23);
  // delay(50);
  motor_enable(31);
  delay(50);
  motor_enable(32);
  delay(50);
  motor_enable(33);
  delay(50);
  motor_enable(41);
  delay(50);
  motor_enable(42);
  delay(50);
  motor_enable(43);
  delay(50);

  // 电机运行模式 电机canid 模式值 1位置模式2速度模式 3 电流模式0运控模式
  //      motor_mode(13, 1);
  //      delay(50);
  //      motor_mode(11, 1);
  //      delay(50);
  //      motor_mode(12, 1);
  //      delay(50);
  //   motor_mode(23, 1);
  //   delay(50);
  //   motor_mode(21, 1);
  //   delay(50);
  //   motor_mode(22, 1);
  //   delay(50);
  motor_mode(31, 1);
  delay(50);
  motor_mode(32, 1);
  delay(50);
  motor_mode(33, 1);
  delay(50);
  motor_mode(41, 1);
  delay(50);
  motor_mode(42, 1);
  delay(50);
  motor_mode(43, 1);
  delay(50);

  // motor_pow_cmd(11, SPD_LIMIT, CUR_LIMIT, KP, KI);  //id, speed_lim, current_lim, kp, ki
  // delay(50);
  // motor_pow_cmd(12, SPD_LIMIT, CUR_LIMIT, KP, KI);
  // delay(50);
  // motor_pow_cmd(13, SPD_LIMIT, CUR_LIMIT, KP, KI);
  // delay(50);
  // motor_pow_cmd(21, SPD_LIMIT, CUR_LIMIT, KP, KI);
  // delay(50);
  // motor_pow_cmd(22, SPD_LIMIT, CUR_LIMIT, KP, KI);
  // delay(50);
  // motor_pow_cmd(23, SPD_LIMIT, CUR_LIMIT, KP, KI);
  // delay(50);
  motor_pow_cmd(31, SPD_LIMIT, CUR_LIMIT, KP, KI);
  delay(50);
  motor_pow_cmd(32, SPD_LIMIT, CUR_LIMIT, KP, KI);
  delay(50);
  motor_pow_cmd(33, SPD_LIMIT, CUR_LIMIT, KP, KI);
  delay(50);
  motor_pow_cmd(41, SPD_LIMIT, CUR_LIMIT, KP, KI);
  delay(50);
  motor_pow_cmd(42, SPD_LIMIT, CUR_LIMIT, KP, KI);
  delay(50);
  motor_pow_cmd(43, SPD_LIMIT, CUR_LIMIT, KP, KI);
  delay(50);

  Serial.println("Setup done");
  delay(2000);
}

BLA::Matrix<3, 3> times3(double t)
{
  BLA::Matrix<3, 3> I;
  I = {t, 0, 0,
       0, t, 0,
       0, 0, t};
  return I;
}

// the function below is used to generate the 4x4 identity matrix multiplied by a scalar
BLA::Matrix<4, 4> times4(double t)
{
  BLA::Matrix<4, 4> I;
  I = {t, 0, 0, 0,
       0, t, 0, 0,
       0, 0, t, 0,
       0, 0, 0, t};
  return I;
}

// the function below is used to generate the 3x3 skew symmetric matrix of a 3x1 vector
BLA::Matrix<3, 3> SkewSemmetric(BLA::Matrix<3, 1> v)
{
  BLA::Matrix<3, 3> S;
  S(0, 0) = 0;
  S(0, 1) = -v(2);
  S(0, 2) = v(1);
  S(1, 0) = v(2);
  S(1, 1) = 0;
  S(1, 2) = -v(0);
  S(2, 0) = -v(1);
  S(2, 1) = v(0);
  S(2, 2) = 0;
  return S;
}

// the function below is used to generate the 3x3 exponential matrix of a 3x1 vector and an angle scalar
BLA::Matrix<3, 3> expm3(BLA::Matrix<3, 1> w, double t)
{
  BLA::Matrix<3, 3> S = SkewSemmetric(w);
  BLA::Matrix<3, 3> E;
  E = times3(1.0) + times3(sin(t)) * S + times3(1.0 - cos(t)) * S * S;
  return E;
}

// the function below is used to generate the 4x4 screw matrix of a 6x1 screw vector
BLA::Matrix<4, 4> Screw44(BLA::Matrix<6, 1> q)
{
  BLA::Matrix<3, 3> w;
  BLA::Matrix<3, 1> v;
  BLA::Matrix<1, 4> Z;
  w = SkewSemmetric(q.Submatrix<3, 1>(0, 0));
  v = q.Submatrix<3, 1>(3, 0);
  Z.Fill(0);
  return (w || v) && Z;
}

// the function below is used to generate the 4x4 exponential matrix of a 6x1 screw vector and an angle scalar
BLA::Matrix<4, 4> expm4(BLA::Matrix<6, 1> S, double t)
{
  BLA::Matrix<3, 1> w = S.Submatrix<3, 1>(0, 0);
  BLA::Matrix<3, 1> v = S.Submatrix<3, 1>(3, 0);
  BLA::Matrix<3, 3> W = expm3(w, t);
  BLA::Matrix<3, 3> SkW = SkewSemmetric(w);
  BLA::Matrix<3, 1> p = (times3(t) + times3(1 - cos(t)) * SkW + times3(t - sin(t)) * SkW * SkW) * v;
  BLA::Matrix<1, 4> Z = {0, 0, 0, 1};
  BLA::Matrix<4, 4> E = (W || p) && Z;
  return E;
}

// the function below is used to generate
BLA::Matrix<6, 6> Adjoint(BLA::Matrix<4, 4> T)
{
  BLA::Matrix<3, 3> R = T.Submatrix<3, 3>(0, 0);
  BLA::Matrix<3, 1> p = T.Submatrix<3, 1>(0, 3);
  BLA::Matrix<3, 3> p_hat = SkewSemmetric(p);
  BLA::Matrix<6, 6> Ad;
  Ad = (R || times3(0.0)) && (p_hat * R || R);
  return Ad;
}

BLA::Matrix<3, 1> IK(BLA::Matrix<3, 1> p)
{
  BLA::Matrix<3, 1> q;
  double a5 = 2 * l1 * (l2 + l4) * sin(bt);
  double b5 = 2 * l3 * (l2 + l4);
  double c5 = l1 * l1 + l3 * l3 + (l2 + l4) * (l2 + l4) - p(0) * p(0) - p(1) * p(1) - p(2) * p(2);

  double q51 = atan2(-c5, sqrt(a5 * a5 + b5 * b5 - c5 * c5)) - atan2(b5, a5);
  double q52 = atan2(-c5, -sqrt(a5 * a5 + b5 * b5 - c5 * c5)) - atan2(b5, a5);
  if (2 * fabs(q51) < PI)
  {
    q(2) = q51;
  }
  else if (2 * fabs(q52) < PI)
  {
    q(2) = q52;
  }
  else
  {
    Serial.println("Pz No solution");
    ErrorFlag = 1;
  }

  double a4 = l3 + (l2 + l4) * cos(q(2));
  double b4 = (l2 + l4) * cos(bt) * sin(q(2));
  double c4 = p(0);

  double q41 = atan2(-c4, sqrt(a4 * a4 + b4 * b4 - c4 * c4)) - atan2(b4, a4);
  double q42 = atan2(-c4, -sqrt(a4 * a4 + b4 * b4 - c4 * c4)) - atan2(b4, a4);

  if (2 * fabs(q41) < PI)
  {
    q(1) = q41;
  }
  else if (2 * fabs(q42) < PI)
  {
    q(1) = q42;
  }
  else
  {
    Serial.println("Py No solution");
    ErrorFlag = 1;
  }

  double a1 = -(l2 + l4) * cos(bt) * sin(q(1)) * sin(q(2)) + (l2 + l4) * cos(q(1)) * cos(q(2)) + l3 * cos(q(1));
  double b1 = l1 + (l2 + l4) * sin(bt) * sin(q(2));

  q(0) = asin((p(1) * a1 + p(2) * b1) / (p(1) * p(1) + p(2) * p(2)));

  return q;
}

BLA::Matrix<6, 3> GeometricalJacobian(double q1, double q4)
{
  BLA::Matrix<6, 1> Jg2 = Adjoint(expm4(S1, q1)) * S2;
  BLA::Matrix<6, 1> Jg3 = Adjoint(expm4(S1, q1) * expm4(S2, q4)) * S3;
  BLA::Matrix<6, 3> Jg = (S1 || Jg2 || Jg3);
  return Jg;
}

BLA::Matrix<3, 3> AnalyticalJacobian(double q1, double q4, BLA::Matrix<3, 1> p)
{
  BLA::Matrix<3, 6> J = (times3(-1.0) * SkewSemmetric(p)) || (times3(1.0));
  BLA::Matrix<3, 3> Ja = J * GeometricalJacobian(q1, q4);
  return Ja;
}

float t = 0.0;
float pos = 0;
static long last = 0;
long now = 0;

BLA::Matrix<3> FL = {0.0, 0.0, 0.0};
BLA::Matrix<3> FR = {0.0, 0.0, 0.0};
BLA::Matrix<3> RL = {0.0, 0.0, 0.0};
BLA::Matrix<3> RR = {0.0, 0.0, 0.0};

void loop()
{
  pos = 0.07 * sin(t);
  t = t + 0.008;
  BLA::Matrix<3> p;
  p(0) = 0.0;
  p(1) = 0.23;
  p(2) = -0.17;

  BLA::Matrix<3> I = {0, 0, pos};
  //  I.Fill(pos);

  //  BLA::Matrix<3> q = IK(p + I) * RAD2DEG;
  RL = IK(p + I);
  RR = RL;
  FL = RL;
  FR = RL;
  // Serial << "RL: " << RL << '\n';
  Serial.print("RL1: ");
  Serial.print(RL(0));
  Serial.print(" ");
  Serial.print("RL2: ");
  Serial.print(RL(1));
  Serial.print(" ");
   Serial.print("RL3: ");
  Serial.print(RL(2));
  Serial.println(" ");

  long now = micros();
  // Serial.println(now - last);
  last = now;
  if (ErrorFlag == 0)
  {
    motor_pos_cmd(11, FL(0) * Sgn11);
    delayMicroseconds(DelayMS);
    motor_pos_cmd(12, FL(1) * Sgn12);
    delayMicroseconds(DelayMS);
    if (FL(2) > -1.2 - Bias13 && FL(2) < -Bias13)
    {
      motor_pos_cmd(13, Sgn13 * (FL(2) + Bias13));
      delayMicroseconds(DelayMS);
    }
    else
    {
      Serial.println("q13 out of range");
    }
    motor_pos_cmd(21, FR(0) * Sgn21);
    delayMicroseconds(DelayMS);
    motor_pos_cmd(22, FR(1) * Sgn22);
    delayMicroseconds(DelayMS);
    if (FR(2) < -Bias23 && FR(2) > -Bias23 - 1.2)
    {
      motor_pos_cmd(23, Sgn23 * (FR(2) + Bias23));
      delayMicroseconds(DelayMS);
    }
    else
    {
      Serial.println("q23 out of range");
    }
    motor_pos_cmd(31, RL(0) * Sgn31);
    delayMicroseconds(DelayMS);
    motor_pos_cmd(32, RL(1) * Sgn32);
    delayMicroseconds(DelayMS);
    if (RL(2) < -Bias33 && RL(2) > -Bias33 - 1.2)
    {
      motor_pos_cmd(33, Sgn33 * (RL(2) + Bias33));
      delayMicroseconds(DelayMS);
    }
    else
    {
      Serial.println("q33 out of range");
    }
    motor_pos_cmd(41, RR(0) * Sgn41);
    delayMicroseconds(DelayMS);
    motor_pos_cmd(42, RR(1) * Sgn42);
    delayMicroseconds(DelayMS);
    if (RR(2) > -1.2 - Bias43 && RR(2) < -Bias43)
    {
      motor_pos_cmd(43, Sgn43 * (RR(2) + Bias43));
      delayMicroseconds(DelayMS);
    }
    else
    {
      Serial.println("q43 out of range");
    }
  }

  //  motor_speed_cmd(1,speed_valueb);
  //  Serial.print(pos);
  //  Serial.print("\t");
  //  now = micros();
  //  Serial.println(now - last);

  //  delayMicroseconds(1000+last-now);
}
