#ifndef RK4_H
#define RK4_H

#include <Arduino.h>

// Generic derivative callback:
// t    : current time
// x    : current state vector, length n
// dxdt : derivative output vector, length n
// ctx  : user context / model parameters
typedef void (*DerivFunc)(float t, const float* x, float* dxdt, void* ctx);

// Small fixed-size workspace version for MCUs.
// Increase RK4_MAX_N if you need more states.
#ifndef RK4_MAX_N
#define RK4_MAX_N 8
#endif

struct RK4Solver {
  DerivFunc f = nullptr;
  void* ctx = nullptr;
  int n = 0;

  bool begin(int n_, DerivFunc f_, void* ctx_) {
    if (n_ <= 0 || n_ > RK4_MAX_N || f_ == nullptr) return false;
    n = n_;
    f = f_;
    ctx = ctx_;
    return true;
  }

  bool step(float t, float h, const float* x, float* xout) {
    if (!f || n <= 0 || n > RK4_MAX_N) return false;

    float k1[RK4_MAX_N];
    float k2[RK4_MAX_N];
    float k3[RK4_MAX_N];
    float k4[RK4_MAX_N];
    float xtmp[RK4_MAX_N];

    f(t, x, k1, ctx);

    for (int i = 0; i < n; i++) {
      xtmp[i] = x[i] + 0.5f * h * k1[i];
    }
    f(t + 0.5f * h, xtmp, k2, ctx);

    for (int i = 0; i < n; i++) {
      xtmp[i] = x[i] + 0.5f * h * k2[i];
    }
    f(t + 0.5f * h, xtmp, k3, ctx);

    for (int i = 0; i < n; i++) {
      xtmp[i] = x[i] + h * k3[i];
    }
    f(t + h, xtmp, k4, ctx);

    const float h6 = h / 6.0f;
    for (int i = 0; i < n; i++) {
      xout[i] = x[i] + h6 * (k1[i] + 2.0f * k2[i] + 2.0f * k3[i] + k4[i]);
    }

    return true;
  }
};

#endif