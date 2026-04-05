#include <Arduino.h>
#include <math.h>

#ifndef DAC_PIN
#define DAC_PIN PA5
#endif

#ifndef ADC_PIN
#define ADC_PIN A0
#endif

struct RDParams {
  float ecg_fs;
  float internal_fs;
  float duration_s;

  float H;
  float HR_bpm;
  float C;
  float beta;

  float alpha1;
  float alpha2;
  float alpha3;
  float alpha4;

  float x1_0;
  float x2_0;
  float x3_0;
  float x4_0;

  int newton_iters;
};

struct RDContext {
  float dt;
  float Gamma_t;
  float x[4];
};

static bool solve4x4(const float Ain[4][4], const float bin[4], float xout[4]) {
  float a[4][5];

  for (int r = 0; r < 4; r++) {
    for (int c = 0; c < 4; c++) a[r][c] = Ain[r][c];
    a[r][4] = bin[r];
  }

  for (int col = 0; col < 4; col++) {
    int piv = col;
    float vmax = fabsf(a[col][col]);
    for (int r = col + 1; r < 4; r++) {
      float v = fabsf(a[r][col]);
      if (v > vmax) {
        vmax = v;
        piv = r;
      }
    }
    if (vmax < 1e-12f) return false;

    if (piv != col) {
      for (int c = col; c < 5; c++) {
        float tmp = a[col][c];
        a[col][c] = a[piv][c];
        a[piv][c] = tmp;
      }
    }

    float div = a[col][col];
    for (int c = col; c < 5; c++) a[col][c] /= div;

    for (int r = 0; r < 4; r++) {
      if (r == col) continue;
      float f = a[r][col];
      for (int c = col; c < 5; c++) a[r][c] -= f * a[col][c];
    }
  }

  for (int i = 0; i < 4; i++) xout[i] = a[i][4];
  return true;
}

static inline void rd_f_raw(const RDParams& p, const float x[4], float f[4]) {
  const float x1 = x[0];
  const float x2 = x[1];
  const float x3 = x[2];
  const float x4 = x[3];

  f[0] = x1 - x2 - p.C * x1 * x2 - x1 * x2 * x2;
  f[1] = p.H * x1 - 3.0f * x2 + p.C * x1 * x2 + x1 * x2 * x2 + p.beta * (x4 - x2);
  f[2] = x3 - x4 - p.C * x3 * x4 - x3 * x4 * x4;
  f[3] = p.H * x3 - 3.0f * x4 + p.C * x3 * x4 + x3 * x4 * x4 + 2.0f * p.beta * (x2 - x4);
}

static inline void rd_f_scaled(const RDParams& p, const RDContext& ctx, const float x[4], float f[4]) {
  rd_f_raw(p, x, f);
  for (int i = 0; i < 4; i++) f[i] *= ctx.Gamma_t;
}

static inline void rd_df_raw(const RDParams& p, const float x[4], float J[4][4]) {
  const float x1 = x[0];
  const float x2 = x[1];
  const float x3 = x[2];
  const float x4 = x[3];

  J[0][0] =  1.0f - p.C * x2 - x2 * x2;
  J[0][1] = -1.0f - p.C * x1 - 2.0f * x1 * x2;
  J[0][2] =  0.0f;
  J[0][3] =  0.0f;

  J[1][0] =  p.H + p.C * x2 + x2 * x2;
  J[1][1] = -3.0f + p.C * x1 + 2.0f * x1 * x2 - p.beta;
  J[1][2] =  0.0f;
  J[1][3] =  p.beta;

  J[2][0] =  0.0f;
  J[2][1] =  0.0f;
  J[2][2] =  1.0f - p.C * x4 - x4 * x4;
  J[2][3] = -1.0f - p.C * x3 - 2.0f * x3 * x4;

  J[3][0] =  0.0f;
  J[3][1] =  2.0f * p.beta;
  J[3][2] =  p.H + p.C * x4 + x4 * x4;
  J[3][3] = -3.0f + p.C * x3 + 2.0f * x3 * x4 - 2.0f * p.beta;
}

static bool rd_tustin_step(const RDParams& p, const RDContext& ctx,
                           const float x_prev[4], float x_next[4]) {
  float xn[4] = { x_prev[0], x_prev[1], x_prev[2], x_prev[3] };
  float fp[4];
  rd_f_scaled(p, ctx, x_prev, fp);

  for (int iter = 0; iter < p.newton_iters; iter++) {
    float fn[4];
    rd_f_scaled(p, ctx, xn, fn);

    float R[4];
    for (int i = 0; i < 4; i++) {
      R[i] = xn[i] - x_prev[i] - 0.5f * ctx.dt * (fn[i] + fp[i]);
    }

    float Df[4][4];
    rd_df_raw(p, xn, Df);

    float J[4][4];
    for (int r = 0; r < 4; r++) {
      for (int c = 0; c < 4; c++) {
        J[r][c] = -0.5f * ctx.dt * ctx.Gamma_t * Df[r][c];
      }
      J[r][r] += 1.0f;
    }

    float dx[4];
    if (!solve4x4(J, R, dx)) return false;

    float max_abs = 0.0f;
    for (int i = 0; i < 4; i++) {
      xn[i] -= dx[i];
      float a = fabsf(dx[i]);
      if (a > max_abs) max_abs = a;
    }

    if (max_abs < 1e-6f) break;
  }

  for (int i = 0; i < 4; i++) x_next[i] = xn[i];
  return true;
}

static inline float rd_ecg_mv(const RDParams& p, const float x[4]) {
  return p.alpha1 * x[0]
       + p.alpha2 * x[1]
       + p.alpha3 * x[2]
       + p.alpha4 * x[3];
}

static bool rd_init(const RDParams& p, RDContext& ctx) {
  if (p.internal_fs <= 0.0f || p.ecg_fs <= 0.0f || p.duration_s <= 0.0f) return false;
  if (fabsf(p.internal_fs - p.ecg_fs) > 1e-6f) return false;

  ctx.dt = 1.0f / p.internal_fs;
  ctx.Gamma_t = 0.08804f * p.HR_bpm - 0.06754f;
  ctx.x[0] = p.x1_0;
  ctx.x[1] = p.x2_0;
  ctx.x[2] = p.x3_0;
  ctx.x[3] = p.x4_0;
  return true;
}

static bool build_rd_block_mv(const RDParams& p, float*& out_mv, int& out_len, RDContext& ctx) {
  out_len = (int)lroundf(p.duration_s * p.ecg_fs);
  if (out_len <= 0) return false;

  out_mv = (float*)malloc((size_t)out_len * sizeof(float));
  if (!out_mv) return false;

  for (int n = 0; n < out_len; n++) {
    float x_next[4];
    if (!rd_tustin_step(p, ctx, ctx.x, x_next)) {
      free(out_mv);
      out_mv = nullptr;
      out_len = 0;
      return false;
    }

    for (int i = 0; i < 4; i++) ctx.x[i] = x_next[i];
    out_mv[n] = rd_ecg_mv(p, ctx.x);
  }

  return true;
}

static uint16_t ecg_mv_to_dac12(float mv, float gain = 1200.0f, float offset_mv = 1700.0f) {
  float v = offset_mv + gain * mv;
  if (v < 0.0f) v = 0.0f;
  if (v > 3300.0f) v = 3300.0f;
  return (uint16_t)lroundf((v / 3300.0f) * 4095.0f);
}

static RDParams g_p;
static RDContext g_ctx;
static float* g_block = nullptr;
static int g_block_len = 0;
volatile int playIndex = 0;

void setup() {
  Serial.begin(921600);
  delay(3000);

  analogWriteResolution(12);
  analogReadResolution(12);
  pinMode(DAC_PIN, OUTPUT);
  pinMode(ADC_PIN, INPUT);

  g_p.ecg_fs = 200.0f;
  g_p.internal_fs = 200.0f;
  g_p.duration_s = 30.0f;

  g_p.H = 3.0f;
  g_p.HR_bpm = 95.0f;
  g_p.C = 1.35f;
  g_p.beta = 4.0f;

  g_p.alpha1 = -0.024f;
  g_p.alpha2 =  0.0216f;
  g_p.alpha3 = -0.0012f;
  g_p.alpha4 =  0.12f;

  g_p.x1_0 = 0.0f;
  g_p.x2_0 = 0.0f;
  g_p.x3_0 = 0.1f;
  g_p.x4_0 = 0.0f;

  g_p.newton_iters = 4;

  Serial.println("Init RD ECG...");
  if (!rd_init(g_p, g_ctx)) {
    Serial.println("rd_init failed");
    while (1) delay(1000);
  }

  uint32_t t0 = micros();
  if (!build_rd_block_mv(g_p, g_block, g_block_len, g_ctx)) {
    Serial.println("build_rd_block_mv failed");
    while (1) delay(1000);
  }
  uint32_t dt_us = micros() - t0;

  Serial.print("Gamma_t = ");
  Serial.println(g_ctx.Gamma_t, 6);
  Serial.print("block_len = ");
  Serial.println(g_block_len);
  Serial.print("build_us = ");
  Serial.println(dt_us);

  Serial.println("Playback start.");
  Serial.println("idx,ecg_mv,dac_code,adc_code");
}

void loop() {
  static uint32_t last_us = 0;
  const uint32_t period_us = (uint32_t)lroundf(1000000.0f / g_p.ecg_fs);
  uint32_t now = micros();

  if ((uint32_t)(now - last_us) >= period_us) {
    last_us += period_us;

    const int idx = playIndex;
    const float ecg_mv = g_block[idx];
    const uint16_t dac_code = ecg_mv_to_dac12(ecg_mv);

    analogWrite(DAC_PIN, dac_code);
    delayMicroseconds(20);
    int adc_code = analogRead(ADC_PIN);

    Serial.print(idx);
    Serial.print(",");
    Serial.print(ecg_mv, 6);
    Serial.print(",");
    Serial.print(dac_code);
    Serial.print(",");
    Serial.println(adc_code);

    playIndex++;
    if (playIndex >= g_block_len) playIndex = 0;
  }
}