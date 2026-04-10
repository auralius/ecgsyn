#include "ecgsyn_model.h"
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#ifndef PI_D
#define PI_D 3.1415926535897932384626433832795
#endif

#define USE_TUSTIN
//#define USE_RK4

static EcgSynContext* g_ctx = NULL;

// ------------------------------------------------------------
// 1-based vector helpers
// ------------------------------------------------------------
static inline float* mallocVect(long n0, long nx)
{
    const long OFFSET = 1;
    float *base = (float *)malloc((size_t)(nx - n0 + 1 + OFFSET) * sizeof(float));
    if (base == NULL) {
        return NULL;
    }
    return base - n0 + OFFSET;
}

static inline void freeVect(float* vect, long n0, long nx)
{
    (void)nx;
    const long OFFSET = 1;
    if (vect != NULL) {
        free((void *)(vect + n0 - OFFSET));
    }
}

static inline float stdev(float* x, int n) {
  if (n <= 1) return 0.0;

  float add = 0.0;
  for (int j = 1; j <= n; j++) add += x[j];
  float mean = add / n;

  float total = 0.0;
  for (int j = 1; j <= n; j++) {
    float diff = x[j] - mean;
    total += diff * diff;
  }
  return sqrtf(total / (n - 1));
}

// ------------------------------------------------------------
// Numerical Recipes ran1
// ------------------------------------------------------------
static inline float ran1(long *idum)
{
#define IA    16807L
#define IM    2147483647L
#define IQ    127773L
#define IR    2836L
#define NTAB  32
#define NDIV  (1L + (IM - 1L) / NTAB)
#define EPS   1.2e-7f
#define RNMX  (1.0f - EPS)
#define AM    (1.0f / (float)IM)

    int j;
    long k;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum = 1;
        else *idum = -(*idum);

        for (j = NTAB + 7; j >= 0; j--) {
            k = (*idum) / IQ;
            *idum = IA * (*idum - k * IQ) - IR * k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy = iv[0];
    }

    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if (*idum < 0) *idum += IM;
    j = (int)(iy / NDIV);
    iy = iv[j];
    iv[j] = *idum;
    temp = AM * (float)iy;

    return (temp > RNMX) ? RNMX : temp;

#undef IA
#undef IM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
#undef AM
}

// ------------------------------------------------------------
// Numerical Recipes dfour1
// ------------------------------------------------------------
static inline void dfour1(float data[], unsigned long nn, int isign) {
  unsigned long n, mmax, m, j, istep, i;
  float wtemp, wr, wpr, wpi, wi, theta;
  float tempr, tempi;

  n = nn << 1;
  j = 1;

  for (i = 1; i < n; i += 2) {
    if (j > i) {
      tempr = data[j];     data[j]     = data[i];     data[i]     = tempr;
      tempr = data[j + 1]; data[j + 1] = data[i + 1]; data[i + 1] = tempr;
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }

  mmax = 2;
  while (n > mmax) {
    istep = mmax << 1;
    theta = isign * (2.0 * PI_D / mmax);
    wtemp = sinf(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sinf(theta);
    wr = 1.0;
    wi = 0.0;

    for (m = 1; m < mmax; m += 2) {
      for (i = m; i <= n; i += istep) {
        j = i + mmax;
        tempr = wr * data[j]     - wi * data[j + 1];
        tempi = wr * data[j + 1] + wi * data[j];
        data[j]     = data[i]     - tempr;
        data[j + 1] = data[i + 1] - tempi;
        data[i]     += tempr;
        data[i + 1] += tempi;
      }
      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
    }
    mmax = istep;
  }
}

// ------------------------------------------------------------
// RR process
// ------------------------------------------------------------
static inline bool rrprocess(
    float* rr,
    float flo, float fhi,
    float flostd, float fhistd, float lfhfratio,
    float hrmean, float hrstd,
    float sf, int n,
    long* rseed)
{
  int i;
  float c1, c2, w1, w2, sig1, sig2, rrmean, rrstd, xstd, ratio;
  float df, dw1, dw2;
  float *w, *Hw, *Sw, *ph0, *ph, *SwC;

  w   = mallocVect(1, n);
  Hw  = mallocVect(1, n);
  Sw  = mallocVect(1, n);
  ph0 = mallocVect(1, n / 2 - 1);
  ph  = mallocVect(1, n);
  SwC = mallocVect(1, 2 * n);

  if (!w || !Hw || !Sw || !ph0 || !ph || !SwC) {
    freeVect(w, 1, n);
    freeVect(Hw, 1, n);
    freeVect(Sw, 1, n);
    freeVect(ph0, 1, n / 2 - 1);
    freeVect(ph, 1, n);
    freeVect(SwC, 1, 2 * n);
    return false;
  }

  w1 = 2.0 * PI_D * flo;
  w2 = 2.0 * PI_D * fhi;
  c1 = 2.0 * PI_D * flostd;
  c2 = 2.0 * PI_D * fhistd;
  sig2 = 1.0;
  sig1 = lfhfratio;
  rrmean = 60.0 / hrmean;
  rrstd  = 60.0 * hrstd / (hrmean * hrmean);

  df = sf / n;
  for (i = 1; i <= n; i++) w[i] = (i - 1) * 2.0 * PI_D * df;

  for (i = 1; i <= n; i++) {
    dw1 = w[i] - w1;
    dw2 = w[i] - w2;
    Hw[i] = sig1 * expf(-dw1 * dw1 / (2.0 * c1 * c1)) / sqrtf(2.0 * PI_D * c1 * c1)
          + sig2 * expf(-dw2 * dw2 / (2.0 * c2 * c2)) / sqrtf(2.0 * PI_D * c2 * c2);
  }

  for (i = 1; i <= n / 2; i++)     Sw[i] = (sf / 2.0) * sqrtf(Hw[i]);
  for (i = n / 2 + 1; i <= n; i++) Sw[i] = (sf / 2.0) * sqrtf(Hw[n - i + 1]);

  for (i = 1; i <= n / 2 - 1; i++) ph0[i] = 2.0 * PI_D * ran1(rseed);
  ph[1] = 0.0;
  for (i = 1; i <= n / 2 - 1; i++) ph[i + 1] = ph0[i];
  ph[n / 2 + 1] = 0.0;
  for (i = 1; i <= n / 2 - 1; i++) ph[n - i + 1] = -ph0[i];

  for (i = 1; i <= n; i++) SwC[2 * i - 1] = Sw[i] * cosf(ph[i]);
  for (i = 1; i <= n; i++) SwC[2 * i]     = Sw[i] * sinf(ph[i]);

  dfour1(SwC, n, -1);

  for (i = 1; i <= n; i++) rr[i] = (1.0 / n) * SwC[2 * i - 1];

  xstd = stdev(rr, n);
  if (xstd < 1e-12) xstd = 1.0;
  ratio = rrstd / xstd;

  for (i = 1; i <= n; i++) rr[i] *= ratio;
  for (i = 1; i <= n; i++) rr[i] += rrmean;

  freeVect(w, 1, n);
  freeVect(Hw, 1, n);
  freeVect(Sw, 1, n);
  freeVect(ph0, 1, n / 2 - 1);
  freeVect(ph, 1, n);
  freeVect(SwC, 1, 2 * n);
  return true;
}

// ------------------------------------------------------------
// ODE pieces
// ------------------------------------------------------------
static inline float angfreq(float t) {
  int i = 1 + (int)floorf(t / g_ctx->h);
  return 2.0 * PI_D / g_ctx->rrpc[i];
}

static inline void derivspqrst(float t0, float x[], float dxdt[]) {
  int i, k;
  float a0, w0, r0, x0, y0;
  float t, dt, dt2, zbase;

  k = 5;
  w0 = angfreq(t0);
  r0 = 1.0; x0 = 0.0; y0 = 0.0;
  a0 = 1.0 - sqrtf((x[1] - x0) * (x[1] - x0) + (x[2] - y0) * (x[2] - y0)) / r0;

  zbase = 0.005 * sinf(2.0 * PI_D * g_ctx->fhi * t0);

  t = atan2f(x[2], x[1]);
  dxdt[1] = a0 * (x[1] - x0) - w0 * (x[2] - y0);
  dxdt[2] = a0 * (x[2] - y0) + w0 * (x[1] - x0);
  dxdt[3] = 0.0;

  for (i = 1; i <= k; i++) {
    dt = fmodf(t - g_ctx->ti[i], 2.0 * PI_D);
    dt2 = dt * dt;
    dxdt[3] += -g_ctx->ai[i] * dt * expf(-0.5 * dt2 / (g_ctx->bi[i] * g_ctx->bi[i]));
  }
  dxdt[3] += -1.0 * (x[3] - zbase);
}

// Removed unnecessary malloc/free here
static inline void drk4(float y[], int n, float x, float h, float yout[]) {
  if (n != 3) return;

  float dydx[4];
  float dym[4];
  float dyt[4];
  float yt[4];

  float hh = 0.5 * h;
  float h6 = h / 6.0;
  float xh = x + hh;

  derivspqrst(x, y, dydx);
  for (int i = 1; i <= 3; i++) yt[i] = y[i] + hh * dydx[i];

  derivspqrst(xh, yt, dyt);
  for (int i = 1; i <= 3; i++) yt[i] = y[i] + hh * dyt[i];

  derivspqrst(xh, yt, dym);
  for (int i = 1; i <= 3; i++) {
    yt[i] = y[i] + h * dym[i];
    dym[i] += dyt[i];
  }

  derivspqrst(x + h, yt, dyt);
  for (int i = 1; i <= 3; i++) {
    yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
  }
}

// ------------------------------------------------------------
// Build ECG block
// ------------------------------------------------------------
bool build_block_mv(
    const EcgSynParams *p,
    float **out_mv,
    int *out_len,
    EcgSynContext *ctx)
{
    *out_mv = NULL;
    *out_len = 0;

    if (p == NULL || out_mv == NULL || out_len == NULL || ctx == NULL) return false;
    if (p->ecg_fs <= 0 || p->internal_fs <= 0 || p->n_beats <= 0) return false;
    if ((p->internal_fs % p->ecg_fs) != 0) return false;

    const int q = p->internal_fs / p->ecg_fs;
    if (q <= 0) return false;

    ctx->h = 1.0f / (float)p->internal_fs;
    ctx->fhi = p->fhi;
    ctx->rseed = -labs((long)p->seed_init);

    ctx->ti[1] = -60.0f; ctx->ti[2] = -15.0f; ctx->ti[3] = 0.0f;  ctx->ti[4] = 15.0f; ctx->ti[5] = 90.0f;
    ctx->ai[1] = 1.2f;   ctx->ai[2] = -5.0f;  ctx->ai[3] = 30.0f; ctx->ai[4] = -7.5f; ctx->ai[5] = 0.75f;
    ctx->bi[1] = 0.25f;  ctx->bi[2] = 0.1f;   ctx->bi[3] = 0.1f;  ctx->bi[4] = 0.1f;  ctx->bi[5] = 0.4f;

    for (int i = 1; i <= 5; i++) ctx->ti[i] *= (float)(PI_D / 180.0);

    float hrfact = sqrtf(p->hr_mean / 60.0f);
    float hrfact2 = sqrtf(hrfact);
    for (int i = 1; i <= 5; i++) ctx->bi[i] *= hrfact;
    ctx->ti[1] *= hrfact2;
    ctx->ti[2] *= hrfact;
    ctx->ti[3] *= 1.0f;
    ctx->ti[4] *= hrfact;
    ctx->ti[5] *= 1.0f;

    float rrmean = 60.0f / p->hr_mean;
    int Nrr = 1;
    while (Nrr < (int)(p->n_beats * rrmean * p->internal_fs)) Nrr <<= 1;

    ctx->rr = mallocVect(1, Nrr);
    if (!ctx->rr) return false;

    if (!rrprocess(ctx->rr, p->flo, p->fhi, p->flo_std, p->fhi_std,
                   p->lf_hf_ratio, p->hr_mean, p->hr_std,
                   (float)p->internal_fs, Nrr, &ctx->rseed)) {
        return false;
    }

    float tecg = 0.0f;
    for (int beat = 1; beat <= p->n_beats; beat++) tecg += ctx->rr[beat];

    int Nt = (int)llround(tecg / ctx->h);
    if (Nt <= 0) return false;

    ctx->rrpc = mallocVect(1, Nt);
    if (!ctx->rrpc) return false;

    tecg = 0.0f;
    int sample_idx = 1;
    for (int beat = 1; beat <= p->n_beats; beat++) {
        tecg += ctx->rr[beat];
        int sample_end = (int)llround(tecg / ctx->h);
        if (sample_end > Nt) sample_end = Nt;
        for (int k = sample_idx; k <= sample_end; k++) ctx->rrpc[k] = ctx->rr[beat];
        sample_idx = sample_end + 1;
    }
    while (sample_idx <= Nt) {
        ctx->rrpc[sample_idx] = ctx->rr[p->n_beats];
        sample_idx++;
    }

    float x[4];
    float xn[4];
    float *zt = mallocVect(1, Nt);
    if (!zt) return false;

    x[1] = p->xinitial;
    x[2] = p->yinitial;
    x[3] = p->zinitial;

    g_ctx = ctx;
    float timev = 0.0f;
#ifdef USE_RK4
    for (int i = 1; i <= Nt; i++) {
        zt[i] = x[3];
        drk4(x, 3, timev, ctx->h, xn);
        x[1] = xn[1];
        x[2] = xn[2];
        x[3] = xn[3];
        timev += ctx->h;
    }
#elif defined(USE_TUSTIN)
    for (int i = 1; i <= Nt; i++) {
        zt[i] = x[3];

        if (!implicit_tustin_step(x, timev, ctx->h, xn)) {
            freeVect(zt, 1, Nt);
            g_ctx = NULL;
            return false;
        }

        x[1] = xn[1];
        x[2] = xn[2];
        x[3] = xn[3];
        timev += ctx->h;
    }
#endif
    g_ctx = NULL;

    *out_len = (Nt + q - 1) / q;
    if (*out_len <= 0) {
        freeVect(zt, 1, Nt);
        return false;
    }

    *out_mv = (float *)malloc((size_t)(*out_len) * sizeof(float));
    if (!(*out_mv)) {
        freeVect(zt, 1, Nt);
        return false;
    }

    float zmin = zt[1], zmax = zt[1];
    for (int i = 2; i <= Nt; i++) {
        if (zt[i] < zmin) zmin = zt[i];
        if (zt[i] > zmax) zmax = zt[i];
    }

    float zrange = zmax - zmin;
    if (zrange < 1e-12f) zrange = 1.0f;

    int j = 0;
    for (int i = 1; i <= Nt; i += q) {
        float z = (zt[i] - zmin) * 1.6f / zrange - 0.4f;
        z += p->noise_mv * (2.0f * ran1(&ctx->rseed) - 1.0f);
        (*out_mv)[j++] = z;
    }

    freeVect(zt, 1, Nt);
    return true;
}

void free_context(EcgSynContext *ctx)
{
    if (ctx == NULL) return;

    if (ctx->rr) {
        freeVect(ctx->rr, 1, 1);
        ctx->rr = NULL;
    }
    if (ctx->rrpc) {
        freeVect(ctx->rrpc, 1, 1);
        ctx->rrpc = NULL;
    }
}

void ecgsyn_init_default_params(EcgSynParams *p)
{
    if (p == NULL) return;

    p->ecg_fs      = 256;
    p->internal_fs = 256;
    p->n_beats     = 16;

    p->noise_mv    = 0.0f;
    p->hr_mean     = 60.0f;
    p->hr_std      = 1.0f;

    p->flo         = 0.1f;
    p->fhi         = 0.25f;
    p->flo_std     = 0.01f;
    p->fhi_std     = 0.01f;
    p->lf_hf_ratio = 0.5f;

    p->seed_init   = 1;

    p->xinitial    = 1.0f;
    p->yinitial    = 0.0f;
    p->zinitial    = 0.04f;
}

void ecgsyn_init_context(EcgSynContext *ctx)
{
    if (ctx == NULL) return;

    ctx->h     = 0.0f;
    ctx->fhi   = 0.0f;
    ctx->rseed = 0;

    ctx->rr    = NULL;
    ctx->rrpc  = NULL;

    for (int i = 0; i < 6; i++) {
        ctx->ti[i] = 0.0f;
        ctx->ai[i] = 0.0f;
        ctx->bi[i] = 0.0f;
    }
}

// ------------------------------------------------------------
// Helpers for implicit Tustin ECGSYN
// ------------------------------------------------------------
static inline float wrap_pm_pi(float x)
{
    return atan2f(sinf(x), cosf(x));
}

static inline void ecgsyn_rhs(float t0, const float x[], float f[])
{
    const float xx = x[1];
    const float yy = x[2];
    const float zz = x[3];

    const float r  = sqrtf(xx * xx + yy * yy);
    const float a0 = 1.0f - r;

    const float w0 = angfreq(t0);
    const float theta = atan2f(yy, xx);
    const float zbase = 0.005f * sinf(2.0f * PI_D * g_ctx->fhi * t0);

    float zsum = 0.0f;
    for (int k = 1; k <= 5; k++) {
        const float dth = wrap_pm_pi(theta - g_ctx->ti[k]);
        const float b   = g_ctx->bi[k];
        const float e   = expf(-(dth * dth) / (2.0f * b * b));
        zsum += g_ctx->ai[k] * dth * e;
    }

    f[1] = a0 * xx - w0 * yy;
    f[2] = a0 * yy + w0 * xx;
    f[3] = -zsum - (zz - zbase);
}

static inline void ecgsyn_jacobian(float t0, const float x[], float J[4][4])
{
    (void)t0;

    const float xx = x[1];
    const float yy = x[2];

    float r2 = xx * xx + yy * yy;
    float r  = sqrtf(r2);

    if (r < 1e-12f) {
        r  = 1e-12f;
        r2 = r * r;
    }

    const float w0 = angfreq(t0);

    // x,y rows
    const float df1_dx = 1.0f - r - (xx * xx) / r;
    const float df1_dy = -w0 - (xx * yy) / r;
    const float df2_dx =  w0 - (xx * yy) / r;
    const float df2_dy = 1.0f - r - (yy * yy) / r;

    // z row
    const float theta = atan2f(yy, xx);
    const float dtheta_dx = -yy / r2;
    const float dtheta_dy =  xx / r2;

    float sum_dx = 0.0f;
    float sum_dy = 0.0f;

    for (int k = 1; k <= 5; k++) {
        const float dth = wrap_pm_pi(theta - g_ctx->ti[k]);
        const float b   = g_ctx->bi[k];
        const float b2  = b * b;
        const float e   = expf(-(dth * dth) / (2.0f * b2));
        const float gp  = e * (1.0f - (dth * dth) / b2);

        sum_dx += g_ctx->ai[k] * gp * dtheta_dx;
        sum_dy += g_ctx->ai[k] * gp * dtheta_dy;
    }

    // zero all
    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            J[i][j] = 0.0f;
        }
    }

    J[1][1] = df1_dx;
    J[1][2] = df1_dy;
    J[1][3] = 0.0f;

    J[2][1] = df2_dx;
    J[2][2] = df2_dy;
    J[2][3] = 0.0f;

    J[3][1] = -sum_dx;
    J[3][2] = -sum_dy;
    J[3][3] = -1.0f;
}

static inline bool solve3x3(float A[4][4], float b[4], float x[4])
{
    // Gaussian elimination with partial pivoting, 1-based indexing
    float M[4][5];

    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) M[i][j] = A[i][j];
        M[i][4] = b[i];
    }

    for (int k = 1; k <= 3; k++) {
        int piv = k;
        float amax = fabsf(M[k][k]);
        for (int i = k + 1; i <= 3; i++) {
            float v = fabsf(M[i][k]);
            if (v > amax) {
                amax = v;
                piv = i;
            }
        }

        if (amax < 1e-12f) return false;

        if (piv != k) {
            for (int j = k; j <= 4; j++) {
                float tmp = M[k][j];
                M[k][j] = M[piv][j];
                M[piv][j] = tmp;
            }
        }

        for (int i = k + 1; i <= 3; i++) {
            float m = M[i][k] / M[k][k];
            for (int j = k; j <= 4; j++) {
                M[i][j] -= m * M[k][j];
            }
        }
    }

    for (int i = 3; i >= 1; i--) {
        float s = M[i][4];
        for (int j = i + 1; j <= 3; j++) s -= M[i][j] * x[j];
        x[i] = s / M[i][i];
    }

    return true;
}

static inline bool implicit_tustin_step(float y[], float t0, float h, float yout[])
{
    const int max_iter = 6;
    const float tol = 1e-6f;

    float fp[4], fn[4];
    float xn[4];
    float R[4];
    float Jf[4][4];
    float J[4][4];
    float dx[4];
    float minusR[4];

    ecgsyn_rhs(t0, y, fp);

    // initial guess: previous state
    for (int i = 1; i <= 3; i++) {
        xn[i] = y[i];
    }

    for (int iter = 0; iter < max_iter; iter++) {
        ecgsyn_rhs(t0 + h, xn, fn);
        ecgsyn_jacobian(t0 + h, xn, Jf);

        for (int i = 1; i <= 3; i++) {
            R[i] = xn[i] - y[i] - 0.5f * h * (fn[i] + fp[i]);
        }

        for (int i = 1; i <= 3; i++) {
            for (int j = 1; j <= 3; j++) {
                J[i][j] = -0.5f * h * Jf[i][j];
            }
            J[i][i] += 1.0f;
        }

        for (int i = 1; i <= 3; i++) {
            minusR[i] = -R[i];
            dx[i] = 0.0f;
        }

        if (!solve3x3(J, minusR, dx)) {
            return false;
        }

        float err = 0.0f;
        for (int i = 1; i <= 3; i++) {
            xn[i] += dx[i];
            float adx = fabsf(dx[i]);
            if (adx > err) err = adx;
        }

        if (err < tol) {
            for (int i = 1; i <= 3; i++) yout[i] = xn[i];
            return true;
        }
    }

    // still return the last Newton iterate
    for (int i = 1; i <= 3; i++) yout[i] = xn[i];
    return true;
}
