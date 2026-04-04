#ifndef ECGSYN_MODEL_H
#define ECGSYN_MODEL_H

#include <Arduino.h>
#include <math.h>
#include <stdlib.h>

#ifndef PI_D
#define PI_D 3.1415926535897932384626433832795
#endif

namespace ecgsyn {

// ============================================================
// Parameters: keep these aligned with ecgsyn.c defaults
// ============================================================
struct EcgSynParams {
  int    ecg_fs        = 256;   // sfecg
  int    internal_fs   = 256;   // sf
  int    n_beats       = 16;    // N (approximate)
  double noise_mv      = 0.0;   // Anoise
  double hr_mean       = 60.0;  // hrmean
  double hr_std        = 1.0;   // hrstd
  double flo           = 0.1;
  double fhi           = 0.25;
  double flo_std       = 0.01;
  double fhi_std       = 0.01;
  double lf_hf_ratio   = 0.5;
  int    seed_init     = 1;

  double xinitial      = 1.0;
  double yinitial      = 0.0;
  double zinitial      = 0.04;
};

struct EcgSynContext {
  double h = 1.0 / 256.0;
  double fhi = 0.25;
  long   rseed = -1;

  // 1-based vectors, like ecgsyn.c
  double* rr   = nullptr;
  double* rrpc = nullptr;
  double* ti   = nullptr;
  double* ai   = nullptr;
  double* bi   = nullptr;
};

static EcgSynContext* g_ctx = nullptr;

// ============================================================
// 1-based vector helpers
// ============================================================
static inline double* mallocVect(long n0, long nx) {
  const long OFFSET = 1;
  double* vect = (double*)malloc((size_t)((nx - n0 + 1 + OFFSET) * sizeof(double)));
  return vect ? (vect - n0 + OFFSET) : nullptr;
}

static inline void freeVect(double* vect, long n0, long /*nx*/) {
  const long OFFSET = 1;
  if (vect) free((void*)(vect + n0 - OFFSET));
}

static inline double stdev(double* x, int n) {
  double add = 0.0;
  for (int j = 1; j <= n; j++) add += x[j];
  double mean = add / n;

  double total = 0.0;
  for (int j = 1; j <= n; j++) {
    double diff = x[j] - mean;
    total += diff * diff;
  }
  return sqrt(total / (n - 1));
}

// ============================================================
// Numerical Recipes ran1
// ============================================================
static inline float ran1(long* idum) {
  const long IA   = 16807;
  const long IM   = 2147483647;
  const double AM = (1.0 / IM);
  const long IQ   = 127773;
  const long IR   = 2836;
  const int NTAB  = 32;
  const long NDIV = (1 + (IM - 1) / NTAB);
  const double EPS  = 1.2e-7;
  const double RNMX = (1.0 - EPS);

  int j;
  long k;
  static long iy = 0;
  static long iv[NTAB];
  double temp;

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
  j = iy / NDIV;
  iy = iv[j];
  iv[j] = *idum;
  temp = AM * iy;
  return (float)((temp > RNMX) ? RNMX : temp);
}

// ============================================================
// Numerical Recipes dfour1 (1-based interleaved complex array)
// data[1],data[2] = real,imag of sample 1, etc.
// ============================================================
static inline void dfour1(double data[], unsigned long nn, int isign) {
  unsigned long n, mmax, m, j, istep, i;
  double wtemp, wr, wpr, wpi, wi, theta;
  double tempr, tempi;

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
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
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

// ============================================================
// RR process exactly in ecgsyn.c style
// ============================================================
static inline bool rrprocess(
    double* rr,
    double flo, double fhi,
    double flostd, double fhistd, double lfhfratio,
    double hrmean, double hrstd,
    double sf, int n,
    long* rseed)
{
  int i;
  double c1, c2, w1, w2, sig1, sig2, rrmean, rrstd, xstd, ratio;
  double df, dw1, dw2;
  double *w, *Hw, *Sw, *ph0, *ph, *SwC;

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
    Hw[i] = sig1 * exp(-dw1 * dw1 / (2.0 * c1 * c1)) / sqrt(2.0 * PI_D * c1 * c1)
          + sig2 * exp(-dw2 * dw2 / (2.0 * c2 * c2)) / sqrt(2.0 * PI_D * c2 * c2);
  }

  for (i = 1; i <= n / 2; i++)     Sw[i] = (sf / 2.0) * sqrt(Hw[i]);
  for (i = n / 2 + 1; i <= n; i++) Sw[i] = (sf / 2.0) * sqrt(Hw[n - i + 1]);

  for (i = 1; i <= n / 2 - 1; i++) ph0[i] = 2.0 * PI_D * ran1(rseed);
  ph[1] = 0.0;
  for (i = 1; i <= n / 2 - 1; i++) ph[i + 1] = ph0[i];
  ph[n / 2 + 1] = 0.0;
  for (i = 1; i <= n / 2 - 1; i++) ph[n - i + 1] = -ph0[i];

  for (i = 1; i <= n; i++) SwC[2 * i - 1] = Sw[i] * cos(ph[i]);
  for (i = 1; i <= n; i++) SwC[2 * i]     = Sw[i] * sin(ph[i]);

  dfour1(SwC, n, -1);

  for (i = 1; i <= n; i++) rr[i] = (1.0 / n) * SwC[2 * i - 1];

  xstd = stdev(rr, n);
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

// ============================================================
// ODE pieces
// ============================================================
static inline double angfreq(double t) {
  int i = 1 + (int)floor(t / g_ctx->h);
  return 2.0 * PI_D / g_ctx->rrpc[i];
}

static inline void derivspqrst(double t0, double x[], double dxdt[]) {
  int i, k;
  double a0, w0, r0, x0, y0;
  double t, dt, dt2, zbase;

  k = 5;
  w0 = angfreq(t0);
  r0 = 1.0; x0 = 0.0; y0 = 0.0;
  a0 = 1.0 - sqrt((x[1] - x0) * (x[1] - x0) + (x[2] - y0) * (x[2] - y0)) / r0;

  zbase = 0.005 * sin(2.0 * PI_D * g_ctx->fhi * t0);

  t = atan2(x[2], x[1]);
  dxdt[1] = a0 * (x[1] - x0) - w0 * (x[2] - y0);
  dxdt[2] = a0 * (x[2] - y0) + w0 * (x[1] - x0);
  dxdt[3] = 0.0;

  for (i = 1; i <= k; i++) {
    dt = fmod(t - g_ctx->ti[i], 2.0 * PI_D);
    dt2 = dt * dt;
    dxdt[3] += -g_ctx->ai[i] * dt * exp(-0.5 * dt2 / (g_ctx->bi[i] * g_ctx->bi[i]));
  }
  dxdt[3] += -1.0 * (x[3] - zbase);
}

static inline void drk4(double y[], int n, double x, double h, double yout[]) {
  int i;
  double xh, hh, h6, *dydx, *dym, *dyt, *yt;

  dydx = mallocVect(1, n);
  dym  = mallocVect(1, n);
  dyt  = mallocVect(1, n);
  yt   = mallocVect(1, n);

  if (!dydx || !dym || !dyt || !yt) {
    freeVect(dydx, 1, n);
    freeVect(dym, 1, n);
    freeVect(dyt, 1, n);
    freeVect(yt, 1, n);
    return;
  }

  hh = h * 0.5;
  h6 = h / 6.0;
  xh = x + hh;

  derivspqrst(x, y, dydx);
  for (i = 1; i <= n; i++) yt[i] = y[i] + hh * dydx[i];

  derivspqrst(xh, yt, dyt);
  for (i = 1; i <= n; i++) yt[i] = y[i] + hh * dyt[i];

  derivspqrst(xh, yt, dym);
  for (i = 1; i <= n; i++) {
    yt[i] = y[i] + h * dym[i];
    dym[i] += dyt[i];
  }

  derivspqrst(x + h, yt, dyt);
  for (i = 1; i <= n; i++) {
    yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
  }

  freeVect(dydx, 1, n);
  freeVect(dym, 1, n);
  freeVect(dyt, 1, n);
  freeVect(yt, 1, n);
}

// ============================================================
// Build ECG block, faithful to ecgsyn.c with one safety repair:
// rrpc is filled using a safe two-pass cumulative method.
// ============================================================
static inline bool build_block_mv(
    const EcgSynParams& p,
    double*& out_mv,
    int& out_len,
    EcgSynContext& ctx)
{
  const int q = (int)llround((double)p.internal_fs / (double)p.ecg_fs);
  if (q <= 0 || q != (p.internal_fs / p.ecg_fs) || (p.internal_fs % p.ecg_fs) != 0) {
    return false;
  }

  ctx.h = 1.0 / p.internal_fs;
  ctx.fhi = p.fhi;
  ctx.rseed = -p.seed_init;

  ctx.ti = mallocVect(1, 5);
  ctx.ai = mallocVect(1, 5);
  ctx.bi = mallocVect(1, 5);
  if (!ctx.ti || !ctx.ai || !ctx.bi) return false;

  // C defaults
  ctx.ti[1] = -60.0; ctx.ti[2] = -15.0; ctx.ti[3] = 0.0;  ctx.ti[4] = 15.0; ctx.ti[5] = 90.0;
  ctx.ai[1] = 1.2;   ctx.ai[2] = -5.0;  ctx.ai[3] = 30.0; ctx.ai[4] = -7.5; ctx.ai[5] = 0.75;
  ctx.bi[1] = 0.25;  ctx.bi[2] = 0.1;   ctx.bi[3] = 0.1;  ctx.bi[4] = 0.1;  ctx.bi[5] = 0.4;

  for (int i = 1; i <= 5; i++) ctx.ti[i] *= PI_D / 180.0;

  // C heart-rate scaling
  double hrfact = sqrt(p.hr_mean / 60.0);
  double hrfact2 = sqrt(hrfact);
  for (int i = 1; i <= 5; i++) ctx.bi[i] *= hrfact;
  ctx.ti[1] *= hrfact2;
  ctx.ti[2] *= hrfact;
  ctx.ti[3] *= 1.0;
  ctx.ti[4] *= hrfact;
  ctx.ti[5] *= 1.0;

  double rrmean = 60.0 / p.hr_mean;
  int Nrr = 1;
  while (Nrr < (int)(p.n_beats * rrmean * p.internal_fs)) Nrr <<= 1;

  ctx.rr = mallocVect(1, Nrr);
  if (!ctx.rr) return false;

  if (!rrprocess(ctx.rr, p.flo, p.fhi, p.flo_std, p.fhi_std,
                 p.lf_hf_ratio, p.hr_mean, p.hr_std,
                 (double)p.internal_fs, Nrr, &ctx.rseed)) {
    return false;
  }

  // Safe replacement for the original rrpc expansion
  double tecg = 0.0;
  for (int beat = 1; beat <= Nrr; beat++) tecg += ctx.rr[beat];
  int Nt = (int)llround(tecg / ctx.h);
  if (Nt <= 0) return false;

  ctx.rrpc = mallocVect(1, Nt);
  if (!ctx.rrpc) return false;

  tecg = 0.0;
  int sample_idx = 1;
  for (int beat = 1; beat <= Nrr; beat++) {
    tecg += ctx.rr[beat];
    int sample_end = (int)llround(tecg / ctx.h);
    if (sample_end > Nt) sample_end = Nt;
    for (int k = sample_idx; k <= sample_end; k++) ctx.rrpc[k] = ctx.rr[beat];
    sample_idx = sample_end + 1;
  }

  double* x   = mallocVect(1, 3);
  double* xn  = mallocVect(1, 3);
  double* zt  = mallocVect(1, Nt);
  if (!x || !xn || !zt) {
    freeVect(x, 1, 3);
    freeVect(xn, 1, 3);
    freeVect(zt, 1, Nt);
    return false;
  }

  x[1] = p.xinitial;
  x[2] = p.yinitial;
  x[3] = p.zinitial;

  g_ctx = &ctx;
  double timev = 0.0;
  for (int i = 1; i <= Nt; i++) {
    zt[i] = x[3];
    drk4(x, 3, timev, ctx.h, xn);
    x[1] = xn[1];
    x[2] = xn[2];
    x[3] = xn[3];
    timev += ctx.h;
  }
  g_ctx = nullptr;

  out_len = (Nt + q - 1) / q;
  if (out_len <= 0) return false;

  out_mv = (double*)malloc((size_t)out_len * sizeof(double));
  if (!out_mv) return false;

  double zmin = zt[1], zmax = zt[1];
  for (int i = 2; i <= Nt; i++) {
    if (zt[i] < zmin) zmin = zt[i];
    else if (zt[i] > zmax) zmax = zt[i];
  }

  double zrange = zmax - zmin;
  if (zrange < 1e-12) zrange = 1.0;

  int j = 0;
  for (int i = 1; i <= Nt; i += q) {
    double z = (zt[i] - zmin) * (1.6) / zrange - 0.4;
    z += p.noise_mv * (2.0 * ran1(&ctx.rseed) - 1.0);
    out_mv[j++] = z;
  }

  freeVect(x, 1, 3);
  freeVect(xn, 1, 3);
  freeVect(zt, 1, Nt);
  return true;
}

static inline void free_context(EcgSynContext& ctx) {
  if (ctx.rr)   { freeVect(ctx.rr,   1, 1); ctx.rr   = nullptr; }
  if (ctx.rrpc) { freeVect(ctx.rrpc, 1, 1); ctx.rrpc = nullptr; }
  if (ctx.ti)   { freeVect(ctx.ti,   1, 5); ctx.ti   = nullptr; }
  if (ctx.ai)   { freeVect(ctx.ai,   1, 5); ctx.ai   = nullptr; }
  if (ctx.bi)   { freeVect(ctx.bi,   1, 5); ctx.bi   = nullptr; }
}

} // namespace ecgsyn_faithful

#endif