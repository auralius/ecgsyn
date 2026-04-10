#include "rd_model.h"
#include <math.h>
#include <stdlib.h>

/* ------------------------------------------------------------
 * Select one solver
 * ------------------------------------------------------------ */
#define RD_SOLVER_TUSTIN   1
/* #define RD_SOLVER_RK4   1 */

#if defined(RD_SOLVER_TUSTIN) && defined(RD_SOLVER_RK4)
#error "Select only one RD solver"
#endif

#if !defined(RD_SOLVER_TUSTIN) && !defined(RD_SOLVER_RK4)
#error "Select one RD solver"
#endif

/* ------------------------------------------------------------
 * Core model
 * ------------------------------------------------------------ */
static inline void rd_f_raw(const RDParams *p, const float x[4], float f[4])
{
    const float x1 = x[0];
    const float x2 = x[1];
    const float x3 = x[2];
    const float x4 = x[3];

    f[0] = x1 - x2 - p->C * x1 * x2 - x1 * x2 * x2;
    f[1] = p->H * x1 - 3.0f * x2 + p->C * x1 * x2 + x1 * x2 * x2 + p->beta * (x4 - x2);
    f[2] = x3 - x4 - p->C * x3 * x4 - x3 * x4 * x4;
    f[3] = p->H * x3 - 3.0f * x4 + p->C * x3 * x4 + x3 * x4 * x4 + 2.0f * p->beta * (x2 - x4);
}

static inline void rd_f_scaled(const RDParams *p, const RDContext *ctx, const float x[4], float f[4])
{
    rd_f_raw(p, x, f);
    for (int i = 0; i < 4; i++) {
        f[i] *= ctx->Gamma_t;
    }
}

static inline float rd_ecg_mv(const RDParams *p, const float x[4])
{
    return p->alpha1 * x[0]
         + p->alpha2 * x[1]
         + p->alpha3 * x[2]
         + p->alpha4 * x[3];
}

/* ------------------------------------------------------------
 * Tustin helpers
 * ------------------------------------------------------------ */
#if defined(RD_SOLVER_TUSTIN)

static bool solve4x4(const float Ain[4][4], const float bin[4], float xout[4])
{
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

static inline void rd_df_raw(const RDParams *p, const float x[4], float J[4][4])
{
    const float x1 = x[0];
    const float x2 = x[1];
    const float x3 = x[2];
    const float x4 = x[3];

    J[0][0] =  1.0f - p->C * x2 - x2 * x2;
    J[0][1] = -1.0f - p->C * x1 - 2.0f * x1 * x2;
    J[0][2] =  0.0f;
    J[0][3] =  0.0f;

    J[1][0] =  p->H + p->C * x2 + x2 * x2;
    J[1][1] = -3.0f + p->C * x1 + 2.0f * x1 * x2 - p->beta;
    J[1][2] =  0.0f;
    J[1][3] =  p->beta;

    J[2][0] =  0.0f;
    J[2][1] =  0.0f;
    J[2][2] =  1.0f - p->C * x4 - x4 * x4;
    J[2][3] = -1.0f - p->C * x3 - 2.0f * x3 * x4;

    J[3][0] =  0.0f;
    J[3][1] =  2.0f * p->beta;
    J[3][2] =  p->H + p->C * x4 + x4 * x4;
    J[3][3] = -3.0f + p->C * x3 + 2.0f * x3 * x4 - 2.0f * p->beta;
}

static bool rd_tustin_step(const RDParams *p, const RDContext *ctx,
                           const float x_prev[4], float x_next[4])
{
    float xn[4] = { x_prev[0], x_prev[1], x_prev[2], x_prev[3] };
    float fp[4];
    rd_f_scaled(p, ctx, x_prev, fp);

    for (int iter = 0; iter < p->newton_iters; iter++) {
        float fn[4];
        float R[4];
        float Df[4][4];
        float J[4][4];
        float dx[4];

        rd_f_scaled(p, ctx, xn, fn);

        for (int i = 0; i < 4; i++) {
            R[i] = xn[i] - x_prev[i] - 0.5f * ctx->dt * (fn[i] + fp[i]);
        }

        rd_df_raw(p, xn, Df);

        for (int r = 0; r < 4; r++) {
            for (int c = 0; c < 4; c++) {
                J[r][c] = -0.5f * ctx->dt * ctx->Gamma_t * Df[r][c];
            }
            J[r][r] += 1.0f;
        }

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

#endif

/* ------------------------------------------------------------
 * RK4 helper
 * ------------------------------------------------------------ */
#if defined(RD_SOLVER_RK4)

static bool rd_rk4_step(const RDParams *p, const RDContext *ctx,
                        const float x_prev[4], float x_next[4])
{
    const float h = ctx->dt;
    float k1[4], k2[4], k3[4], k4[4];
    float xtmp[4];

    rd_f_scaled(p, ctx, x_prev, k1);

    for (int i = 0; i < 4; i++) {
        xtmp[i] = x_prev[i] + 0.5f * h * k1[i];
    }
    rd_f_scaled(p, ctx, xtmp, k2);

    for (int i = 0; i < 4; i++) {
        xtmp[i] = x_prev[i] + 0.5f * h * k2[i];
    }
    rd_f_scaled(p, ctx, xtmp, k3);

    for (int i = 0; i < 4; i++) {
        xtmp[i] = x_prev[i] + h * k3[i];
    }
    rd_f_scaled(p, ctx, xtmp, k4);

    for (int i = 0; i < 4; i++) {
        x_next[i] = x_prev[i] + (h / 6.0f) * (k1[i] + 2.0f * k2[i] + 2.0f * k3[i] + k4[i]);
    }

    return true;
}

#endif

/* ------------------------------------------------------------
 * Public API
 * ------------------------------------------------------------ */
void rd_init_default_params(RDParams *p)
{
    if (p == NULL) return;

    p->ecg_fs = 200.0f;
    p->internal_fs = 200.0f;
    p->duration_s = 30.0f;

    p->H = 3.0f;
    p->HR_bpm = 95.0f;
    p->C = 1.35f;
    p->beta = 4.0f;

    p->alpha1 = -0.024f;
    p->alpha2 =  0.0216f;
    p->alpha3 = -0.0012f;
    p->alpha4 =  0.12f;

    p->x1_0 = 0.0f;
    p->x2_0 = 0.0f;
    p->x3_0 = 0.1f;
    p->x4_0 = 0.0f;

    p->newton_iters = 4;
}

bool rd_init(const RDParams *p, RDContext *ctx)
{
    if (p == NULL || ctx == NULL) return false;
    if (p->internal_fs <= 0.0f || p->ecg_fs <= 0.0f || p->duration_s <= 0.0f) return false;

    /* keep same-rate assumption for now */
    if (fabsf(p->internal_fs - p->ecg_fs) > 1e-6f) return false;

    ctx->dt = 1.0f / p->internal_fs;
    ctx->Gamma_t = 0.08804f * p->HR_bpm - 0.06754f;

    ctx->x[0] = p->x1_0;
    ctx->x[1] = p->x2_0;
    ctx->x[2] = p->x3_0;
    ctx->x[3] = p->x4_0;

    return true;
}

bool build_rd_block_mv(const RDParams *p, float **out_mv, int *out_len, RDContext *ctx)
{
    if (p == NULL || out_mv == NULL || out_len == NULL || ctx == NULL) return false;

    *out_len = (int)lroundf(p->duration_s * p->ecg_fs);
    if (*out_len <= 0) return false;

    *out_mv = (float*)malloc((size_t)(*out_len) * sizeof(float));
    if (*out_mv == NULL) return false;

    for (int n = 0; n < *out_len; n++) {
        float x_next[4];
        bool ok = false;

#if defined(RD_SOLVER_TUSTIN)
        ok = rd_tustin_step(p, ctx, ctx->x, x_next);
#elif defined(RD_SOLVER_RK4)
        ok = rd_rk4_step(p, ctx, ctx->x, x_next);
#endif

        if (!ok) {
            free(*out_mv);
            *out_mv = NULL;
            *out_len = 0;
            return false;
        }

        for (int i = 0; i < 4; i++) {
            ctx->x[i] = x_next[i];
        }

        (*out_mv)[n] = rd_ecg_mv(p, ctx->x);
    }

    return true;
}
