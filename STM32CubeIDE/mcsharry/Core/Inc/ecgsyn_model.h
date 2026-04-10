#ifndef ECGSYN_MODEL_H
#define ECGSYN_MODEL_H

/*
 * ecgsyn_model.h
 *
 * Program description:
 * This module implements the McSharry synthetic ECG model for embedded
 * waveform generation. The model is used to generate ECG samples with
 * configurable heart-rate variability and PQRST morphology.
 *
 * The state equations can be solved using either:
 * - RK4 (4th-order Runge-Kutta)
 * - Implicit Tustin
 *
 * The generated output is returned in millivolts, so it can be used
 * directly for analysis, visualization, or DAC-based signal playback.
 */

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef PI_D
#define PI_D 3.1415926535897932384626433832795
#endif

typedef struct {
    int   ecg_fs;
    int   internal_fs;
    int   n_beats;
    float noise_mv;
    float hr_mean;
    float hr_std;
    float flo;
    float fhi;
    float flo_std;
    float fhi_std;
    float lf_hf_ratio;
    int   seed_init;

    float xinitial;
    float yinitial;
    float zinitial;
} EcgSynParams;

typedef struct {
    float h;
    float fhi;
    long  rseed;

    float *rr;
    float *rrpc;

    float ti[6];
    float ai[6];
    float bi[6];
} EcgSynContext;

void ecgsyn_init_default_params(EcgSynParams *p);
void ecgsyn_init_context(EcgSynContext *ctx);

bool build_block_mv(const EcgSynParams *p, float **out_mv, int *out_len, EcgSynContext *ctx);
void free_context(EcgSynContext *ctx);

// Helpers for implicit Tustin ECGSYN
static inline bool implicit_tustin_step(float y[], float t0, float h, float yout[]);

#ifdef __cplusplus
}
#endif

#endif
