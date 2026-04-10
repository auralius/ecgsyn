#ifndef RD_MODEL_H
#define RD_MODEL_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
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

    int newton_iters;   /* used by Tustin */
} RDParams;

typedef struct {
    float dt;
    float Gamma_t;
    float x[4];
} RDContext;

void rd_init_default_params(RDParams *p);
bool rd_init(const RDParams *p, RDContext *ctx);
bool build_rd_block_mv(const RDParams *p, float **out_mv, int *out_len, RDContext *ctx);

#ifdef __cplusplus
}
#endif

#endif
