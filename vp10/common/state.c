#include "vp10/common/odintrin.h"
#include "vp10/common/state.h"

void od_adapt_ctx_reset(od_adapt_ctx *state, int is_keyframe) {
  int i;
  int pli;
  od_adapt_pvq_ctx_reset(&state->pvq, is_keyframe);
  state->skip_increment = 128;
  OD_CDFS_INIT(state->skip_cdf, state->skip_increment >> 2);
  for (pli = 0; pli < OD_NPLANES_MAX; pli++) {
    generic_model_init(&state->model_dc[pli]);
    for (i = 0; i < OD_NBSIZES; i++) {
      state->ex_g[pli][i] = 8;
    }
    for (i = 0; i < 4; i++) {
      int j;
      for (j = 0; j < 3; j++) {
        state->ex_dc[pli][i][j] = pli > 0 ? 8 : 32768;
      }
    }
  }
}

void od_init_skipped_coeffs(int16_t *d, int16_t *pred, int is_keyframe,
 int bo, int n, int w) {
  int i;
  int j;
  if (is_keyframe) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        /* skip DC */
        if (i || j) d[bo + i*w + j] = 0;
      }
    }
  } else {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        d[bo + i*w + j] = pred[i*n + j];
      }
    }
  }
}
