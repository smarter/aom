/*Daala video codec
Copyright (c) 2006-2010 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#if !defined(_state_H)
# define _state_H (1)

typedef struct od_state          od_state;
typedef struct od_adapt_ctx      od_adapt_ctx;

# include <stdio.h>
# include "odintrin.h"
# include "pvq.h"
# include "generic_code.h"

/*Adaptation speed of scalar Laplace encoding.*/
# define OD_SCALAR_ADAPT_SPEED (4)

struct od_adapt_ctx {
  /* Support for PVQ encode/decode */
  od_pvq_adapt_ctx pvq;

  generic_encoder model_dc[OD_NPLANES_MAX];

  int ex_dc[OD_NPLANES_MAX][OD_NBSIZES][3];
  int ex_g[OD_NPLANES_MAX][OD_NBSIZES];

  /* Joint skip flag for DC and AC */
  uint16_t skip_cdf[OD_NBSIZES*2][4];
  int skip_increment;
};

struct od_state{
  od_adapt_ctx        adapt;
  //unsigned char pvq_qm_q4[OD_NPLANES_MAX][OD_QM_SIZE];
  int16_t *qm;
  int16_t *qm_inv;
};

void od_adapt_ctx_reset(od_adapt_ctx *state, int is_keyframe);
void od_init_skipped_coeffs(int16_t *d, int16_t *pred, int is_keyframe,
 int bo, int n, int w);

#endif
