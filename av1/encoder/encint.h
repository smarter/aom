/*Daala video codec
Copyright (c) 2006-2013 Daala project contributors.  All rights reserved.

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

#if !defined(_encint_H)
# define _encint_H (1)

typedef struct daala_enc_ctx od_enc_ctx;
typedef struct od_params_ctx od_params_ctx;
typedef struct od_rollback_buffer od_rollback_buffer;

# include "av1/common/odintrin.h"
# include "av1/common/state.h"
# include "aom_dsp/entenc.h"

struct daala_enc_ctx{
  od_state state;
  od_ec_enc ec;
  int use_activity_masking;
  int qm;
  /*Normalized PVQ lambda for use where we've already performed
     quantization.*/
  double pvq_norm_lambda;
};

// from daalaenc.h
/**The encoder context.*/
typedef struct daala_enc_ctx daala_enc_ctx;

/** Holds important encoder information so we can roll back decisions */
struct od_rollback_buffer {
  od_ec_enc ec;
  od_adapt_ctx adapt;
};

void od_encode_checkpoint(const daala_enc_ctx *enc, od_rollback_buffer *rbuf);
void od_encode_rollback(daala_enc_ctx *enc, const od_rollback_buffer *rbuf);

#endif
