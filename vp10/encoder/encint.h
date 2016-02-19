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

# include "vp10/common/state.h"
# include "vpx_dsp/entenc.h"
# include "vp10/common/daala.h"
//# include "block_size_enc.h"

/*Constants for the packet state machine specific to the encoder.*/
/*No packet currently ready to output.*/
# define OD_PACKET_EMPTY       (0)
/*A packet ready to output.*/
# define OD_PACKET_READY       (1)
/*The number of fractional bits of precision in our \lambda values.*/
# define OD_LAMBDA_SCALE       (2)
/*The number of bits of precision to add to distortion values to match
   \lambda*R.*/
# define OD_ERROR_SCALE        (OD_LAMBDA_SCALE + OD_BITRES)


/*Unsanitized user parameters*/
struct od_params_ctx {
  /*Set using OD_SET_MV_LEVEL_MIN*/
  int mv_level_min;
  /*Set using OD_SET_MV_LEVEL_MAX*/
  int mv_level_max;
};

struct daala_enc_ctx{
  od_state state;
  od_ec_enc ec;
  int use_activity_masking;
  int qm;
};

/** Holds important encoder information so we can roll back decisions */
struct od_rollback_buffer {
  od_ec_enc ec;
  od_adapt_ctx adapt;
};

void od_encode_checkpoint(const daala_enc_ctx *enc, od_rollback_buffer *rbuf);
void od_encode_rollback(daala_enc_ctx *enc, const od_rollback_buffer *rbuf);

#endif
