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
# include "daala.h"
# include "pvq.h"
# include "generic_code.h"

extern const od_coeff OD_DC_QM[OD_NBSIZES - 1][2];

extern const int OD_HAAR_QM[2][OD_LOG_BSIZE_MAX];

/*Adaptation speed of scalar Laplace encoding.*/
# define OD_SCALAR_ADAPT_SPEED (4)


/*This should be a power of 2, and at least 8.*/
# define OD_UMV_CLAMP (32)
/*Half the value of OD_SUBPEL_FILTER_TAP_SIZE in mc.h.*/
# define OD_RESAMPLE_PADDING (3)
/*Buffer padding alignment dictated by SIMD.*/
# define OD_PADDING_ALIGN (32)
/*Actual size of padding on all sides of reference buffers.*/
# define OD_BUFFER_PADDING \
 ((OD_UMV_CLAMP + OD_RESAMPLE_PADDING + OD_PADDING_ALIGN - 1) \
 /OD_PADDING_ALIGN*OD_PADDING_ALIGN)


struct od_adapt_ctx {
  /* Support for PVQ encode/decode */
  od_pvq_adapt_ctx pvq;

  generic_encoder model_dc[OD_NPLANES_MAX];

  int ex_sb_dc[OD_NPLANES_MAX];
  int ex_dc[OD_NPLANES_MAX][OD_NBSIZES][3];
  int ex_g[OD_NPLANES_MAX][OD_NBSIZES];

  /* Joint skip flag for DC and AC */
  uint16_t skip_cdf[OD_NBSIZES*2][5];
  int skip_increment;
  /* 4 possible values for the sblock above (or skip), 4 for the sblock to the
     left (or skip). */
  uint16_t q_cdf[4*4][4];
  int q_increment;
};

struct od_state{
  od_adapt_ctx        adapt;
  int32_t         frame_width;
  int32_t         frame_height;

  /** Number of horizontal motion-vector blocks. */
  int                 nhmvbs;
  /** Number of vertical motion-vector blocks. */
  int                 nvmvbs;
  int                 nhsb;
  int                 nvsb;
  /** Each 8x8 block of pixels in the image (+ one superblock of
      padding on each side) has a corresponding byte in this array, and
      every 64x64 superblock is represented by 64 (8 by 8) entries
      ((8 * 8) * (8 * 8) == 64 * 64) that encode the block size decisions
      for the superblock. The entry format is:
      - 0 means the 8x8 block has been split into 4x4 blocks
      - 1 means the 8x8 block is an 8x8 block
      - 2 means the 8x8 block is part of a 16x16 block
      - 3 means the 8x8 block is part of a 32x32 block.
      - 4 means the 8x8 block is part of a 64x64 block.
      The padding is filled as though it consisted of 64x64 blocks.

      E.g., `state->bsize[j * state->bstride + i]` accesses the i'th 8x8
      block in the j'th row of 8x8 blocks.

      The `bstride` member has the distance between vertically adjacent
      entries (horizontally adjacent entries are adjacent in memory). */
  unsigned char *bsize;
  int                 bstride;
  unsigned char *bskip[3];
  int                 skip_stride;
  od_coeff           *(sb_dc_mem[OD_NPLANES_MAX]);
  int quantizer[OD_NPLANES_MAX];
  int coded_quantizer[OD_NPLANES_MAX];
  unsigned char pvq_qm_q4[OD_NPLANES_MAX][OD_QM_SIZE];
  /*This provides context for the quantizer CDF.*/
  unsigned char *sb_q_scaling;
  /*Magnitude compensated quantization matrices and their inverses.
   1 per block-size and decimation factor (i.e. OD_NBSIZES*2*(OD_BSIZE_MAX^2)),
   assuming 2 possible decimation values (see OD_BASIS_MAG).*/
  int16_t *qm;
  int16_t *qm_inv;
};

void od_adapt_ctx_reset(od_adapt_ctx *state, int is_keyframe);
void od_init_skipped_coeffs(od_coeff *d, od_coeff *pred, int is_keyframe,
 int bo, int n, int w);


#endif
