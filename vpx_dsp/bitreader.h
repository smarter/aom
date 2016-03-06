/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VPX_DSP_BITREADER_H_
#define VPX_DSP_BITREADER_H_

#include <stddef.h>
#include <limits.h>

#include "./vpx_config.h"
#include "vpx_ports/mem.h"
#include "vpx/vp8dx.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/prob.h"
#if CONFIG_DAALA_EC
#include "vpx_dsp/entdec.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef size_t BD_VALUE;

#define BD_VALUE_SIZE ((int)sizeof(BD_VALUE) * CHAR_BIT)

// This is meant to be a large, positive constant that can still be efficiently
// loaded as an immediate (on platforms like ARM, for example).
// Even relatively modest values like 100 would work fine.
#define LOTS_OF_BITS 0x40000000

typedef struct {
  // Be careful when reordering this struct, it may impact the cache negatively.
  BD_VALUE value;
  unsigned int range;
  int count;
  const uint8_t *buffer_end;
  const uint8_t *buffer;
  vpx_decrypt_cb decrypt_cb;
  void *decrypt_state;
  uint8_t clear_buffer[sizeof(BD_VALUE) + 1];
#if CONFIG_DAALA_EC
  od_ec_dec ec;
#endif
} vpx_reader;

int vpx_reader_init(vpx_reader *r, const uint8_t *buffer, size_t size,
                    vpx_decrypt_cb decrypt_cb, void *decrypt_state);

void vpx_reader_fill(vpx_reader *r);

const uint8_t *vpx_reader_find_end(vpx_reader *r);

static INLINE int vpx_reader_has_error(vpx_reader *r) {
  // Check if we have reached the end of the buffer.
  //
  // Variable 'count' stores the number of bits in the 'value' buffer, minus
  // 8. The top byte is part of the algorithm, and the remainder is buffered
  // to be shifted into it. So if count == 8, the top 16 bits of 'value' are
  // occupied, 8 for the algorithm and 8 in the buffer.
  //
  // When reading a byte from the user's buffer, count is filled with 8 and
  // one byte is filled into the value buffer. When we reach the end of the
  // data, count is additionally filled with LOTS_OF_BITS. So when
  // count == LOTS_OF_BITS - 1, the user's data has been exhausted.
  //
  // 1 if we have tried to decode bits after the end of stream was encountered.
  // 0 No error.
  return r->count > BD_VALUE_SIZE && r->count < LOTS_OF_BITS;
}

static INLINE int vpx_read(vpx_reader *r, int prob) {
#if CONFIG_DAALA_EC
  if (prob == 128) {
    return od_ec_dec_bits(&r->ec, 1, "vpx_bits");
  }
  else {
    int p = (32768*prob + (256 - prob)) >> 8;
    return od_ec_decode_bool_q15(&r->ec, p, "vpx");
  }
#else
  unsigned int bit = 0;
  BD_VALUE value;
  BD_VALUE bigsplit;
  int count;
  unsigned int range;
  unsigned int split = (r->range * prob + (256 - prob)) >> CHAR_BIT;

  if (r->count < 0) vpx_reader_fill(r);

  value = r->value;
  count = r->count;

  bigsplit = (BD_VALUE)split << (BD_VALUE_SIZE - CHAR_BIT);

  range = split;

  if (value >= bigsplit) {
    range = r->range - split;
    value = value - bigsplit;
    bit = 1;
  }

  {
    register int shift = vpx_norm[range];
    range <<= shift;
    value <<= shift;
    count -= shift;
  }
  r->value = value;
  r->count = count;
  r->range = range;

  return bit;
#endif
}

static INLINE int vpx_read_bit(vpx_reader *r) {
  return vpx_read(r, 128);  // vpx_prob_half
}

static INLINE int vpx_read_literal(vpx_reader *r, int bits) {
  int literal = 0, bit;

  for (bit = bits - 1; bit >= 0; bit--) literal |= vpx_read_bit(r) << bit;

  return literal;
}

static INLINE int vpx_read_tree(vpx_reader *r, const vpx_tree_index *tree,
                                const vpx_prob *probs) {
  vpx_tree_index i = 0;
#if CONFIG_DAALA_EC
  do {
    uint16_t cdf[16];
    vpx_tree_index index[16];
    int path[16];
    int dist[16];
    int nsymbs;
    int symb;
    nsymbs = tree_to_cdf(tree, probs, i, cdf, index, path, dist);
    symb = od_ec_decode_cdf_q15(&r->ec, cdf, nsymbs, "vpx");
    OD_ASSERT(symb >= 0 && symb < nsymbs);
    i = index[symb];
  }
  while (i > 0);
#else
  while ((i = tree[i + vpx_read(r, probs[i >> 1])]) > 0) continue;
#endif
  return -i;
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VPX_DSP_BITREADER_H_
