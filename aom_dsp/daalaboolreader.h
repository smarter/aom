/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_DSP_DAALABOOLREADER_H_
#define AOM_DSP_DAALABOOLREADER_H_

#include "aom_dsp/entdec.h"
#include "aom_dsp/prob.h"
#if CONFIG_ACCOUNTING
#include "av1/common/accounting.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct daala_reader {
  const uint8_t *buffer;
  const uint8_t *buffer_end;
  od_ec_dec ec;
};

typedef struct daala_reader daala_reader;

#if CONFIG_ACCOUNTING
# define aom_daala_read(r, prob, str) aom_daala_read_(r, prob, str)
# define aom_daala_read_bit(r, str) aom_daala_read_bit_(r, str)
# define daala_read_tree_cdf(r, tree, probs, str) daala_read_tree_cdf_(r, cdf, nsymbs, str)
# define daala_read_tree_bits(r, tree, probs, str) daala_read_tree_bits_(r, tree, probs, str)
#else
# define aom_daala_read(r, prob, str) aom_daala_read_(r, prob)
# define aom_daala_read_bit(r, str) aom_daala_read_bit_(r)
# define daala_read_tree_cdf(r, tree, probs, str) daala_read_tree_cdf_(r, cdf, nsymbs)
# define daala_read_tree_bits(r, tree, probs, str) daala_read_tree_bits_(r, tree, probs)
#endif

int aom_daala_reader_init(daala_reader *r, const uint8_t *buffer, int size);
const uint8_t *aom_daala_reader_find_end(daala_reader *r);
ptrdiff_t aom_daala_reader_tell(const daala_reader *r);

static INLINE int aom_daala_read_(daala_reader *r, int prob AOM_ACCT_STR_PARAM) {
  if (prob == 128) {
    return od_ec_dec_bits(&r->ec, 1, acct_str);
  } else {
    int p = ((prob << 15) + (256 - prob)) >> 8;
    return od_ec_decode_bool_q15(&r->ec, p, acct_str);
  }
}

static INLINE int aom_daala_read_bit_(daala_reader *r AOM_ACCT_STR_PARAM) {
  return aom_daala_read(r, 128, acct_str);
}

static INLINE int aom_daala_reader_has_error(daala_reader *r) {
  return r->ec.error;
}

static INLINE int daala_read_tree_bits_(daala_reader *r,
                                       const aom_tree_index *tree,
                                       const aom_prob *probs AOM_ACCT_STR_PARAM) {
  aom_tree_index i = 0;
  do {
    uint16_t cdf[16];
    aom_tree_index index[16];
    int path[16];
    int dist[16];
    int nsymbs;
    int symb;
    nsymbs = tree_to_cdf(tree, probs, i, cdf, index, path, dist);
    symb = od_ec_decode_cdf_q15(&r->ec, cdf, nsymbs, acct_str);
    OD_ASSERT(symb >= 0 && symb < nsymbs);
    i = index[symb];
  } while (i > 0);
  return -i;
}

static INLINE int daala_read_tree_cdf_(daala_reader *r, const uint16_t *cdf,
                                      int nsymbs AOM_ACCT_STR_PARAM) {
  return od_ec_decode_cdf_q15(&r->ec, cdf, nsymbs, acct_str);
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif
