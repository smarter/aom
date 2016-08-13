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

#ifndef AOM_DSP_BITREADER_H_
#define AOM_DSP_BITREADER_H_

#include <stddef.h>
#include <limits.h>

#include "./aom_config.h"
#include "aom/aomdx.h"
#include "aom/aom_integer.h"
#if CONFIG_ANS
#include "aom_dsp/ansreader.h"
#elif CONFIG_DAALA_EC
#include "aom_dsp/daalaboolreader.h"
#else
#include "aom_dsp/dkboolreader.h"
#endif
#include "aom_dsp/prob.h"
#include "av1/common/odintrin.h"
#if CONFIG_ACCOUNTING
#include "av1/common/accounting.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#if CONFIG_ANS
typedef struct AnsDecoder aom_reader;
#elif CONFIG_DAALA_EC
typedef struct daala_reader aom_reader;
#else
typedef struct aom_dk_reader aom_reader;
#endif

#if CONFIG_ACCOUNTING
# define aom_read(r, prob, str) aom_read_(r, prob, str)
# define aom_read_bit(r, str) aom_read_bit_(r, str)
# define aom_read_tree(r, tree, probs, str) aom_read_tree_(r, tree, probs, str)
# define aom_read_literal(r, bits, str) aom_read_literal_(r, bits, str)
# define aom_read_tree_bits(r, tree, probs, str) aom_read_tree_bits_(r, tree, probs, str)
# define aom_read_tree_cdf(r, cdf, nsymbs, str) aom_read_tree_cdf_(r, cdf, nsymbs, str)
#else
# define AOM_ACCT_STR_PARAM
# define aom_read(r, prob, str) aom_read_(r, prob)
# define aom_read_bit(r, str) aom_read_bit_(r)
# define aom_read_tree(r, tree, probs, str) aom_read_tree_(r, tree, probs)
# define aom_read_literal(r, bits, str) aom_read_literal_(r, bits)
# define aom_read_tree_bits(r, tree, probs, str) aom_read_tree_bits_(r, tree, probs)
# define aom_read_tree_cdf(r, cdf, nsymbs, str) aom_read_tree_cdf_(r, cdf, nsymbs)
#endif

static INLINE int aom_reader_init(aom_reader *r, const uint8_t *buffer,
                                  size_t size, aom_decrypt_cb decrypt_cb,
                                  void *decrypt_state) {
#if CONFIG_ANS
  (void)decrypt_cb;
  (void)decrypt_state;
  assert(size <= INT_MAX);
  return ans_read_init(r, buffer, size);
#elif CONFIG_DAALA_EC
  (void)decrypt_cb;
  (void)decrypt_state;
  return aom_daala_reader_init(r, buffer, size);
#else
  return aom_dk_reader_init(r, buffer, size, decrypt_cb, decrypt_state);
#endif
}

static INLINE const uint8_t *aom_reader_find_end(aom_reader *r) {
#if CONFIG_ANS
  (void)r;
  assert(0 && "Use the raw buffer size with ANS");
  return NULL;
#elif CONFIG_DAALA_EC
  return aom_daala_reader_find_end(r);
#else
  return aom_dk_reader_find_end(r);
#endif
}

static INLINE int aom_reader_has_error(aom_reader *r) {
#if CONFIG_ANS
  return ans_reader_has_error(r);
#elif CONFIG_DAALA_EC
  return aom_daala_reader_has_error(r);
#else
  return aom_dk_reader_has_error(r);
#endif
}

static INLINE ptrdiff_t aom_reader_tell(const aom_reader *r) {
#if CONFIG_DAALA_EC
  return aom_daala_reader_tell(r);
#else
  return aom_dk_reader_tell(r);
#endif
}

#if CONFIG_ACCOUNTING && CONFIG_DAALA_EC
static INLINE void od_process_accounting(od_ec_dec *dec, const char *str) {
  if (dec->accounting != NULL) {
    uint32_t tell;
    tell = od_ec_dec_tell_frac(dec);
    OD_ASSERT(tell >= dec->accounting->last_tell);
    aom_accounting_record(dec->accounting, str, tell - dec->accounting->last_tell);
    dec->accounting->last_tell = tell;
  }
}
#endif

#if CONFIG_ACCOUNTING && !CONFIG_DAALA_EC
static INLINE void aom_process_accounting(aom_reader *r, const char *str) {
  if (r->accounting != NULL) {
    uint32_t tell;
    tell = aom_reader_tell(r);
    aom_accounting_record(r->accounting, str, tell - r->accounting->last_tell);
    r->accounting->last_tell = tell;
  }
}
#endif

static INLINE int aom_read_(aom_reader *r, int prob AOM_ACCT_STR_PARAM) {
  int ret;
#if CONFIG_ANS
  ret = uabs_read(r, prob);
#elif CONFIG_DAALA_EC
  ret = aom_daala_read(r, prob, acct_str);
#else
  ret = aom_dk_read(r, prob);
#endif
#if CONFIG_ACCOUNTING && !CONFIG_DAALA_EC
  aom_process_accounting(r, acct_str);
#endif
  return ret;
}

static INLINE int aom_read_bit_(aom_reader *r AOM_ACCT_STR_PARAM) {
#if CONFIG_ANS
  return uabs_read_bit(r);  // Non trivial optimization at half probability
#else
  return aom_read(r, 128, acct_str);  // aom_prob_half
#endif
}

static INLINE int aom_read_literal_(aom_reader *r, int bits AOM_ACCT_STR_PARAM) {
  int literal = 0, bit;

  for (bit = bits - 1; bit >= 0; bit--) literal |= aom_read_bit(r, acct_str) << bit;

  return literal;
}

static INLINE int aom_read_tree_bits_(aom_reader *r, const aom_tree_index *tree,
                                     const aom_prob *probs AOM_ACCT_STR_PARAM) {
  aom_tree_index i = 0;

  while ((i = tree[i + aom_read(r, probs[i >> 1], acct_str)]) > 0) continue;

  return -i;
}

static INLINE int aom_read_tree_(aom_reader *r, const aom_tree_index *tree,
                                const aom_prob *probs AOM_ACCT_STR_PARAM) {
#if CONFIG_DAALA_EC
  return daala_read_tree_bits(r, tree, probs, acct_str);
#else
  return aom_read_tree_bits(r, tree, probs, acct_str);
#endif
}

static INLINE int aom_read_tree_cdf_(aom_reader *r, const uint16_t *cdf,
                                    int nsymbs AOM_ACCT_STR_PARAM) {
#if CONFIG_RANS
  (void)nsymbs;
  return rans_read(r, cdf);
#elif CONFIG_DAALA_EC
  return daala_read_tree_cdf(r, cdf, nsymbs, acct_str);
#else
  (void)r;
  (void)cdf;
  (void)nsymbs;
  assert(0 && "Unsupported bitreader operation");
  return -1;
#endif
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_DSP_BITREADER_H_
