/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef VP10_ENCODER_TREEWRITER_H_
#define VP10_ENCODER_TREEWRITER_H_

#include "vpx_dsp/bitwriter.h"
#if CONFIG_DAALA_EC
#include "vpx_dsp/prob.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

void vp10_tree_probs_from_distribution(vpx_tree tree,
                                       unsigned int branch_ct[/* n - 1 */][2],
                                       const unsigned int num_events[/* n */]);

struct vp10_token {
  int value;
  int len;
};

void vp10_tokens_from_tree(struct vp10_token *, const vpx_tree_index *);

#if CONFIG_DAALA_EC
static INLINE void daala_write_tree(vpx_writer *w, const vpx_tree_index *tree,
                                    const vpx_prob *probs, int bits, int len,
                                    vpx_tree_index i) {
  vpx_tree_index root;
  root = i;
  do {
    uint16_t cdf[16];
    vpx_tree_index index[16];
    int path[16];
    int dist[16];
    int nsymbs;
    int symb;
    int j;
    /* Compute the CDF of the binary tree using the given probabilities. */
    nsymbs = tree_to_cdf(tree, probs, root, cdf, index, path, dist);
    /* Find the symbol to code. */
    symb = -1;
    for (j = 0; j < nsymbs; j++) {
      /* If this symbol codes a leaf node,  */
      if (index[j] <= 0) {
        if (len == dist[j] && path[j] == bits) {
          symb = j;
          break;
        }
      }
      else {
        if (len > dist[j] && path[j] == bits >> (len - dist[j])) {
          symb = j;
          break;
        }
      }
    }
    OD_ASSERT(symb != -1);
    od_ec_encode_cdf_q15(&w->ec, symb, cdf, nsymbs);
    bits &= (1 << (len - dist[symb])) - 1;
    len -= dist[symb];
  }
  while (len);
}
#endif

static INLINE void vp10_write_tree(vpx_writer *w, const vpx_tree_index *tree,
                                   const vpx_prob *probs, int bits, int len,
                                   vpx_tree_index i) {
  do {
    const int bit = (bits >> --len) & 1;
    vpx_write(w, bit, probs[i >> 1]);
    i = tree[i + bit];
  } while (len);
}

static INLINE void vp10_write_token(vpx_writer *w, const vpx_tree_index *tree,
                                    const vpx_prob *probs,
                                    const struct vp10_token *token) {
#if CONFIG_DAALA_EC
  daala_write_tree(w, tree, probs, token->value, token->len, 0);
#else
  vp10_write_tree(w, tree, probs, token->value, token->len, 0);
#endif
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // VP10_ENCODER_TREEWRITER_H_
