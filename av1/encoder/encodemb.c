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

#include "./av1_rtcd.h"
#include "./aom_config.h"
#include "./aom_dsp_rtcd.h"

#include "aom_dsp/quantize.h"
#include "aom_mem/aom_mem.h"
#include "aom_ports/mem.h"

#include "av1/common/idct.h"
#include "av1/common/reconinter.h"
#include "av1/common/reconintra.h"
#include "av1/common/scan.h"
#include "av1/common/partition.h"

#include "av1/encoder/encodemb.h"
#include "av1/encoder/rd.h"
#include "av1/encoder/tokenize.h"

#if CONFIG_PVQ
#include "av1/encoder/encint.h"
#include "av1/encoder/pvq_encoder.h"
#endif

struct optimize_ctx {
  ENTROPY_CONTEXT ta[MAX_MB_PLANE][16];
  ENTROPY_CONTEXT tl[MAX_MB_PLANE][16];
};

void av1_subtract_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane) {
  struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &x->e_mbd.plane[plane];
  const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize, pd);
  const int bw = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  const int bh = 4 * num_4x4_blocks_high_lookup[plane_bsize];

#if CONFIG_AOM_HIGHBITDEPTH
  if (x->e_mbd.cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    aom_highbd_subtract_block(bh, bw, p->src_diff, bw, p->src.buf,
                              p->src.stride, pd->dst.buf, pd->dst.stride,
                              x->e_mbd.bd);
    return;
  }
#endif  // CONFIG_AOM_HIGHBITDEPTH
  aom_subtract_block(bh, bw, p->src_diff, bw, p->src.buf, p->src.stride,
                     pd->dst.buf, pd->dst.stride);
}

#define RDTRUNC(RM, DM, R, D)                        \
  (((1 << (AV1_PROB_COST_SHIFT - 1)) + (R) * (RM)) & \
   ((1 << AV1_PROB_COST_SHIFT) - 1))

typedef struct av1_token_state {
  int rate;
  int error;
  int next;
  int16_t token;
  short qc;
} av1_token_state;

// TODO(jimbankoski): experiment to find optimal RD numbers.
static const int plane_rd_mult[PLANE_TYPES] = { 4, 2 };

#define UPDATE_RD_COST()                                \
  {                                                     \
    rd_cost0 = RDCOST(rdmult, rddiv, rate0, error0);    \
    rd_cost1 = RDCOST(rdmult, rddiv, rate1, error1);    \
    if (rd_cost0 == rd_cost1) {                         \
      rd_cost0 = RDTRUNC(rdmult, rddiv, rate0, error0); \
      rd_cost1 = RDTRUNC(rdmult, rddiv, rate1, error1); \
    }                                                   \
  }

#if !CONFIG_PVQ
// This function is a place holder for now but may ultimately need
// to scan previous tokens to work out the correct context.
static int trellis_get_coeff_context(const int16_t *scan, const int16_t *nb,
                                     int idx, int token, uint8_t *token_cache) {
  int bak = token_cache[scan[idx]], pt;
  token_cache[scan[idx]] = av1_pt_energy_class[token];
  pt = get_coef_context(nb, token_cache, idx + 1);
  token_cache[scan[idx]] = bak;
  return pt;
}

static int optimize_b(MACROBLOCK *mb, int plane, int block, TX_SIZE tx_size,
                      int ctx) {
  MACROBLOCKD *const xd = &mb->e_mbd;
  struct macroblock_plane *const p = &mb->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  const int ref = is_inter_block(&xd->mi[0]->mbmi);
  av1_token_state tokens[1025][2];
  unsigned best_index[1025][2];
  uint8_t token_cache[1024];
  const tran_low_t *const coeff = BLOCK_OFFSET(mb->plane[plane].coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  const int eob = p->eobs[block];
  const PLANE_TYPE type = pd->plane_type;
  const int default_eob = 16 << (tx_size << 1);
  const int mul = 1 + (tx_size == TX_32X32);
#if CONFIG_AOM_QM
  int seg_id = xd->mi[0]->mbmi.segment_id;
  int is_intra = !is_inter_block(&xd->mi[0]->mbmi);
  const qm_val_t *iqmatrix = pd->seg_iqmatrix[seg_id][is_intra][tx_size];
#endif
  const int16_t *dequant_ptr = pd->dequant;
  const uint8_t *const band_translate = get_band_translate(tx_size);
  TX_TYPE tx_type = get_tx_type(type, xd, block);
  const scan_order *const so = get_scan(tx_size, tx_type);
  const int16_t *const scan = so->scan;
  const int16_t *const nb = so->neighbors;
  int next = eob, sz = 0;
  int64_t rdmult = mb->rdmult * plane_rd_mult[type], rddiv = mb->rddiv;
  int64_t rd_cost0, rd_cost1;
  int rate0, rate1, error0, error1;
  int16_t t0, t1;
  EXTRABIT e0;
  int best, band, pt, i, final_eob;
#if CONFIG_AOM_HIGHBITDEPTH
  const int *cat6_high_cost = av1_get_high_cost_table(xd->bd);
#else
  const int *cat6_high_cost = av1_get_high_cost_table(8);
#endif

  assert((!type && !plane) || (type && plane));
  assert(eob <= default_eob);

  /* Now set up a Viterbi trellis to evaluate alternative roundings. */
  if (!ref) rdmult = (rdmult * 9) >> 4;

  /* Initialize the sentinel node of the trellis. */
  tokens[eob][0].rate = 0;
  tokens[eob][0].error = 0;
  tokens[eob][0].next = default_eob;
  tokens[eob][0].token = EOB_TOKEN;
  tokens[eob][0].qc = 0;
  tokens[eob][1] = tokens[eob][0];

  for (i = 0; i < eob; i++)
    token_cache[scan[i]] = av1_pt_energy_class[av1_get_token(qcoeff[scan[i]])];

  for (i = eob; i-- > 0;) {
    int base_bits, d2, dx;

    const int rc = scan[i];
#if CONFIG_AOM_QM
    int iwt = iqmatrix[rc];
#endif
    int x = qcoeff[rc];
    /* Only add a trellis state for non-zero coefficients. */
    if (x) {
      int shortcut = 0;
      error0 = tokens[next][0].error;
      error1 = tokens[next][1].error;
      /* Evaluate the first possibility for this state. */
      rate0 = tokens[next][0].rate;
      rate1 = tokens[next][1].rate;
      av1_get_token_extra(x, &t0, &e0);
      /* Consider both possible successor states. */
      if (next < default_eob) {
        band = band_translate[i + 1];
        pt = trellis_get_coeff_context(scan, nb, i, t0, token_cache);
        rate0 +=
            mb->token_costs[tx_size][type][ref][band][0][pt][tokens[next][0]
                                                                 .token];
        rate1 +=
            mb->token_costs[tx_size][type][ref][band][0][pt][tokens[next][1]
                                                                 .token];
      }
      UPDATE_RD_COST();
      /* And pick the best. */
      best = rd_cost1 < rd_cost0;
      base_bits = av1_get_cost(t0, e0, cat6_high_cost);
      dx = mul * (dqcoeff[rc] - coeff[rc]);
#if CONFIG_AOM_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        dx >>= xd->bd - 8;
      }
#endif  // CONFIG_AOM_HIGHBITDEPTH
      d2 = dx * dx;
      tokens[i][0].rate = base_bits + (best ? rate1 : rate0);
      tokens[i][0].error = d2 + (best ? error1 : error0);
      tokens[i][0].next = next;
      tokens[i][0].token = t0;
      tokens[i][0].qc = x;
      best_index[i][0] = best;

      /* Evaluate the second possibility for this state. */
      rate0 = tokens[next][0].rate;
      rate1 = tokens[next][1].rate;

#if CONFIG_AOM_QM
      if ((abs(x) * dequant_ptr[rc != 0] * iwt >
           ((abs(coeff[rc]) * mul) << AOM_QM_BITS)) &&
          (abs(x) * dequant_ptr[rc != 0] * iwt <
           ((abs(coeff[rc]) * mul + dequant_ptr[rc != 0]) << AOM_QM_BITS)))
#else
      if ((abs(x) * dequant_ptr[rc != 0] > abs(coeff[rc]) * mul) &&
          (abs(x) * dequant_ptr[rc != 0] <
           abs(coeff[rc]) * mul + dequant_ptr[rc != 0]))
#endif
        shortcut = 1;
      else
        shortcut = 0;

      if (shortcut) {
        sz = -(x < 0);
        x -= 2 * sz + 1;
      }

      /* Consider both possible successor states. */
      if (!x) {
        /* If we reduced this coefficient to zero, check to see if
         *  we need to move the EOB back here.
         */
        t0 = tokens[next][0].token == EOB_TOKEN ? EOB_TOKEN : ZERO_TOKEN;
        t1 = tokens[next][1].token == EOB_TOKEN ? EOB_TOKEN : ZERO_TOKEN;
        e0 = 0;
      } else {
        av1_get_token_extra(x, &t0, &e0);
        t1 = t0;
      }
      if (next < default_eob) {
        band = band_translate[i + 1];
        if (t0 != EOB_TOKEN) {
          pt = trellis_get_coeff_context(scan, nb, i, t0, token_cache);
          rate0 +=
              mb->token_costs[tx_size][type][ref][band][!x][pt][tokens[next][0]
                                                                    .token];
        }
        if (t1 != EOB_TOKEN) {
          pt = trellis_get_coeff_context(scan, nb, i, t1, token_cache);
          rate1 +=
              mb->token_costs[tx_size][type][ref][band][!x][pt][tokens[next][1]
                                                                    .token];
        }
      }

      UPDATE_RD_COST();
      /* And pick the best. */
      best = rd_cost1 < rd_cost0;
      base_bits = av1_get_cost(t0, e0, cat6_high_cost);

      if (shortcut) {
#if CONFIG_AOM_HIGHBITDEPTH
        if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
          dx -= ((dequant_ptr[rc != 0] >> (xd->bd - 8)) + sz) ^ sz;
        } else {
          dx -= (dequant_ptr[rc != 0] + sz) ^ sz;
        }
#else
        dx -= (dequant_ptr[rc != 0] + sz) ^ sz;
#endif  // CONFIG_AOM_HIGHBITDEPTH
        d2 = dx * dx;
      }

      tokens[i][1].rate = base_bits + (best ? rate1 : rate0);
      tokens[i][1].error = d2 + (best ? error1 : error0);
      tokens[i][1].next = next;
      tokens[i][1].token = best ? t1 : t0;
      tokens[i][1].qc = x;
      best_index[i][1] = best;
      /* Finally, make this the new head of the trellis. */
      next = i;
    } else {
      /* There's no choice to make for a zero coefficient, so we don't
       *  add a new trellis node, but we do need to update the costs.
       */
      band = band_translate[i + 1];
      t0 = tokens[next][0].token;
      t1 = tokens[next][1].token;
      /* Update the cost of each path if we're past the EOB token. */
      if (t0 != EOB_TOKEN) {
        tokens[next][0].rate +=
            mb->token_costs[tx_size][type][ref][band][1][0][t0];
        tokens[next][0].token = ZERO_TOKEN;
      }
      if (t1 != EOB_TOKEN) {
        tokens[next][1].rate +=
            mb->token_costs[tx_size][type][ref][band][1][0][t1];
        tokens[next][1].token = ZERO_TOKEN;
      }
      best_index[i][0] = best_index[i][1] = 0;
      /* Don't update next, because we didn't add a new node. */
    }
  }

  /* Now pick the best path through the whole trellis. */
  band = band_translate[i + 1];
  rate0 = tokens[next][0].rate;
  rate1 = tokens[next][1].rate;
  error0 = tokens[next][0].error;
  error1 = tokens[next][1].error;
  t0 = tokens[next][0].token;
  t1 = tokens[next][1].token;
  rate0 += mb->token_costs[tx_size][type][ref][band][0][ctx][t0];
  rate1 += mb->token_costs[tx_size][type][ref][band][0][ctx][t1];
  UPDATE_RD_COST();
  best = rd_cost1 < rd_cost0;
  final_eob = -1;
  memset(qcoeff, 0, sizeof(*qcoeff) * (16 << (tx_size * 2)));
  memset(dqcoeff, 0, sizeof(*dqcoeff) * (16 << (tx_size * 2)));
  for (i = next; i < eob; i = next) {
    const int x = tokens[i][best].qc;
    const int rc = scan[i];
#if CONFIG_AOM_QM
    const int iwt = iqmatrix[rc];
    const int dequant =
        (dequant_ptr[rc != 0] * iwt + (1 << (AOM_QM_BITS - 1))) >> AOM_QM_BITS;
#endif
    if (x) {
      final_eob = i;
    }

    qcoeff[rc] = x;
#if CONFIG_AOM_QM
    dqcoeff[rc] = (x * dequant) / mul;
#else
    dqcoeff[rc] = (x * dequant_ptr[rc != 0]) / mul;
#endif

    next = tokens[i][best].next;
    best = best_index[i][best];
  }
  final_eob++;

  mb->plane[plane].eobs[block] = final_eob;
  return final_eob;
}
#endif

static INLINE void fdct32x32(int rd_transform, const int16_t *src,
                             tran_low_t *dst, int src_stride) {
  if (rd_transform)
    aom_fdct32x32_rd(src, dst, src_stride);
  else
    aom_fdct32x32(src, dst, src_stride);
}

#if CONFIG_AOM_HIGHBITDEPTH
static INLINE void highbd_fdct32x32(int rd_transform, const int16_t *src,
                                    tran_low_t *dst, int src_stride) {
  if (rd_transform)
    aom_highbd_fdct32x32_rd(src, dst, src_stride);
  else
    aom_highbd_fdct32x32(src, dst, src_stride);
}
#endif  // CONFIG_AOM_HIGHBITDEPTH

void av1_fwd_txfm_4x4(const int16_t *src_diff, tran_low_t *coeff,
                      int diff_stride, TX_TYPE tx_type, int lossless) {
  if (lossless) {
    av1_fwht4x4(src_diff, coeff, diff_stride);
  } else {
    switch (tx_type) {
      case DCT_DCT: aom_fdct4x4(src_diff, coeff, diff_stride); break;
      case ADST_DCT:
      case DCT_ADST:
      case ADST_ADST: av1_fht4x4(src_diff, coeff, diff_stride, tx_type); break;
      default: assert(0); break;
    }
  }
}

void fwd_txfm_8x8(const int16_t *src_diff, tran_low_t *coeff,
                         int diff_stride, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST: av1_fht8x8(src_diff, coeff, diff_stride, tx_type); break;
    default: assert(0); break;
  }
}

void fwd_txfm_16x16(const int16_t *src_diff, tran_low_t *coeff,
                           int diff_stride, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST: av1_fht16x16(src_diff, coeff, diff_stride, tx_type); break;
    default: assert(0); break;
  }
}

void fwd_txfm_32x32(int rd_transform, const int16_t *src_diff,
                           tran_low_t *coeff, int diff_stride,
                           TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT: fdct32x32(rd_transform, src_diff, coeff, diff_stride); break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST: assert(0); break;
    default: assert(0); break;
  }
}

#if CONFIG_AOM_HIGHBITDEPTH
void av1_highbd_fwd_txfm_4x4(const int16_t *src_diff, tran_low_t *coeff,
                             int diff_stride, TX_TYPE tx_type, int lossless) {
  if (lossless) {
    assert(tx_type == DCT_DCT);
    av1_highbd_fwht4x4(src_diff, coeff, diff_stride);
  } else {
    switch (tx_type) {
      case DCT_DCT: aom_highbd_fdct4x4(src_diff, coeff, diff_stride); break;
      case ADST_DCT:
      case DCT_ADST:
      case ADST_ADST:
        av1_highbd_fht4x4(src_diff, coeff, diff_stride, tx_type);
        break;
      default: assert(0); break;
    }
  }
}

static void highbd_fwd_txfm_8x8(const int16_t *src_diff, tran_low_t *coeff,
                                int diff_stride, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT: aom_highbd_fdct8x8(src_diff, coeff, diff_stride); break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      av1_highbd_fht8x8(src_diff, coeff, diff_stride, tx_type);
      break;
    default: assert(0); break;
  }
}

static void highbd_fwd_txfm_16x16(const int16_t *src_diff, tran_low_t *coeff,
                                  int diff_stride, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT: aom_highbd_fdct16x16(src_diff, coeff, diff_stride); break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      av1_highbd_fht16x16(src_diff, coeff, diff_stride, tx_type);
      break;
    default: assert(0); break;
  }
}

static void highbd_fwd_txfm_32x32(int rd_transform, const int16_t *src_diff,
                                  tran_low_t *coeff, int diff_stride,
                                  TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
      highbd_fdct32x32(rd_transform, src_diff, coeff, diff_stride);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST: assert(0); break;
    default: assert(0); break;
  }
}
#endif  // CONFIG_AOM_HIGHBITDEPTH

void av1_xform_quant_fp(MACROBLOCK *x, int plane, int block, int blk_row,
                        int blk_col, BLOCK_SIZE plane_bsize, TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
#if !CONFIG_PVQ
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
#else
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
#endif
  PLANE_TYPE plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
  TX_TYPE tx_type = get_tx_type(plane_type, xd, block);
  const scan_order *const scan_order = get_scan(tx_size, tx_type);
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  int seg_id = xd->mi[0]->mbmi.segment_id;
#if CONFIG_AOM_QM
  int is_intra = !is_inter_block(&xd->mi[0]->mbmi);
  const qm_val_t *qmatrix = pd->seg_qmatrix[seg_id][is_intra][tx_size];
  const qm_val_t *iqmatrix = pd->seg_iqmatrix[seg_id][is_intra][tx_size];
#endif

#if !CONFIG_PVQ
  const int16_t *src_diff;
  src_diff = &p->src_diff[4 * (blk_row * diff_stride + blk_col)];
#else
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  tran_low_t *ref_coeff = BLOCK_OFFSET(pd->pvq_ref_coeff, block);
  uint8_t *src, *dst;
  int16_t *src_int16, *pred;
  const int src_stride = p->src.stride;
  const int dst_stride = pd->dst.stride;
  int tx_blk_size;
  int i, j;
  int skip = 1;
  PVQ_INFO *pvq_info;

  (void) scan_order;
  (void) qcoeff;
  pvq_info = &x->pvq[block][plane];
  dst = &pd->dst.buf[4 * (blk_row * dst_stride + blk_col)];
  src = &p->src.buf[4 * (blk_row * src_stride + blk_col)];
  src_int16 = &p->src_int16[4 * (blk_row * diff_stride + blk_col)];
  pred = &pd->pred[4 * (blk_row * diff_stride + blk_col)];
  // transform block size in pixels
  tx_blk_size = 1 << (tx_size + 2);

  // copy uint8 orig and predicted block to int16 buffer
  // in order to use existing VP10 transform functions
  for (j = 0; j < tx_blk_size; j++)
    for (i = 0; i < tx_blk_size; i++) {
      src_int16[diff_stride * j + i] = src[src_stride * j + i];
      pred[diff_stride * j + i] = dst[dst_stride * j + i];
    }
#endif

#if CONFIG_AOM_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
      case TX_32X32:
        highbd_fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
        av1_highbd_quantize_fp_32x32(
            coeff, 1024, x->skip_block, p->zbin, p->round_fp, p->quant_fp,
            p->quant_shift, qcoeff, dqcoeff, pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
            scan_order->iscan);
#else
            scan_order->iscan, qmatrix, iqmatrix);
#endif
        break;
      case TX_16X16:
        aom_highbd_fdct16x16(src_diff, coeff, diff_stride);
        av1_highbd_quantize_fp(coeff, 256, x->skip_block, p->zbin, p->round_fp,
                               p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                               pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                               scan_order->iscan);
#else
                               scan_order->iscan, qmatrix, iqmatrix);
#endif
        break;
      case TX_8X8:
        aom_highbd_fdct8x8(src_diff, coeff, diff_stride);
        av1_highbd_quantize_fp(coeff, 64, x->skip_block, p->zbin, p->round_fp,
                               p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                               pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                               scan_order->iscan);
#else
                               scan_order->iscan, qmatrix, iqmatrix);
#endif
        break;
      case TX_4X4:
        if (xd->lossless[seg_id]) {
          av1_highbd_fwht4x4(src_diff, coeff, diff_stride);
        } else {
          aom_highbd_fdct4x4(src_diff, coeff, diff_stride);
        }
        av1_highbd_quantize_fp(coeff, 16, x->skip_block, p->zbin, p->round_fp,
                               p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                               pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                               scan_order->iscan);
#else
                               scan_order->iscan, qmatrix, iqmatrix);
#endif
        break;
      default: assert(0);
    }
    return;
  }
#endif  // CONFIG_AOM_HIGHBITDEPTH

#if !CONFIG_PVQ
  switch (tx_size) {
    case TX_32X32:
      fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
      av1_quantize_fp_32x32(coeff, 1024, x->skip_block, p->zbin, p->round_fp,
                            p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                            pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                            scan_order->iscan);
#else
                            scan_order->iscan, qmatrix, iqmatrix);
#endif
      break;
    case TX_16X16:
      aom_fdct16x16(src_diff, coeff, diff_stride);
      av1_quantize_fp(coeff, 256, x->skip_block, p->zbin, p->round_fp,
                      p->quant_fp, p->quant_shift, qcoeff, dqcoeff, pd->dequant,
                      eob, scan_order->scan,
#if !CONFIG_AOM_QM
                      scan_order->iscan);
#else
                      scan_order->iscan, qmatrix, iqmatrix);
#endif
      break;
    case TX_8X8:
      av1_fdct8x8_quant(src_diff, diff_stride, coeff, 64, x->skip_block,
                        p->zbin, p->round_fp, p->quant_fp, p->quant_shift,
                        qcoeff, dqcoeff, pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                        scan_order->iscan);
#else
                        scan_order->iscan, qmatrix, iqmatrix);
#endif
      break;
    case TX_4X4:
      if (xd->lossless[seg_id]) {
        av1_fwht4x4(src_diff, coeff, diff_stride);
      } else {
        aom_fdct4x4(src_diff, coeff, diff_stride);
      }
      av1_quantize_fp(coeff, 16, x->skip_block, p->zbin, p->round_fp,
                      p->quant_fp, p->quant_shift, qcoeff, dqcoeff, pd->dequant,
                      eob, scan_order->scan,
#if !CONFIG_AOM_QM
                      scan_order->iscan);
#else
                      scan_order->iscan, qmatrix, iqmatrix);
#endif
      break;
    default: assert(0); break;
  }
#else//#if !CONFIG_PVQ
  switch (tx_size) {
    case TX_32X32:
      //forward transform of predicted image.
      fdct32x32(0, pred, ref_coeff, diff_stride);
      //forward transform of original image.
      fdct32x32(0, src_int16, coeff, diff_stride);
      break;
    case TX_16X16:
      aom_fdct16x16(pred, ref_coeff, diff_stride);
      aom_fdct16x16(src_int16, coeff, diff_stride);
      break;
    case TX_8X8:
      aom_fdct8x8(pred, ref_coeff, diff_stride);
      aom_fdct8x8(src_int16, coeff, diff_stride);
      break;
    case TX_4X4:
      if (xd->lossless[seg_id]) {
        av1_fwht4x4(pred, ref_coeff, diff_stride);
        av1_fwht4x4(src_int16, coeff, diff_stride);
      } else {
        aom_fdct4x4(pred, ref_coeff, diff_stride);
        aom_fdct4x4(src_int16, coeff, diff_stride);
      }
      break;
    default: assert(0); break;
  }

  // PVQ for inter mode block
  if (!x->skip_block)
    skip = pvq_encode_helper(&x->daala_enc,
                             coeff,          // target original vector
                             ref_coeff,      // reference vector
                             dqcoeff,        // de-quantized vector
                             eob,            // End of Block marker
                             pd->dequant,    // aom's quantizers
                             plane,          // image plane
                             tx_size,        // block size in log_2 - 2
                             &x->rate,       // rate measured
                             pvq_info);      // PVQ info for a block

  if (!skip)
    mbmi->skip = 0;
#endif//#if !CONFIG_PVQ
}

#if !CONFIG_PVQ
void av1_xform_quant_dc(MACROBLOCK *x, int plane, int block, int blk_row,
                        int blk_col, BLOCK_SIZE plane_bsize, TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  int seg_id = xd->mi[0]->mbmi.segment_id;
#if CONFIG_AOM_QM
  int is_intra = !is_inter_block(&xd->mi[0]->mbmi);
  const qm_val_t *qmatrix = pd->seg_qmatrix[seg_id][is_intra][tx_size];
  const qm_val_t *iqmatrix = pd->seg_iqmatrix[seg_id][is_intra][tx_size];
#endif

  const int16_t *src_diff;
  src_diff = &p->src_diff[4 * (blk_row * diff_stride + blk_col)];

#if CONFIG_AOM_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
      case TX_32X32:
        aom_highbd_fdct32x32_1(src_diff, coeff, diff_stride);
        aom_highbd_quantize_dc_32x32(coeff, x->skip_block, p->round,
                                     p->quant_fp[0], qcoeff, dqcoeff,
                                     pd->dequant[0],
#if !CONFIG_AOM_QM
                                     eob);
#else
                                     eob, qmatrix, iqmatrix);
#endif
        break;
      case TX_16X16:
        aom_highbd_fdct16x16_1(src_diff, coeff, diff_stride);
        aom_highbd_quantize_dc(coeff, 256, x->skip_block, p->round,
                               p->quant_fp[0], qcoeff, dqcoeff, pd->dequant[0],
#if !CONFIG_AOM_QM
                               eob);
#else
                               eob, qmatrix, iqmatrix);
#endif
        break;
      case TX_8X8:
        aom_highbd_fdct8x8_1(src_diff, coeff, diff_stride);
        aom_highbd_quantize_dc(coeff, 64, x->skip_block, p->round,
                               p->quant_fp[0], qcoeff, dqcoeff, pd->dequant[0],
#if !CONFIG_AOM_QM
                               eob);
#else
                               eob, qmatrix, iqmatrix);
#endif
        break;
      case TX_4X4:
        if (xd->lossless[seg_id]) {
          av1_highbd_fwht4x4(src_diff, coeff, diff_stride);
        } else {
          aom_highbd_fdct4x4(src_diff, coeff, diff_stride);
        }
        aom_highbd_quantize_dc(coeff, 16, x->skip_block, p->round,
                               p->quant_fp[0], qcoeff, dqcoeff, pd->dequant[0],
#if !CONFIG_AOM_QM
                               eob);
#else
                               eob, qmatrix, iqmatrix);
#endif
        break;
      default: assert(0);
    }
    return;
  }
#endif  // CONFIG_AOM_HIGHBITDEPTH

  switch (tx_size) {
    case TX_32X32:
      aom_fdct32x32_1(src_diff, coeff, diff_stride);
      aom_quantize_dc_32x32(coeff, x->skip_block, p->round, p->quant_fp[0],
                            qcoeff, dqcoeff, pd->dequant[0],
#if !CONFIG_AOM_QM
                            eob);
#else
                            eob, qmatrix, iqmatrix);
#endif
      break;
    case TX_16X16:
      aom_fdct16x16_1(src_diff, coeff, diff_stride);
      aom_quantize_dc(coeff, 256, x->skip_block, p->round, p->quant_fp[0],
                      qcoeff, dqcoeff, pd->dequant[0],
#if !CONFIG_AOM_QM
                      eob);
#else
                      eob, qmatrix, iqmatrix);
#endif
      break;
    case TX_8X8:
      aom_fdct8x8_1(src_diff, coeff, diff_stride);
      aom_quantize_dc(coeff, 64, x->skip_block, p->round, p->quant_fp[0],
                      qcoeff, dqcoeff, pd->dequant[0],
#if !CONFIG_AOM_QM
                      eob);
#else
                      eob, qmatrix, iqmatrix);
#endif
      break;
    case TX_4X4:
      if (xd->lossless[seg_id]) {
        av1_fwht4x4(src_diff, coeff, diff_stride);
      } else {
        aom_fdct4x4(src_diff, coeff, diff_stride);
      }
      aom_quantize_dc(coeff, 16, x->skip_block, p->round, p->quant_fp[0],
                      qcoeff, dqcoeff, pd->dequant[0],
#if !CONFIG_AOM_QM
                      eob);
#else
                      eob, qmatrix, iqmatrix);
#endif
      break;
    default: assert(0); break;
  }
}
#endif

void av1_xform_quant(MACROBLOCK *x, int plane, int block, int blk_row,
                     int blk_col, BLOCK_SIZE plane_bsize, TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
#if !CONFIG_PVQ
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
#else
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
#endif
  PLANE_TYPE plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
  TX_TYPE tx_type = get_tx_type(plane_type, xd, block);
  const scan_order *const scan_order = get_scan(tx_size, tx_type);
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  int seg_id = xd->mi[0]->mbmi.segment_id;
#if CONFIG_AOM_QM
  int is_intra = !is_inter_block(&xd->mi[0]->mbmi);
  const qm_val_t *qmatrix = pd->seg_qmatrix[seg_id][is_intra][tx_size];
  const qm_val_t *iqmatrix = pd->seg_iqmatrix[seg_id][is_intra][tx_size];
#endif

#if !CONFIG_PVQ
  const int16_t *src_diff;
  src_diff = &p->src_diff[4 * (blk_row * diff_stride + blk_col)];
#else
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  tran_low_t *ref_coeff = BLOCK_OFFSET(pd->pvq_ref_coeff, block);
  uint8_t *src, *dst;
  int16_t *src_int16, *pred;
  const int src_stride = p->src.stride;
  const int dst_stride = pd->dst.stride;
  int tx_blk_size;
  int i, j;
  int skip = 1;
  PVQ_INFO *pvq_info;

  (void) scan_order;
  (void) qcoeff;
  pvq_info = &x->pvq[block][plane];
  dst = &pd->dst.buf[4 * (blk_row * dst_stride + blk_col)];
  src = &p->src.buf[4 * (blk_row * src_stride + blk_col)];
  src_int16 = &p->src_int16[4 * (blk_row * diff_stride + blk_col)];
  pred = &pd->pred[4 * (blk_row * diff_stride + blk_col)];

  // transform block size in pixels
  tx_blk_size = 1 << (tx_size + 2);

  // copy uint8 orig and predicted block to int16 buffer
  // in order to use existing VP10 transform functions
  for (j = 0; j < tx_blk_size; j++)
    for (i = 0; i < tx_blk_size; i++) {
      src_int16[diff_stride * j + i] = src[src_stride * j + i];
      pred[diff_stride * j + i] = dst[dst_stride * j + i];
    }
#endif

#if CONFIG_AOM_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
      case TX_32X32:
        highbd_fwd_txfm_32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride,
                              tx_type);
        aom_highbd_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin,
                                    p->round, p->quant, p->quant_shift, qcoeff,
                                    dqcoeff, pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                                    scan_order->iscan);
#else
                                    scan_order->iscan, qmatrix, iqmatrix);
#endif
        break;
      case TX_16X16:
        highbd_fwd_txfm_16x16(src_diff, coeff, diff_stride, tx_type);
        aom_highbd_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                              scan_order->iscan);
#else
                              scan_order->iscan, qmatrix, iqmatrix);
#endif
        break;
      case TX_8X8:
        highbd_fwd_txfm_8x8(src_diff, coeff, diff_stride, tx_type);
        aom_highbd_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                              scan_order->iscan);
#else
                              scan_order->iscan, qmatrix, iqmatrix);
#endif
        break;
      case TX_4X4:
        av1_highbd_fwd_txfm_4x4(src_diff, coeff, diff_stride, tx_type,
                                xd->lossless[seg_id]);
        aom_highbd_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                              scan_order->iscan);
#else
                              scan_order->iscan, qmatrix, iqmatrix);
#endif
        break;
      default: assert(0);
    }
    return;
  }
#endif  // CONFIG_AOM_HIGHBITDEPTH

#if !CONFIG_PVQ
  switch (tx_size) {
    case TX_32X32:
      fwd_txfm_32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride, tx_type);
      aom_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin, p->round,
                           p->quant, p->quant_shift, qcoeff, dqcoeff,
                           pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                           scan_order->iscan);
#else
                           scan_order->iscan, qmatrix, iqmatrix);
#endif
      break;
    case TX_16X16:
      fwd_txfm_16x16(src_diff, coeff, diff_stride, tx_type);
      aom_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round, p->quant,
                     p->quant_shift, qcoeff, dqcoeff, pd->dequant, eob,
                     scan_order->scan,
#if !CONFIG_AOM_QM
                     scan_order->iscan);
#else
                     scan_order->iscan, qmatrix, iqmatrix);
#endif
      break;
    case TX_8X8:
      fwd_txfm_8x8(src_diff, coeff, diff_stride, tx_type);
      aom_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round, p->quant,
                     p->quant_shift, qcoeff, dqcoeff, pd->dequant, eob,
                     scan_order->scan,
#if !CONFIG_AOM_QM
                     scan_order->iscan);
#else
                     scan_order->iscan, qmatrix, iqmatrix);
#endif
      break;
    case TX_4X4:
      av1_fwd_txfm_4x4(src_diff, coeff, diff_stride, tx_type,
                       xd->lossless[seg_id]);
      aom_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round, p->quant,
                     p->quant_shift, qcoeff, dqcoeff, pd->dequant, eob,
                     scan_order->scan,
#if !CONFIG_AOM_QM
                     scan_order->iscan);
#else
                     scan_order->iscan, qmatrix, iqmatrix);
#endif
      break;
    default: assert(0); break;
  }
#else//#if !CONFIG_PVQ
  switch (tx_size) {
    case TX_32X32:
      //forward transform of predicted image.
      fwd_txfm_32x32(0, pred, ref_coeff, diff_stride, tx_type);
      //forward transform of original image.
      fwd_txfm_32x32(0, src_int16, coeff, diff_stride, tx_type);
      break;
    case TX_16X16:
      fwd_txfm_16x16(pred, ref_coeff, diff_stride, tx_type);
      fwd_txfm_16x16(src_int16, coeff, diff_stride, tx_type);
      break;
    case TX_8X8:
      fwd_txfm_8x8(pred, ref_coeff, diff_stride, tx_type);
      fwd_txfm_8x8(src_int16, coeff, diff_stride, tx_type);
      break;
    case TX_4X4:
      av1_fwd_txfm_4x4(pred, ref_coeff, diff_stride, tx_type,
                        xd->lossless[seg_id]);
      av1_fwd_txfm_4x4(src_int16, coeff, diff_stride, tx_type,
                        xd->lossless[seg_id]);
      break;
    default: assert(0); break;
  }

  // PVQ for inter mode block
  if (!x->skip_block)
    skip = pvq_encode_helper(&x->daala_enc,
                             coeff,          // target original vector
                             ref_coeff,      // reference vector
                             dqcoeff,        // de-quantized vector
                             eob,            // End of Block marker
                             pd->dequant,    // aom's quantizers
                             plane,          // image plane
                             tx_size,        // block size in log_2 - 2
                             &x->rate,       // rate measured
                             pvq_info);      // PVQ info for a block

  if (!skip)
    mbmi->skip = 0;
#endif//#if !CONFIG_PVQ
}

static void encode_block(int plane, int block, int blk_row, int blk_col,
                         BLOCK_SIZE plane_bsize, TX_SIZE tx_size, void *arg) {
  struct encode_b_args *const args = arg;
  MACROBLOCK *const x = args->x;
  MACROBLOCKD *const xd = &x->e_mbd;
  struct optimize_ctx *const ctx = args->ctx;
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint8_t *dst;
  ENTROPY_CONTEXT *a, *l;
  TX_TYPE tx_type = get_tx_type(pd->plane_type, xd, block);
#if CONFIG_PVQ
  int tx_blk_size;
  int i, j;
  PVQ_INFO *pvq_info = &x->pvq[block][plane];
#endif

  dst = &pd->dst.buf[4 * blk_row * pd->dst.stride + 4 * blk_col];
  a = &ctx->ta[plane][blk_col];
  l = &ctx->tl[plane][blk_row];

  // TODO(jingning): per transformed block zero forcing only enabled for
  // luma component. will integrate chroma components as well.
  if (x->zcoeff_blk[tx_size][block] && plane == 0) {
    p->eobs[block] = 0;
    *a = *l = 0;
    return;
  }

  if (x->quant_fp) {
    // Encoding process for rtc mode
#if !CONFIG_PVQ
    if (x->skip_txfm[0] == SKIP_TXFM_AC_DC && plane == 0) {
      // skip forward transform
      p->eobs[block] = 0;
      *a = *l = 0;
      return;
    } else {
      av1_xform_quant_fp(x, plane, block, blk_row, blk_col, plane_bsize,
                         tx_size);
    }
#else
      av1_xform_quant_fp(x, plane, block, blk_row, blk_col, plane_bsize,
                         tx_size);
#endif
  } else {
#if !CONFIG_PVQ
    if (max_txsize_lookup[plane_bsize] == tx_size) {
      int txfm_blk_index = (plane << 2) + (block >> (tx_size << 1));
      if (x->skip_txfm[txfm_blk_index] == SKIP_TXFM_NONE) {
        // full forward transform and quantization
        av1_xform_quant(x, plane, block, blk_row, blk_col, plane_bsize,
                        tx_size);
      } else if (x->skip_txfm[txfm_blk_index] == SKIP_TXFM_AC_ONLY) {
        // fast path forward transform and quantization
        av1_xform_quant_dc(x, plane, block, blk_row, blk_col, plane_bsize,
                           tx_size);
      } else {
        // skip forward transform
        p->eobs[block] = 0;
        *a = *l = 0;
        return;
      }
    } else {
      av1_xform_quant(x, plane, block, blk_row, blk_col, plane_bsize, tx_size);
    }
#else
      av1_xform_quant(x, plane, block, blk_row, blk_col, plane_bsize, tx_size);
#endif
  }

#if !CONFIG_PVQ
  if (x->optimize) {
    const int ctx = combine_entropy_contexts(*a, *l);
    *a = *l = optimize_b(x, plane, block, tx_size, ctx) > 0;
  } else {
    *a = *l = p->eobs[block] > 0;
  }

  if (p->eobs[block]) *(args->skip) = 0;

  if (p->eobs[block] == 0) return;
#else
  if (pvq_info->ac_dc_coded)
    *(args->skip) = 0;

  if (!pvq_info->ac_dc_coded) return;

  // transform block size in pixels
  tx_blk_size = 1 << (tx_size + 2);

  // Since av1 does not have separate function which does inverse transform
  // but av1_inv_txfm_add_*x*() also does addition of predicted image to
  // inverse transformed image,
  // pass blank dummy image to av1_inv_txfm_add_*x*(), i.e. set dst as zeros
    for (j=0; j < tx_blk_size; j++)
      for (i = 0; i < tx_blk_size; i++)
        dst[j * pd->dst.stride + i] -= dst[j *  pd->dst.stride + i];
#endif

#if CONFIG_AOM_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
      case TX_32X32:
        av1_highbd_inv_txfm_add_32x32(dqcoeff, dst, pd->dst.stride,
                                      p->eobs[block], xd->bd, tx_type);
        break;
      case TX_16X16:
        av1_highbd_inv_txfm_add_16x16(dqcoeff, dst, pd->dst.stride,
                                      p->eobs[block], xd->bd, tx_type);
        break;
      case TX_8X8:
        av1_highbd_inv_txfm_add_8x8(dqcoeff, dst, pd->dst.stride,
                                    p->eobs[block], xd->bd, tx_type);
        break;
      case TX_4X4:
        // this is like av1_short_idct4x4 but has a special case around eob<=1
        // which is significant (not just an optimization) for the lossless
        // case.
        av1_highbd_inv_txfm_add_4x4(dqcoeff, dst, pd->dst.stride,
                                    p->eobs[block], xd->bd, tx_type,
                                    xd->lossless[xd->mi[0]->mbmi.segment_id]);
        break;
      default: assert(0 && "Invalid transform size"); break;
    }

    return;
  }
#endif  // CONFIG_AOM_HIGHBITDEPTH
  switch (tx_size) {
    case TX_32X32:
      av1_inv_txfm_add_32x32(dqcoeff, dst, pd->dst.stride, p->eobs[block],
                             tx_type);
      break;
    case TX_16X16:
      av1_inv_txfm_add_16x16(dqcoeff, dst, pd->dst.stride, p->eobs[block],
                             tx_type);
      break;
    case TX_8X8:
      av1_inv_txfm_add_8x8(dqcoeff, dst, pd->dst.stride, p->eobs[block],
                           tx_type);
      break;
    case TX_4X4:
      // this is like av1_short_idct4x4 but has a special case around eob<=1
      // which is significant (not just an optimization) for the lossless
      // case.
      av1_inv_txfm_add_4x4(dqcoeff, dst, pd->dst.stride, p->eobs[block],
                           tx_type, xd->lossless[xd->mi[0]->mbmi.segment_id]);
      break;
    default: assert(0 && "Invalid transform size"); break;
  }
}

static void encode_block_pass1(int plane, int block, int blk_row, int blk_col,
                               BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                               void *arg) {
  MACROBLOCK *const x = (MACROBLOCK *)arg;
  MACROBLOCKD *const xd = &x->e_mbd;
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint8_t *dst;
  dst = &pd->dst.buf[4 * blk_row * pd->dst.stride + 4 * blk_col];

  av1_xform_quant(x, plane, block, blk_row, blk_col, plane_bsize, tx_size);

  if (p->eobs[block] > 0) {
#if CONFIG_AOM_HIGHBITDEPTH
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
      if (xd->lossless[0]) {
        av1_highbd_iwht4x4_add(dqcoeff, dst, pd->dst.stride, p->eobs[block],
                               xd->bd);
      } else {
        av1_highbd_idct4x4_add(dqcoeff, dst, pd->dst.stride, p->eobs[block],
                               xd->bd);
      }
      return;
    }
#endif  // CONFIG_AOM_HIGHBITDEPTH
    if (xd->lossless[0]) {
      av1_iwht4x4_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
    } else {
      av1_idct4x4_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
    }
  }
}

void av1_encode_sby_pass1(MACROBLOCK *x, BLOCK_SIZE bsize) {
  av1_subtract_plane(x, bsize, 0);
  av1_foreach_transformed_block_in_plane(&x->e_mbd, bsize, 0,
                                         encode_block_pass1, x);
}

void av1_encode_sb(MACROBLOCK *x, BLOCK_SIZE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  struct optimize_ctx ctx;
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  struct encode_b_args arg = { x, &ctx, &mbmi->skip };
  int plane;

  mbmi->skip = 1;

  if (x->skip) return;

  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
#if !CONFIG_PVQ
    av1_subtract_plane(x, bsize, plane);
#endif
    if (x->optimize) {
      const struct macroblockd_plane *const pd = &xd->plane[plane];
      const TX_SIZE tx_size = plane ? get_uv_tx_size(mbmi, pd) : mbmi->tx_size;
      av1_get_entropy_contexts(bsize, tx_size, pd, ctx.ta[plane],
                               ctx.tl[plane]);
    }

    av1_foreach_transformed_block_in_plane(xd, bsize, plane, encode_block,
                                           &arg);
  }
}

void av1_encode_block_intra(int plane, int block, int blk_row, int blk_col,
                            BLOCK_SIZE plane_bsize, TX_SIZE tx_size,
                            void *arg) {
  struct encode_b_args *const args = arg;
  MACROBLOCK *const x = args->x;
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
  tran_low_t *coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  PLANE_TYPE plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
  TX_TYPE tx_type = get_tx_type(plane_type, xd, block);
  const scan_order *const scan_order = get_scan(tx_size, tx_type);
  PREDICTION_MODE mode;
  const int bwl = b_width_log2_lookup[plane_bsize];
  const int bhl = b_height_log2_lookup[plane_bsize];
  const int diff_stride = 4 * (1 << bwl);
  uint8_t *src, *dst;
  uint16_t *eob = &p->eobs[block];
  int seg_id = xd->mi[0]->mbmi.segment_id;
#if CONFIG_AOM_QM
  int is_intra = !is_inter_block(&xd->mi[0]->mbmi);
  const qm_val_t *qmatrix = pd->seg_qmatrix[seg_id][is_intra][tx_size];
  const qm_val_t *iqmatrix = pd->seg_iqmatrix[seg_id][is_intra][tx_size];
#endif
  const int src_stride = p->src.stride;
  const int dst_stride = pd->dst.stride;
#if !CONFIG_PVQ
  int16_t *src_diff;
  src_diff = &p->src_diff[4 * (blk_row * diff_stride + blk_col)];
#else
   tran_low_t *ref_coeff = BLOCK_OFFSET(pd->pvq_ref_coeff, block);
  int16_t *src_int16;
  int tx_blk_size;
  int i, j;
  int16_t *pred = &pd->pred[4 * (blk_row * diff_stride + blk_col)];
  int skip = 1;
  PVQ_INFO *pvq_info = &x->pvq[block][plane];

  (void) scan_order;
  (void) qcoeff;
  src_int16 = &p->src_int16[4 * (blk_row * diff_stride + blk_col)];
#endif
  dst = &pd->dst.buf[4 * (blk_row * dst_stride + blk_col)];
  src = &p->src.buf[4 * (blk_row * src_stride + blk_col)];
  mode = plane == 0 ? get_y_mode(xd->mi[0], block) : mbmi->uv_mode;
  av1_predict_intra_block(xd, bwl, bhl, tx_size, mode, dst, dst_stride, dst,
                          dst_stride, blk_col, blk_row, plane);

#if CONFIG_AOM_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
      case TX_32X32:
        aom_highbd_subtract_block(32, 32, src_diff, diff_stride, src,
                                  src_stride, dst, dst_stride, xd->bd);
        highbd_fwd_txfm_32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride,
                              tx_type);
        aom_highbd_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin,
                                    p->round, p->quant, p->quant_shift, qcoeff,
                                    dqcoeff, pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                                    scan_order->iscan);
#else
                                    scan_order->iscan, qmatrix, iqmatrix);
#endif
        if (*eob)
          av1_highbd_inv_txfm_add_32x32(dqcoeff, dst, dst_stride, *eob, xd->bd,
                                        tx_type);
        break;
      case TX_16X16:
        aom_highbd_subtract_block(16, 16, src_diff, diff_stride, src,
                                  src_stride, dst, dst_stride, xd->bd);
        highbd_fwd_txfm_16x16(src_diff, coeff, diff_stride, tx_type);
        aom_highbd_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                              scan_order->iscan);
#else
                              scan_order->iscan, qmatrix, iqmatrix);
#endif
        if (*eob)
          av1_highbd_inv_txfm_add_16x16(dqcoeff, dst, dst_stride, *eob, xd->bd,
                                        tx_type);
        break;
      case TX_8X8:
        aom_highbd_subtract_block(8, 8, src_diff, diff_stride, src, src_stride,
                                  dst, dst_stride, xd->bd);
        highbd_fwd_txfm_8x8(src_diff, coeff, diff_stride, tx_type);
        aom_highbd_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                              scan_order->iscan);
#else
                              scan_order->iscan, qmatrix, iqmatrix);
#endif
        if (*eob)
          av1_highbd_inv_txfm_add_8x8(dqcoeff, dst, dst_stride, *eob, xd->bd,
                                      tx_type);
        break;
      case TX_4X4:
        aom_highbd_subtract_block(4, 4, src_diff, diff_stride, src, src_stride,
                                  dst, dst_stride, xd->bd);
        av1_highbd_fwd_txfm_4x4(src_diff, coeff, diff_stride, tx_type,
                                xd->lossless[seg_id]);
        aom_highbd_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                              scan_order->iscan);
#else
                              scan_order->iscan, qmatrix, iqmatrix);
#endif
        if (*eob)
          // this is like av1_short_idct4x4 but has a special case around
          // eob<=1 which is significant (not just an optimization) for the
          // lossless case.
          av1_highbd_inv_txfm_add_4x4(dqcoeff, dst, dst_stride, *eob, xd->bd,
                                      tx_type, xd->lossless[seg_id]);
        break;
      default: assert(0); return;
    }
    if (*eob) *(args->skip) = 0;
    return;
  }
#endif  // CONFIG_AOM_HIGHBITDEPTH

#if !CONFIG_PVQ
  switch (tx_size) {
    case TX_32X32:
      aom_subtract_block(32, 32, src_diff, diff_stride, src, src_stride, dst,
                         dst_stride);
      fwd_txfm_32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride, tx_type);
      aom_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin, p->round,
                           p->quant, p->quant_shift, qcoeff, dqcoeff,
                           pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                           scan_order->iscan);
#else
                           scan_order->iscan, qmatrix, iqmatrix);
#endif
      if (*eob) av1_inv_txfm_add_32x32(dqcoeff, dst, dst_stride, *eob, tx_type);
      break;
    case TX_16X16:
      aom_subtract_block(16, 16, src_diff, diff_stride, src, src_stride, dst,
                         dst_stride);
      fwd_txfm_16x16(src_diff, coeff, diff_stride, tx_type);
      aom_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round, p->quant,
                     p->quant_shift, qcoeff, dqcoeff, pd->dequant, eob,
                     scan_order->scan,
#if !CONFIG_AOM_QM
                     scan_order->iscan);
#else
                     scan_order->iscan, qmatrix, iqmatrix);
#endif
      if (*eob) av1_inv_txfm_add_16x16(dqcoeff, dst, dst_stride, *eob, tx_type);
      break;
    case TX_8X8:
      aom_subtract_block(8, 8, src_diff, diff_stride, src, src_stride, dst,
                         dst_stride);
      fwd_txfm_8x8(src_diff, coeff, diff_stride, tx_type);
      aom_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round, p->quant,
                     p->quant_shift, qcoeff, dqcoeff, pd->dequant, eob,
                     scan_order->scan,
#if !CONFIG_AOM_QM
                     scan_order->iscan);
#else
                     scan_order->iscan, qmatrix, iqmatrix);
#endif
      if (*eob) av1_inv_txfm_add_8x8(dqcoeff, dst, dst_stride, *eob, tx_type);
      break;
    case TX_4X4:
      aom_subtract_block(4, 4, src_diff, diff_stride, src, src_stride, dst,
                         dst_stride);
      av1_fwd_txfm_4x4(src_diff, coeff, diff_stride, tx_type,
                       xd->lossless[seg_id]);
      aom_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round, p->quant,
                     p->quant_shift, qcoeff, dqcoeff, pd->dequant, eob,
                     scan_order->scan,
#if !CONFIG_AOM_QM
                     scan_order->iscan);
#else
                     scan_order->iscan, qmatrix, iqmatrix);
#endif

      if (*eob) {
        // this is like av1_short_idct4x4 but has a special case around eob<=1
        // which is significant (not just an optimization) for the lossless
        // case.
        av1_inv_txfm_add_4x4(dqcoeff, dst, dst_stride, *eob, tx_type,
                             xd->lossless[seg_id]);
      }
      break;
    default: assert(0); break;
  }
#else//#if !CONFIG_PVQ
  // transform block size in pixels
  tx_blk_size = 1 << (tx_size + 2);

  // copy uint8 orig and predicted block to int16 buffer
  // in order to use existing VP10 transform functions
  for (j = 0; j < tx_blk_size; j++)
    for (i = 0; i < tx_blk_size; i++) {
      src_int16[diff_stride * j + i] = src[src_stride * j + i];
      pred[diff_stride * j + i] = dst[dst_stride * j + i];
    }

  switch (tx_size) {
    case TX_32X32:
      //forward transform of predicted image.
      fwd_txfm_32x32(0, pred, ref_coeff, diff_stride,
                     tx_type);
      //forward transform of original image.
      fwd_txfm_32x32(0, src_int16, coeff, diff_stride,
                     tx_type);
      break;
    case TX_16X16:
      fwd_txfm_16x16(pred, ref_coeff, diff_stride, tx_type);
      fwd_txfm_16x16(src_int16, coeff, diff_stride, tx_type);
      break;
    case TX_8X8:
      fwd_txfm_8x8(pred, ref_coeff, diff_stride, tx_type);
      fwd_txfm_8x8(src_int16, coeff, diff_stride, tx_type);
      break;
    case TX_4X4:
      av1_fwd_txfm_4x4(pred, ref_coeff, diff_stride, tx_type,
                        xd->lossless[seg_id]);
      av1_fwd_txfm_4x4(src_int16, coeff, diff_stride, tx_type,
                        xd->lossless[seg_id]);
      break;
    default: assert(0); break;
  }

  // PVQ for intra mode block
  if (!x->skip_block)
    skip = pvq_encode_helper(&x->daala_enc,
                             coeff,          // target original vector
                             ref_coeff,      // reference vector
                             dqcoeff,        // de-quantized vector
                             eob,            // End of Block marker
                             pd->dequant,    // aom's quantizers
                             plane,          // image plane
                             tx_size,        // block size in log_2 - 2
                             &x->rate,       // rate measured
                             pvq_info);       // PVQ info for a block

  if (!skip)
    mbmi->skip = 0;

  // Since av1 does not have separate function which does inverse transform
  // but av1_inv_txfm_add_*x*() also does addition of predicted image to
  // inverse transformed image,
  // pass blank dummy image to av1_inv_txfm_add_*x*(), i.e. set dst as zeros

  if (!skip) {
    for (j=0; j < tx_blk_size; j++)
      for (i = 0; i < tx_blk_size; i++)
        dst[j * dst_stride + i] -= dst[j * dst_stride + i];

    switch (tx_size) {
      case TX_32X32:
       av1_inv_txfm_add_32x32(dqcoeff, dst, dst_stride, *eob, tx_type);
        break;
      case TX_16X16:
        av1_inv_txfm_add_16x16(dqcoeff, dst, dst_stride, *eob, tx_type);
        break;
      case TX_8X8:
        av1_inv_txfm_add_8x8(dqcoeff, dst, dst_stride, *eob, tx_type);
        break;
      case TX_4X4:
        // this is like av1_short_idct4x4 but has a special case around eob<=1
        // which is significant (not just an optimization) for the lossless
        // case.
        av1_inv_txfm_add_4x4(dqcoeff, dst, dst_stride, *eob, tx_type,
                              xd->lossless[seg_id]);
        break;
      default: assert(0); break;
    }
  }
#endif//#if !CONFIG_PVQ

#if !CONFIG_PVQ
  if (*eob) *(args->skip) = 0;
#else
  // Note : *(args->skip) == mbmi->skip
#endif
}

void av1_encode_intra_block_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane) {
  const MACROBLOCKD *const xd = &x->e_mbd;
  struct encode_b_args arg = { x, NULL, &xd->mi[0]->mbmi.skip };

  av1_foreach_transformed_block_in_plane(xd, bsize, plane,
                                         av1_encode_block_intra, &arg);
}

#if CONFIG_PVQ
int pvq_encode_helper(daala_enc_ctx *daala_enc,
    tran_low_t *const coeff, tran_low_t *ref_coeff,
    tran_low_t *const dqcoeff,
    uint16_t *eob, const int16_t *quant,
    int plane, int tx_size, int *rate, PVQ_INFO *pvq_info) {
  const int tx_blk_size = 1 << (tx_size + 2);
  int skip;
  // TODO: Enable this later, if pvq_qm_q4 is available in AOM.
  //int pvq_dc_quant = OD_MAXI(1,
  //  quant * daala_enc->state.pvq_qm_q4[plane][od_qm_get_index(tx_size, 0)] >> 4);
  int quant_shift = tx_size == TX_32X32 ? 1 : 0;
  int pvq_dc_quant = OD_MAXI(1, quant[0] >> quant_shift);
  int tell;
  int has_dc_skip = 1;
  int i;
  int off = od_qm_offset(tx_size, plane ? 1 : 0);

  DECLARE_ALIGNED(16, int16_t, coeff_pvq[64 * 64]);
  DECLARE_ALIGNED(16, int16_t, ref_coeff_pvq[64 * 64]);
  DECLARE_ALIGNED(16, int16_t, dqcoeff_pvq[64 * 64]);

  DECLARE_ALIGNED(16, int32_t, in_int32[64 * 64]);
  DECLARE_ALIGNED(16, int32_t, ref_int32[64 * 64]);
  DECLARE_ALIGNED(16, int32_t, out_int32[64 * 64]);

  *eob = 0;

  tell = od_ec_enc_tell(&daala_enc->ec);

  /*{// for debugging, to match with decoder side ref vector.
  int i;
  for (i=0; i<tx_blk_size*tx_blk_size; i++)
    pvq_info->ref_coeff[i] = ref_coeff[i];
  }*/

  // Change coefficient ordering for pvq encoding.
  od_raster_to_coding_order(coeff_pvq, tx_blk_size, coeff, tx_blk_size);
  od_raster_to_coding_order(ref_coeff_pvq, tx_blk_size, ref_coeff, tx_blk_size);

  //copy int16 inputs to int32
  for (i=0; i < tx_blk_size * tx_blk_size; i++) {
    ref_int32[i] = ref_coeff_pvq[i];
    in_int32[i] = coeff_pvq[i];
  }

  if (abs(in_int32[0] - ref_int32[0]) < pvq_dc_quant * 141/256) { /* 0.55 */
    out_int32[0] = 0;
  }
  else {
    out_int32[0] = OD_DIV_R0(in_int32[0] - ref_int32[0], pvq_dc_quant);
  }

  skip = od_pvq_encode(daala_enc, ref_int32, in_int32, out_int32,
          (int)quant[0] >> quant_shift,//scale/quantizer
          (int)quant[1] >> quant_shift,//scale/quantizer
          plane, tx_size,
          OD_PVQ_BETA[0/*use_activity_masking*/][plane][tx_size],
          1,//OD_ROBUST_STREAM
          0, //is_keyframe,
          0, 0, 0, //q_scaling, bx, by,
          daala_enc->state.qm + off, daala_enc->state.qm_inv + off, pvq_info);

  if (skip)
    assert(pvq_info->ac_dc_coded == 0);

  if (!skip)
    assert(pvq_info->ac_dc_coded > 0);

  // Encode residue of DC coeff, if required.
  if (!has_dc_skip || out_int32[0]) {
    generic_encode(&daala_enc->ec, &daala_enc->state.adapt.model_dc[plane],
     abs(out_int32[0]) - has_dc_skip, -1,
     &daala_enc->state.adapt.ex_dc[plane][tx_size][0], 2);
  }
  if (out_int32[0]) {
    od_ec_enc_bits(&daala_enc->ec, out_int32[0] < 0, 1);
    skip = 0;
  }

  // need to save quantized residue of DC coeff
  // so that final pvq bitstream writing can know whether DC is coded.
  pvq_info->dq_dc_residue = out_int32[0];

  out_int32[0] = out_int32[0] * pvq_dc_quant;
  out_int32[0] += ref_int32[0];

  //copy int32 result back to int16
  for (i=0; i < tx_blk_size * tx_blk_size; i++)
    dqcoeff_pvq[i] = out_int32[i];

  // Safely initialize dqcoeff since some coeffs (band size > 128 coeffs)
  // are skipped by PVQ.
  //od_init_skipped_coeffs(dqcoeff, ref_coeff, 0, 0, tx_blk_size, tx_blk_size);

  // Back to original coefficient order
  od_coding_order_to_raster(dqcoeff, tx_blk_size, dqcoeff_pvq, tx_blk_size);

  // Mark last nonzero coeff position.
  //for (j = 0; j < tx_blk_size * tx_blk_size; j++)
  //  if (dqcoeff[j]) *eob = j + 1;

  *eob = tx_blk_size * tx_blk_size;

  pvq_info->eob = *eob;

  *rate = (od_ec_enc_tell(&daala_enc->ec) - tell) << AV1_PROB_COST_SHIFT;
  assert(*rate >= 0);

  return skip;
}

void store_pvq_enc_info(PVQ_INFO *pvq_info,
                        int *qg,
                        int *theta,
                        int *max_theta,
                        int *k,
                        od_coeff *y,
                        int nb_bands,
                        const int *off,
                        int *size,
                        int skip_rest,
                        int skip_dir,
                        int bs) {  // block size in log_2 -2
  int i;

  for (i=0; i < nb_bands; i++) {
    pvq_info->qg[i] = qg[i];
    pvq_info->theta[i] = theta[i];
    pvq_info->max_theta[i] = max_theta[i];
    pvq_info->k[i] = k[i];
    pvq_info->off[i] = off[i];
    pvq_info->size[i] = size[i];
  }
  // TODO: just copying block size should be fine
  for (i=0; i < OD_BSIZE_MAX*OD_BSIZE_MAX; i++) {
    pvq_info->y[i] = y[i];
  }
  pvq_info->nb_bands = nb_bands;
  pvq_info->skip_rest = skip_rest;
  pvq_info->skip_dir = skip_dir;
  pvq_info->bs = bs;
}
#endif
