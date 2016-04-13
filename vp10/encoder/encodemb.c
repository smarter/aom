/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include "./vp10_rtcd.h"
#include "./vpx_config.h"
#include "./vpx_dsp_rtcd.h"

#include "vpx_dsp/quantize.h"
#include "vpx_mem/vpx_mem.h"
#include "vpx_ports/mem.h"

#include "vp10/common/idct.h"
#include "vp10/common/reconinter.h"
#include "vp10/common/reconintra.h"
#include "vp10/common/scan.h"
#include "vp10/common/partition.h"

#include "vp10/encoder/encodemb.h"
#include "vp10/encoder/rd.h"
#include "vp10/encoder/tokenize.h"

#if CONFIG_PVQ
#include "vp10/encoder/encint.h"
#include "vp10/encoder/pvq_encoder.h"

extern daala_enc_ctx daala_enc;
#endif

struct optimize_ctx {
  ENTROPY_CONTEXT ta[MAX_MB_PLANE][16];
  ENTROPY_CONTEXT tl[MAX_MB_PLANE][16];
};

void vp10_subtract_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane) {
  struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &x->e_mbd.plane[plane];
  const BLOCK_SIZE plane_bsize = get_plane_block_size(bsize, pd);
  const int bw = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  const int bh = 4 * num_4x4_blocks_high_lookup[plane_bsize];

#if CONFIG_VPX_HIGHBITDEPTH
  if (x->e_mbd.cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    vpx_highbd_subtract_block(bh, bw, p->src_diff, bw, p->src.buf,
                              p->src.stride, pd->dst.buf, pd->dst.stride,
                              x->e_mbd.bd);
    return;
  }
#endif  // CONFIG_VPX_HIGHBITDEPTH
  vpx_subtract_block(bh, bw, p->src_diff, bw, p->src.buf, p->src.stride,
                     pd->dst.buf, pd->dst.stride);
}

#define RDTRUNC(RM, DM, R, D)                        \
  (((1 << (VP9_PROB_COST_SHIFT - 1)) + (R) * (RM)) & \
   ((1 << VP9_PROB_COST_SHIFT) - 1))

typedef struct vp10_token_state {
  int rate;
  int error;
  int next;
  int16_t token;
  short qc;
} vp10_token_state;

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

// This function is a place holder for now but may ultimately need
// to scan previous tokens to work out the correct context.
static int trellis_get_coeff_context(const int16_t *scan, const int16_t *nb,
                                     int idx, int token, uint8_t *token_cache) {
  int bak = token_cache[scan[idx]], pt;
  token_cache[scan[idx]] = vp10_pt_energy_class[token];
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
  vp10_token_state tokens[1025][2];
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
#if CONFIG_VPX_HIGHBITDEPTH
  const int *cat6_high_cost = vp10_get_high_cost_table(xd->bd);
#else
  const int *cat6_high_cost = vp10_get_high_cost_table(8);
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
    token_cache[scan[i]] =
        vp10_pt_energy_class[vp10_get_token(qcoeff[scan[i]])];

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
      vp10_get_token_extra(x, &t0, &e0);
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
      base_bits = vp10_get_cost(t0, e0, cat6_high_cost);
      dx = mul * (dqcoeff[rc] - coeff[rc]);
#if CONFIG_VPX_HIGHBITDEPTH
      if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
        dx >>= xd->bd - 8;
      }
#endif  // CONFIG_VPX_HIGHBITDEPTH
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
        vp10_get_token_extra(x, &t0, &e0);
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
      base_bits = vp10_get_cost(t0, e0, cat6_high_cost);

      if (shortcut) {
#if CONFIG_VPX_HIGHBITDEPTH
        if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
          dx -= ((dequant_ptr[rc != 0] >> (xd->bd - 8)) + sz) ^ sz;
        } else {
          dx -= (dequant_ptr[rc != 0] + sz) ^ sz;
        }
#else
        dx -= (dequant_ptr[rc != 0] + sz) ^ sz;
#endif  // CONFIG_VPX_HIGHBITDEPTH
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

static INLINE void fdct32x32(int rd_transform, const int16_t *src,
                             tran_low_t *dst, int src_stride) {
  if (rd_transform)
    vpx_fdct32x32_rd(src, dst, src_stride);
  else
    vpx_fdct32x32(src, dst, src_stride);
}

#if CONFIG_VPX_HIGHBITDEPTH
static INLINE void highbd_fdct32x32(int rd_transform, const int16_t *src,
                                    tran_low_t *dst, int src_stride) {
  if (rd_transform)
    vpx_highbd_fdct32x32_rd(src, dst, src_stride);
  else
    vpx_highbd_fdct32x32(src, dst, src_stride);
}
#endif  // CONFIG_VPX_HIGHBITDEPTH

void vp10_fwd_txfm_4x4(const int16_t *src_diff, tran_low_t *coeff,
                       int diff_stride, TX_TYPE tx_type, int lossless) {
  if (lossless) {
    vp10_fwht4x4(src_diff, coeff, diff_stride);
  } else {
    switch (tx_type) {
      case DCT_DCT:
        vpx_fdct4x4(src_diff, coeff, diff_stride);
        break;
      case ADST_DCT:
      case DCT_ADST:
      case ADST_ADST:
        vp10_fht4x4(src_diff, coeff, diff_stride, tx_type);
        break;
      default:
        assert(0);
        break;
    }
  }
}

void fwd_txfm_8x8(const int16_t *src_diff, tran_low_t *coeff,
                         int diff_stride, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_fht8x8(src_diff, coeff, diff_stride, tx_type);
      break;
    default:
      assert(0);
      break;
  }
}

void fwd_txfm_16x16(const int16_t *src_diff, tran_low_t *coeff,
                           int diff_stride, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_fht16x16(src_diff, coeff, diff_stride, tx_type);
      break;
    default:
      assert(0);
      break;
  }
}

void fwd_txfm_32x32(int rd_transform, const int16_t *src_diff,
                           tran_low_t *coeff, int diff_stride,
                           TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
      fdct32x32(rd_transform, src_diff, coeff, diff_stride);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      assert(0);
      break;
    default:
      assert(0);
      break;
  }
}

#if CONFIG_VPX_HIGHBITDEPTH
void vp10_highbd_fwd_txfm_4x4(const int16_t *src_diff, tran_low_t *coeff,
                              int diff_stride, TX_TYPE tx_type, int lossless) {
  if (lossless) {
    assert(tx_type == DCT_DCT);
    vp10_highbd_fwht4x4(src_diff, coeff, diff_stride);
  } else {
    switch (tx_type) {
      case DCT_DCT:
        vpx_highbd_fdct4x4(src_diff, coeff, diff_stride);
        break;
      case ADST_DCT:
      case DCT_ADST:
      case ADST_ADST:
        vp10_highbd_fht4x4(src_diff, coeff, diff_stride, tx_type);
        break;
      default:
        assert(0);
        break;
    }
  }
}

static void highbd_fwd_txfm_8x8(const int16_t *src_diff, tran_low_t *coeff,
                                int diff_stride, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
      vpx_highbd_fdct8x8(src_diff, coeff, diff_stride);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_highbd_fht8x8(src_diff, coeff, diff_stride, tx_type);
      break;
    default:
      assert(0);
      break;
  }
}

static void highbd_fwd_txfm_16x16(const int16_t *src_diff, tran_low_t *coeff,
                                  int diff_stride, TX_TYPE tx_type) {
  switch (tx_type) {
    case DCT_DCT:
      vpx_highbd_fdct16x16(src_diff, coeff, diff_stride);
      break;
    case ADST_DCT:
    case DCT_ADST:
    case ADST_ADST:
      vp10_highbd_fht16x16(src_diff, coeff, diff_stride, tx_type);
      break;
    default:
      assert(0);
      break;
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
    case ADST_ADST:
      assert(0);
      break;
    default:
      assert(0);
      break;
  }
}
#endif  // CONFIG_VPX_HIGHBITDEPTH

void vp10_xform_quant_fp(MACROBLOCK *x, int plane, int block, int blk_row,
                         int blk_col, BLOCK_SIZE plane_bsize, TX_SIZE tx_size) {
  MACROBLOCKD *const xd = &x->e_mbd;
#if !CONFIG_PVQ
  const struct macroblock_plane *const p = &x->plane[plane];
  const struct macroblockd_plane *const pd = &xd->plane[plane];
#else
  struct macroblock_plane *const p = &x->plane[plane];
  struct macroblockd_plane *const pd = &xd->plane[plane];
#endif
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  PLANE_TYPE plane_type = (plane == 0) ? PLANE_TYPE_Y : PLANE_TYPE_UV;
  TX_TYPE tx_type = get_tx_type(plane_type, xd, block);
  const scan_order *const scan_order = get_scan(tx_size, tx_type);
  tran_low_t *const coeff = BLOCK_OFFSET(p->coeff, block);
  tran_low_t *const qcoeff = BLOCK_OFFSET(p->qcoeff, block);
  tran_low_t *const dqcoeff = BLOCK_OFFSET(pd->dqcoeff, block);
  uint16_t *const eob = &p->eobs[block];
  const int diff_stride = 4 * num_4x4_blocks_wide_lookup[plane_bsize];
  int seg_id = mbmi->segment_id;
#if CONFIG_AOM_QM
  int is_intra = !is_inter_block(&xd->mi[0]->mbmi);
  const qm_val_t *qmatrix = pd->seg_qmatrix[seg_id][is_intra][tx_size];
  const qm_val_t *iqmatrix = pd->seg_iqmatrix[seg_id][is_intra][tx_size];
#endif

#if !CONFIG_PVQ
  const int16_t *src_diff;
  src_diff = &p->src_diff[4 * (blk_row * diff_stride + blk_col)];
#else
  uint8_t *src, *dst;
  tran_low_t *ref_coeff = BLOCK_OFFSET(pd->pvq_ref_coeff, block);
  int16_t *src_int16, *pred;
  const int src_stride = p->src.stride;
  const int dst_stride = pd->dst.stride;
  int tx_blk_size;
  int i, j;
  int skip;

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

#if CONFIG_VPX_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
      case TX_32X32:
        highbd_fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
        vp10_highbd_quantize_fp_32x32(coeff, 1024, x->skip_block, p->zbin,
                                      p->round_fp, p->quant_fp, p->quant_shift,
                                      qcoeff, dqcoeff, pd->dequant, eob,
                                      scan_order->scan,
#if !CONFIG_AOM_QM
                                      scan_order->iscan);
#else
                                      scan_order->iscan, qmatrix, iqmatrix);
#endif
        break;
      case TX_16X16:
        vpx_highbd_fdct16x16(src_diff, coeff, diff_stride);
        vp10_highbd_quantize_fp(coeff, 256, x->skip_block, p->zbin, p->round_fp,
                                p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                                pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                                scan_order->iscan);
#else
                                scan_order->iscan, qmatrix, iqmatrix);
#endif
        break;
      case TX_8X8:
        vpx_highbd_fdct8x8(src_diff, coeff, diff_stride);
        vp10_highbd_quantize_fp(coeff, 64, x->skip_block, p->zbin, p->round_fp,
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
          vp10_highbd_fwht4x4(src_diff, coeff, diff_stride);
        } else {
          vpx_highbd_fdct4x4(src_diff, coeff, diff_stride);
        }
        vp10_highbd_quantize_fp(coeff, 16, x->skip_block, p->zbin, p->round_fp,
                                p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                                pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                                scan_order->iscan);
#else
                                scan_order->iscan, qmatrix, iqmatrix);
#endif
        break;
      default:
        assert(0);
    }
    return;
  }
#endif  // CONFIG_VPX_HIGHBITDEPTH

#if !CONFIG_PVQ
  switch (tx_size) {
    case TX_32X32:
      fdct32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride);
      vp10_quantize_fp_32x32(coeff, 1024, x->skip_block, p->zbin, p->round_fp,
                             p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                             pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                             scan_order->iscan);
#else
                             scan_order->iscan, qmatrix, iqmatrix);
#endif
      break;
    case TX_16X16:
      vpx_fdct16x16(src_diff, coeff, diff_stride);
      vp10_quantize_fp(coeff, 256, x->skip_block, p->zbin, p->round_fp,
                       p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                       pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                       scan_order->iscan);
#else
                       scan_order->iscan, qmatrix, iqmatrix);
#endif
      break;
    case TX_8X8:
      /*vp10_fdct8x8_quant(src_diff, diff_stride, coeff, 64, x->skip_block,
                         p->zbin, p->round_fp, p->quant_fp, p->quant_shift,
                         qcoeff, dqcoeff, pd->dequant, eob, scan_order->scan,
                         scan_order->iscan);*/
      vpx_fdct8x8(src_diff, coeff, diff_stride);
      vp10_quantize_fp(coeff, 64, x->skip_block, p->zbin, p->round_fp,
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
        vp10_fwht4x4(src_diff, coeff, diff_stride);
      } else {
        vpx_fdct4x4(src_diff, coeff, diff_stride);
      }
      vp10_quantize_fp(coeff, 16, x->skip_block, p->zbin, p->round_fp,
                       p->quant_fp, p->quant_shift, qcoeff, dqcoeff,
                       pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                       scan_order->iscan);
#else
                       scan_order->iscan, qmatrix, iqmatrix);
#endif
      break;
    default:
      assert(0);
      break;
  }
#else//#if !CONFIG_PVQ
  switch (tx_size) {
    case TX_32X32:
      //forward transform of predicted image.
      fdct32x32(x->use_lp32x32fdct, pred, ref_coeff, diff_stride);
      //forward transform of original image.
      fdct32x32(x->use_lp32x32fdct, src_int16, coeff, diff_stride);
      break;
    case TX_16X16:
      vpx_fdct16x16(pred, ref_coeff, diff_stride);
      vpx_fdct16x16(src_int16, coeff, diff_stride);
      break;
    case TX_8X8:
      vpx_fdct8x8(pred, ref_coeff, diff_stride);
      vpx_fdct8x8(src_int16, coeff, diff_stride);
      break;
    case TX_4X4:
      if (xd->lossless[seg_id]) {
        vp10_fwht4x4(pred, ref_coeff, diff_stride);
        vp10_fwht4x4(src_int16, coeff, diff_stride);
      } else {
        vpx_fdct4x4(pred, ref_coeff, diff_stride);
        vpx_fdct4x4(src_int16, coeff, diff_stride);
      }
      break;
    default: assert(0); break;
  }

  // pvq of daala will be called here for inter mode block
  skip = pvq_encode_helper2(coeff,          // target original vector
                            ref_coeff,      // reference vector
                            dqcoeff,        // de-quantized vector
                            eob,             // End of Block marker
                            pd->dequant[0], // vpx's DC quantization step size
                            0,              // keyframe (daala's definition)? 0 for now
                            tx_size,        // block size in log_2 - 2, 0 for 4x4.
                            &x->rate,       // rate measured
                            &mbmi->pvq[plane]); // PVQ info for a block

  if (!skip)
    mbmi->skip = 0;
#endif//#if !CONFIG_PVQ
}

void vp10_xform_quant_dc(MACROBLOCK *x, int plane, int block, int blk_row,
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

#if CONFIG_VPX_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
      case TX_32X32:
        vpx_highbd_fdct32x32_1(src_diff, coeff, diff_stride);
        vpx_highbd_quantize_dc_32x32(coeff, x->skip_block, p->round,
                                     p->quant_fp[0], qcoeff, dqcoeff,
                                     pd->dequant[0],
#if !CONFIG_AOM_QM
                                     eob);
#else
                                     eob, qmatrix, iqmatrix);
#endif
        break;
      case TX_16X16:
        vpx_highbd_fdct16x16_1(src_diff, coeff, diff_stride);
        vpx_highbd_quantize_dc(coeff, 256, x->skip_block, p->round,
                               p->quant_fp[0], qcoeff, dqcoeff, pd->dequant[0],
#if !CONFIG_AOM_QM
                               eob);
#else
                               eob, qmatrix, iqmatrix);
#endif
        break;
      case TX_8X8:
        vpx_highbd_fdct8x8_1(src_diff, coeff, diff_stride);
        vpx_highbd_quantize_dc(coeff, 64, x->skip_block, p->round,
                               p->quant_fp[0], qcoeff, dqcoeff, pd->dequant[0],
#if !CONFIG_AOM_QM
                               eob);
#else
                               eob, qmatrix, iqmatrix);
#endif
        break;
      case TX_4X4:
        if (xd->lossless[seg_id]) {
          vp10_highbd_fwht4x4(src_diff, coeff, diff_stride);
        } else {
          vpx_highbd_fdct4x4(src_diff, coeff, diff_stride);
        }
        vpx_highbd_quantize_dc(coeff, 16, x->skip_block, p->round,
                               p->quant_fp[0], qcoeff, dqcoeff, pd->dequant[0],
#if !CONFIG_AOM_QM
                               eob);
#else
                               eob, qmatrix, iqmatrix);
#endif
        break;
      default:
        assert(0);
    }
    return;
  }
#endif  // CONFIG_VPX_HIGHBITDEPTH

  switch (tx_size) {
    case TX_32X32:
      vpx_fdct32x32_1(src_diff, coeff, diff_stride);
      vpx_quantize_dc_32x32(coeff, x->skip_block, p->round, p->quant_fp[0],
                            qcoeff, dqcoeff, pd->dequant[0],
#if !CONFIG_AOM_QM
                            eob);
#else
                            eob, qmatrix, iqmatrix);
#endif
      break;
    case TX_16X16:
      vpx_fdct16x16_1(src_diff, coeff, diff_stride);
      vpx_quantize_dc(coeff, 256, x->skip_block, p->round, p->quant_fp[0],
                      qcoeff, dqcoeff, pd->dequant[0],
#if !CONFIG_AOM_QM
                      eob);
#else
                      eob, qmatrix, iqmatrix);
#endif
      break;
    case TX_8X8:
      vpx_fdct8x8_1(src_diff, coeff, diff_stride);
      vpx_quantize_dc(coeff, 64, x->skip_block, p->round, p->quant_fp[0],
                      qcoeff, dqcoeff, pd->dequant[0],
#if !CONFIG_AOM_QM
                      eob);
#else
                      eob, qmatrix, iqmatrix);
#endif
      break;
    case TX_4X4:
      if (xd->lossless[seg_id]) {
        vp10_fwht4x4(src_diff, coeff, diff_stride);
      } else {
        vpx_fdct4x4(src_diff, coeff, diff_stride);
      }
      vpx_quantize_dc(coeff, 16, x->skip_block, p->round, p->quant_fp[0],
                      qcoeff, dqcoeff, pd->dequant[0],
#if !CONFIG_AOM_QM
                      eob);
#else
                      eob, qmatrix, iqmatrix);
#endif
      break;
    default:
      assert(0);
      break;
  }
}

void vp10_xform_quant(MACROBLOCK *x, int plane, int block, int blk_row,
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
  int16_t *pred = &pd->pred[4 * (blk_row * diff_stride + blk_col)];
  uint8_t *src, *dst;
  int16_t *src_int16;
  const int src_stride = p->src.stride;
  const int dst_stride = pd->dst.stride;
  int tx_blk_size;
  int i, j;
  int skip;
  DECLARE_ALIGNED(16, int16_t, coeff_pvq[64 * 64]);
  DECLARE_ALIGNED(16, int16_t, dqcoeff_pvq[64 * 64]);
  DECLARE_ALIGNED(16, int16_t, ref_coeff_pvq[64 * 64]);
  dst = &pd->dst.buf[4 * (blk_row * dst_stride + blk_col)];
  src = &p->src.buf[4 * (blk_row * src_stride + blk_col)];
  src_int16 = &p->src_int16[4 * (blk_row * diff_stride + blk_col)];
#endif

#if CONFIG_PVQ
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

#if CONFIG_VPX_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
      case TX_32X32:
        highbd_fwd_txfm_32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride,
                              tx_type);
        vpx_highbd_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin,
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
        vpx_highbd_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round,
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
        vpx_highbd_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                              scan_order->iscan);
#else
                              scan_order->iscan, qmatrix, iqmatrix);
#endif
        break;
      case TX_4X4:
        vp10_highbd_fwd_txfm_4x4(src_diff, coeff, diff_stride, tx_type,
                                 xd->lossless[seg_id]);
        vpx_highbd_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round,
                              p->quant, p->quant_shift, qcoeff, dqcoeff,
                              pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                              scan_order->iscan);
#else
                              scan_order->iscan, qmatrix, iqmatrix);
#endif
        break;
      default:
        assert(0);
    }
    return;
  }
#endif  // CONFIG_VPX_HIGHBITDEPTH

#if !CONFIG_PVQ
  switch (tx_size) {
    case TX_32X32:
      fwd_txfm_32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride, tx_type);
      vpx_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin, p->round,
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
      vpx_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round, p->quant,
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
      vpx_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round, p->quant,
                     p->quant_shift, qcoeff, dqcoeff, pd->dequant, eob,
                     scan_order->scan,
#if !CONFIG_AOM_QM
                     scan_order->iscan);
#else
                     scan_order->iscan, qmatrix, iqmatrix);
#endif
      break;
    case TX_4X4:
      vp10_fwd_txfm_4x4(src_diff, coeff, diff_stride, tx_type,
                        xd->lossless[seg_id]);
      vpx_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round, p->quant,
                     p->quant_shift, qcoeff, dqcoeff, pd->dequant, eob,
                     scan_order->scan,
#if !CONFIG_AOM_QM
                     scan_order->iscan);
#else
                     scan_order->iscan, qmatrix, iqmatrix);
#endif
      break;
    default:
      assert(0);
      break;
  }
#else//#if !CONFIG_PVQ
  switch (tx_size) {
    case TX_32X32:
      //forward transform of predicted image.
      fwd_txfm_32x32(x->use_lp32x32fdct, pred, ref_coeff, diff_stride, tx_type);
      //forward transform of original image.
      fwd_txfm_32x32(x->use_lp32x32fdct, src_int16, coeff, diff_stride, tx_type);
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
      vp10_fwd_txfm_4x4(pred, ref_coeff, diff_stride, tx_type,
                        xd->lossless[seg_id]);
      vp10_fwd_txfm_4x4(src_int16, coeff, diff_stride, tx_type,
                        xd->lossless[seg_id]);
      break;
    default: assert(0); break;
  }

  // pvq of daala will be called here for inter mode block
  if (!x->skip_block)
  skip = pvq_encode_helper2(coeff,          // target original vector
                            ref_coeff,      // reference vector
                            dqcoeff,        // de-quantized vector
                            eob,            // End of Block marker
                            pd->dequant[0], // vpx's DC quantization step size
                            plane,          // image plane
                            tx_size,        // block size in log_2 - 2, 0 for 4x4.
                            &x->rate,       // rate measured
                            &mbmi->pvq[plane]); // PVQ info for a block

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
  int seg_id = xd->mi[0]->mbmi.segment_id;
#if CONFIG_PVQ
  int tx_blk_size;
  int i, j;
  tran_low_t *ref_coeff = BLOCK_OFFSET(pd->pvq_ref_coeff, block);
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  //block skip info from pvq, 0 means both DC and AC are skipped.
  int skip;
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

  if (!x->skip_recode) {
    if (x->quant_fp) {
      // Encoding process for rtc mode
      if (x->skip_txfm[0] == SKIP_TXFM_AC_DC && plane == 0) {
        // skip forward transform
        p->eobs[block] = 0;
        *a = *l = 0;
        return;
      } else {
        vp10_xform_quant_fp(x, plane, block, blk_row, blk_col, plane_bsize,
                            tx_size);
      }
    } else {
      if (max_txsize_lookup[plane_bsize] == tx_size) {
        int txfm_blk_index = (plane << 2) + (block >> (tx_size << 1));
        if (x->skip_txfm[txfm_blk_index] == SKIP_TXFM_NONE) {
          // full forward transform and quantization
          vp10_xform_quant(x, plane, block, blk_row, blk_col, plane_bsize,
                           tx_size);
        } else if (x->skip_txfm[txfm_blk_index] == SKIP_TXFM_AC_ONLY) {
          // fast path forward transform and quantization
          vp10_xform_quant_dc(x, plane, block, blk_row, blk_col, plane_bsize,
                              tx_size);
        } else {
          // skip forward transform
          p->eobs[block] = 0;
          *a = *l = 0;
          return;
        }
      } else {
        vp10_xform_quant(x, plane, block, blk_row, blk_col, plane_bsize,
                         tx_size);
      }
    }
  }

  if (x->optimize && (!x->skip_recode || !x->skip_optimize)) {
    const int ctx = combine_entropy_contexts(*a, *l);
    *a = *l = optimize_b(x, plane, block, tx_size, ctx) > 0;
  } else {
    *a = *l = p->eobs[block] > 0;
  }

#if !CONFIG_PVQ
  if (p->eobs[block]) *(args->skip) = 0;
#else
  //if (p->eobs[block]) *(args->skip) = 0;
  if (mbmi->pvq[plane].ac_dc_coded)
    *(args->skip) = 0;
#endif

  if (p->eobs[block] == 0) return;

#if CONFIG_PVQ
  // transform block size in pixels
  tx_blk_size = 1 << (tx_size + 2);

  // Since vp10 does not have inverse transform only function
  // but contain adding the inverse transform to predicted image,
  // pass blank dummy image to vp10_inv_txfm_add_*x*(), i.e. set dst as zeros
  if (mbmi->pvq[plane].ac_dc_coded)
  for (j=0; j < tx_blk_size; j++)
    memset(dst + j * pd->dst.stride, 0, tx_blk_size);
#endif

#if CONFIG_VPX_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
      case TX_32X32:
        vp10_highbd_inv_txfm_add_32x32(dqcoeff, dst, pd->dst.stride,
                                       p->eobs[block], xd->bd, tx_type);
        break;
      case TX_16X16:
        vp10_highbd_inv_txfm_add_16x16(dqcoeff, dst, pd->dst.stride,
                                       p->eobs[block], xd->bd, tx_type);
        break;
      case TX_8X8:
        vp10_highbd_inv_txfm_add_8x8(dqcoeff, dst, pd->dst.stride,
                                     p->eobs[block], xd->bd, tx_type);
        break;
      case TX_4X4:
        // this is like vp10_short_idct4x4 but has a special case around eob<=1
        // which is significant (not just an optimization) for the lossless
        // case.
        vp10_highbd_inv_txfm_add_4x4(dqcoeff, dst, pd->dst.stride,
                                     p->eobs[block], xd->bd, tx_type,
                                     xd->lossless[seg_id]);
        break;
      default:
        assert(0 && "Invalid transform size");
        break;
    }

    return;
  }
#endif  // CONFIG_VPX_HIGHBITDEPTH

#if CONFIG_PVQ
  if (mbmi->pvq[plane].ac_dc_coded)
#endif
  switch (tx_size) {
    case TX_32X32:
      vp10_inv_txfm_add_32x32(dqcoeff, dst, pd->dst.stride, p->eobs[block],
                              tx_type);
      break;
    case TX_16X16:
      vp10_inv_txfm_add_16x16(dqcoeff, dst, pd->dst.stride, p->eobs[block],
                              tx_type);
      break;
    case TX_8X8:
      vp10_inv_txfm_add_8x8(dqcoeff, dst, pd->dst.stride, p->eobs[block],
                            tx_type);
      break;
    case TX_4X4:
      // this is like vp10_short_idct4x4 but has a special case around eob<=1
      // which is significant (not just an optimization) for the lossless
      // case.
      vp10_inv_txfm_add_4x4(dqcoeff, dst, pd->dst.stride, p->eobs[block],
                            tx_type, xd->lossless[seg_id]);
      break;
    default:
      assert(0 && "Invalid transform size");
      break;
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

  vp10_xform_quant(x, plane, block, blk_row, blk_col, plane_bsize, tx_size);

  if (p->eobs[block] > 0) {
#if CONFIG_VPX_HIGHBITDEPTH
    if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
      if (xd->lossless[0]) {
        vp10_highbd_iwht4x4_add(dqcoeff, dst, pd->dst.stride, p->eobs[block],
                                xd->bd);
      } else {
        vp10_highbd_idct4x4_add(dqcoeff, dst, pd->dst.stride, p->eobs[block],
                                xd->bd);
      }
      return;
    }
#endif  // CONFIG_VPX_HIGHBITDEPTH
    if (xd->lossless[0]) {
      vp10_iwht4x4_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
    } else {
      vp10_idct4x4_add(dqcoeff, dst, pd->dst.stride, p->eobs[block]);
    }
  }
}

void vp10_encode_sby_pass1(MACROBLOCK *x, BLOCK_SIZE bsize) {
  vp10_subtract_plane(x, bsize, 0);
  vp10_foreach_transformed_block_in_plane(&x->e_mbd, bsize, 0,
                                          encode_block_pass1, x);
}

void vp10_encode_sb(MACROBLOCK *x, BLOCK_SIZE bsize) {
  MACROBLOCKD *const xd = &x->e_mbd;
  struct optimize_ctx ctx;
  MB_MODE_INFO *mbmi = &xd->mi[0]->mbmi;
  struct encode_b_args arg = { x, &ctx, &mbmi->skip };
  int plane;

  mbmi->skip = 1;

  if (x->skip) return;

  for (plane = 0; plane < MAX_MB_PLANE; ++plane) {
    if (!x->skip_recode) vp10_subtract_plane(x, bsize, plane);

    if (x->optimize && (!x->skip_recode || !x->skip_optimize)) {
      const struct macroblockd_plane *const pd = &xd->plane[plane];
      const TX_SIZE tx_size = plane ? get_uv_tx_size(mbmi, pd) : mbmi->tx_size;
      vp10_get_entropy_contexts(bsize, tx_size, pd, ctx.ta[plane],
                                ctx.tl[plane]);
    }

    vp10_foreach_transformed_block_in_plane(xd, bsize, plane, encode_block,
                                            &arg);
  }
}

void vp10_encode_block_intra(int plane, int block, int blk_row, int blk_col,
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
  int skip;
  DECLARE_ALIGNED(16, int16_t, coeff_pvq[64 * 64]);
  DECLARE_ALIGNED(16, int16_t, dqcoeff_pvq[64 * 64]);
  DECLARE_ALIGNED(16, int16_t, ref_coeff_pvq[64 * 64]);
  src_int16 = &p->src_int16[4 * (blk_row * diff_stride + blk_col)];
#endif
  dst = &pd->dst.buf[4 * (blk_row * dst_stride + blk_col)];
  src = &p->src.buf[4 * (blk_row * src_stride + blk_col)];
  mode = plane == 0 ? get_y_mode(xd->mi[0], block) : mbmi->uv_mode;
  vp10_predict_intra_block(xd, bwl, bhl, tx_size, mode, dst, dst_stride, dst,
                           dst_stride, blk_col, blk_row, plane);

#if CONFIG_VPX_HIGHBITDEPTH
  if (xd->cur_buf->flags & YV12_FLAG_HIGHBITDEPTH) {
    switch (tx_size) {
      case TX_32X32:
        if (!x->skip_recode) {
          vpx_highbd_subtract_block(32, 32, src_diff, diff_stride, src,
                                    src_stride, dst, dst_stride, xd->bd);
          highbd_fwd_txfm_32x32(x->use_lp32x32fdct, src_diff, coeff,
                                diff_stride, tx_type);
          vpx_highbd_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin,
                                      p->round, p->quant, p->quant_shift,
                                      qcoeff, dqcoeff, pd->dequant, eob,
                                      scan_order->scan,
#if !CONFIG_AOM_QM
                                      scan_order->iscan);
#else
                                      scan_order->iscan, qmatrix, iqmatrix);
#endif
        }
        if (*eob)
          vp10_highbd_inv_txfm_add_32x32(dqcoeff, dst, dst_stride, *eob, xd->bd,
                                         tx_type);
        break;
      case TX_16X16:
        if (!x->skip_recode) {
          vpx_highbd_subtract_block(16, 16, src_diff, diff_stride, src,
                                    src_stride, dst, dst_stride, xd->bd);
          highbd_fwd_txfm_16x16(src_diff, coeff, diff_stride, tx_type);
          vpx_highbd_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round,
                                p->quant, p->quant_shift, qcoeff, dqcoeff,
                                pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                                scan_order->iscan);
#else
                                scan_order->iscan, qmatrix, iqmatrix);
#endif
        }
        if (*eob)
          vp10_highbd_inv_txfm_add_16x16(dqcoeff, dst, dst_stride, *eob, xd->bd,
                                         tx_type);
        break;
      case TX_8X8:
        if (!x->skip_recode) {
          vpx_highbd_subtract_block(8, 8, src_diff, diff_stride, src,
                                    src_stride, dst, dst_stride, xd->bd);
          highbd_fwd_txfm_8x8(src_diff, coeff, diff_stride, tx_type);
          vpx_highbd_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round,
                                p->quant, p->quant_shift, qcoeff, dqcoeff,
                                pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                                scan_order->iscan);
#else
                                scan_order->iscan, qmatrix, iqmatrix);
#endif
        }
        if (*eob)
          vp10_highbd_inv_txfm_add_8x8(dqcoeff, dst, dst_stride, *eob, xd->bd,
                                       tx_type);
        break;
      case TX_4X4:
        if (!x->skip_recode) {
          vpx_highbd_subtract_block(4, 4, src_diff, diff_stride, src,
                                    src_stride, dst, dst_stride, xd->bd);
          vp10_highbd_fwd_txfm_4x4(src_diff, coeff, diff_stride, tx_type,
                                   xd->lossless[seg_id]);
          vpx_highbd_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round,
                                p->quant, p->quant_shift, qcoeff, dqcoeff,
                                pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                                scan_order->iscan);
#else
                                scan_order->iscan, qmatrix, iqmatrix);
#endif
        }

        if (*eob)
          // this is like vp10_short_idct4x4 but has a special case around
          // eob<=1 which is significant (not just an optimization) for the
          // lossless case.
          vp10_highbd_inv_txfm_add_4x4(dqcoeff, dst, dst_stride, *eob, xd->bd,
                                       tx_type, xd->lossless[seg_id]);
        break;
      default:
        assert(0);
        return;
    }
    if (*eob) *(args->skip) = 0;
    return;
  }
#endif  // CONFIG_VPX_HIGHBITDEPTH

#if !CONFIG_PVQ
  switch (tx_size) {
    case TX_32X32:
      if (!x->skip_recode) {
        vpx_subtract_block(32, 32, src_diff, diff_stride, src, src_stride, dst,
                           dst_stride);
        fwd_txfm_32x32(x->use_lp32x32fdct, src_diff, coeff, diff_stride,
                       tx_type);
        vpx_quantize_b_32x32(coeff, 1024, x->skip_block, p->zbin, p->round,
                             p->quant, p->quant_shift, qcoeff, dqcoeff,
                             pd->dequant, eob, scan_order->scan,
#if !CONFIG_AOM_QM
                             scan_order->iscan);
#else
                             scan_order->iscan, qmatrix, iqmatrix);
#endif
      }
      if (*eob)
        vp10_inv_txfm_add_32x32(dqcoeff, dst, dst_stride, *eob, tx_type);
      break;
    case TX_16X16:
      if (!x->skip_recode) {
        vpx_subtract_block(16, 16, src_diff, diff_stride, src, src_stride, dst,
                           dst_stride);
        fwd_txfm_16x16(src_diff, coeff, diff_stride, tx_type);
        vpx_quantize_b(coeff, 256, x->skip_block, p->zbin, p->round, p->quant,
                       p->quant_shift, qcoeff, dqcoeff, pd->dequant, eob,
                       scan_order->scan,
#if !CONFIG_AOM_QM
                       scan_order->iscan);
#else
                       scan_order->iscan, qmatrix, iqmatrix);
#endif
      }
      if (*eob)
        vp10_inv_txfm_add_16x16(dqcoeff, dst, dst_stride, *eob, tx_type);
      break;
    case TX_8X8:
      if (!x->skip_recode) {
        vpx_subtract_block(8, 8, src_diff, diff_stride, src, src_stride, dst,
                           dst_stride);
        fwd_txfm_8x8(src_diff, coeff, diff_stride, tx_type);
        vpx_quantize_b(coeff, 64, x->skip_block, p->zbin, p->round, p->quant,
                       p->quant_shift, qcoeff, dqcoeff, pd->dequant, eob,
                       scan_order->scan,
#if !CONFIG_AOM_QM
                       scan_order->iscan);
#else
                       scan_order->iscan, qmatrix, iqmatrix);
#endif
      }
      if (*eob) vp10_inv_txfm_add_8x8(dqcoeff, dst, dst_stride, *eob, tx_type);
      break;
    case TX_4X4:
      if (!x->skip_recode) {
        vpx_subtract_block(4, 4, src_diff, diff_stride, src, src_stride, dst,
                           dst_stride);
        vp10_fwd_txfm_4x4(src_diff, coeff, diff_stride, tx_type,
                          xd->lossless[seg_id]);
        vpx_quantize_b(coeff, 16, x->skip_block, p->zbin, p->round, p->quant,
                       p->quant_shift, qcoeff, dqcoeff, pd->dequant, eob,
                       scan_order->scan,
#if !CONFIG_AOM_QM
                       scan_order->iscan);
#else
                       scan_order->iscan, qmatrix, iqmatrix);
#endif
      }

      if (*eob) {
        // this is like vp10_short_idct4x4 but has a special case around eob<=1
        // which is significant (not just an optimization) for the lossless
        // case.
        vp10_inv_txfm_add_4x4(dqcoeff, dst, dst_stride, *eob, tx_type,
                              xd->lossless[seg_id]);
      }
      break;
    default:
      assert(0);
      break;
  }
#else//#if !CONFIG_PVQ
  // transform block size in pixels
  tx_blk_size = 1 << (tx_size + 2);

  // Instead of computing residue in pixel domain,
  // pvq uses the residue defined in freq domain.
  // For this, forward transform is applied to 1) predicted image
  // and 2) target original image.
  // Note that pvq in decoder also need to apply forward transform
  // to predicted image to obtain reference vector.
  if (!x->skip_recode) {
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
        fwd_txfm_32x32(x->use_lp32x32fdct, pred, ref_coeff, diff_stride,
                       tx_type);
        //forward transform of original image.
        fwd_txfm_32x32(x->use_lp32x32fdct, src_int16, coeff, diff_stride,
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
        vp10_fwd_txfm_4x4(pred, ref_coeff, diff_stride, tx_type,
                          xd->lossless[seg_id]);
        vp10_fwd_txfm_4x4(src_int16, coeff, diff_stride, tx_type,
                          xd->lossless[seg_id]);
        break;
      default: assert(0); break;
    }
    // pvq of daala will be called here for intra mode block
    if (!x->skip_block)
    skip = pvq_encode_helper2(coeff,          // target original vector
                              ref_coeff,      // reference vector
                              dqcoeff,        // de-quantized vector
                              eob,            // End of Block marker
                              pd->dequant[0], // vpx's DC quantization step size
                              plane,          // image plane
                              tx_size,        // block size in log_2 - 2, 0 for 4x4.
                              &x->rate,       // rate measured
                              &mbmi->pvq[plane]); // PVQ info for a block

    if (!skip)
      mbmi->skip = 0;

    // NOTE: *eob == 0 and skip == 0 are not equivalent,
    // since the later means PVQ has not coded both DC and AC,
    // while the former can be false when skip == 0 but the ref vector is nonzero.
    // *eob == 0 is superset to skip == 0.
    // If *eob == 0, decoder does not call pvq_decode().
    // If skip == 0, decoder just use spatial domain prediction.
  }//if (!x->skip_recode) {

  // Since vp10 does not have inverse transform only function
  // but contain adding the inverse transform to predicted image,
  // pass blank dummy image to vp10_inv_txfm_add_*x*(), i.e. set dst as zeros

  if (!skip) {
    for (j=0; j < tx_blk_size; j++)
      memset(dst + j * dst_stride, 0, tx_blk_size);

    switch (tx_size) {
      case TX_32X32:
       vp10_inv_txfm_add_32x32(dqcoeff, dst, dst_stride, *eob, tx_type);
        break;
      case TX_16X16:
        vp10_inv_txfm_add_16x16(dqcoeff, dst, dst_stride, *eob, tx_type);
        break;
      case TX_8X8:
        vp10_inv_txfm_add_8x8(dqcoeff, dst, dst_stride, *eob, tx_type);
        break;
      case TX_4X4:
        // this is like vp10_short_idct4x4 but has a special case around eob<=1
        // which is significant (not just an optimization) for the lossless
        // case.
        vp10_inv_txfm_add_4x4(dqcoeff, dst, dst_stride, *eob, tx_type,
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

void vp10_encode_intra_block_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane) {
  const MACROBLOCKD *const xd = &x->e_mbd;
  struct encode_b_args arg = { x, NULL, &xd->mi[0]->mbmi.skip };

  vp10_foreach_transformed_block_in_plane(xd, bsize, plane,
                                          vp10_encode_block_intra, &arg);
}

#if CONFIG_PVQ
/** Helper function to daala's PVQ encode function
 *
 * @param [in]     ref     'reference' (prediction) vector
 * @param [in]     in      coefficient block to quantize and encode
 * @param [out]    out     quantized coefficient block
 * @param [in]     quant   scale/quantizer
 * @param [in]     pli     plane index
 * @param [in]     bs      log of the block size minus two
 * @param [in]     is_keyframe whether we're encoding a keyframe
 * @return         Returns 1 if both DC and AC coefficients are skipped,
 *                 zero otherwise
 */
int pvq_encode_helper(daala_enc_ctx *daala_enc,
                   int16_t *ref,
                   const int16_t *in,
                   int16_t *out,
                   int quant,
                   int pli,
                   int bs,
                   int is_keyframe,
                   PVQ_INFO *pvq_info) {
  int off = od_qm_offset(bs, pli ? 1 : 0);
  int skip;
  int i;
  od_coeff ref_int32[64*64];
  od_coeff in_int32[64*64];
  od_coeff out_int32[64*64];
  int blk_size = 1 << (bs + 2);

  //copy int16 inputs to int32
  for (i=0; i < blk_size*blk_size; i++) {
    ref_int32[i] = ref[i];
    in_int32[i] = in[i];
    //out_int32[i] = out[i];
  }
  //out_int32[0] = out[0];

  skip = od_pvq_encode(daala_enc, ref_int32, in_int32, out_int32,
          quant,//scale/quantizer
          pli, bs,
          OD_PVQ_BETA[0/*use_activity_masking*/][pli][bs],
          1,//OD_ROBUST_STREAM
          is_keyframe,
          0, 0, 0, //q_scaling, bx, by,
          daala_enc->state.qm + off, daala_enc->state.qm_inv + off, pvq_info);

  if (skip)
    assert(pvq_info->ac_dc_coded == 0);

  if (!skip)
    assert(pvq_info->ac_dc_coded > 0);

  //copy int32 result back to int16
  for (i=0; i < blk_size*blk_size; i++)
    out[i] = out_int32[i];

  return skip;
}

int pvq_encode_helper2(tran_low_t *const coeff, tran_low_t *ref_coeff,
    tran_low_t *const dqcoeff,
    uint16_t *eob, int quant,
    int plane, int tx_size, int *rate, PVQ_INFO *pvq_info) {
  const int tx_blk_size = 1 << (tx_size + 2);
  int skip;
  int j;

  if (plane == 0)
    assert(tx_size > TX_4X4);
  // TODO: Enable this later, if pvq_qm_q4 is available in AOM.
  //int pvq_dc_quant = OD_MAXI(1,
  //  quant * daala_enc.state.pvq_qm_q4[plane][od_qm_get_index(tx_size, 0)] >> 4);
  int pvq_dc_quant = OD_MAXI(1, quant);
  int tell;
  int has_dc_skip = 1;

  DECLARE_ALIGNED(16, int16_t, coeff_pvq[64 * 64]);
  DECLARE_ALIGNED(16, int16_t, dqcoeff_pvq[64 * 64]);
  DECLARE_ALIGNED(16, int16_t, ref_coeff_pvq[64 * 64]);

  *eob = 0;

  // Change coefficient ordering for pvq encoding.
  od_raster_to_coding_order(coeff_pvq, tx_blk_size, coeff, tx_blk_size);
  od_raster_to_coding_order(ref_coeff_pvq, tx_blk_size, ref_coeff, tx_blk_size);

  if (abs(coeff_pvq[0] - ref_coeff_pvq[0]) < pvq_dc_quant * 141/256) { /* 0.55 */
    dqcoeff_pvq[0] = 0;
  }
  else {
    dqcoeff_pvq[0] = OD_DIV_R0(coeff_pvq[0] - ref_coeff_pvq[0], pvq_dc_quant);
  }

  {// for debugging, to match with decoder side ref vector.
  int i;
  for (i=0; i<tx_blk_size*tx_blk_size; i++)
    pvq_info->ref_coeff[i] = ref_coeff[i];
  }

  tell = od_ec_enc_tell(&daala_enc.ec);

  skip = pvq_encode_helper(&daala_enc,    // daala encoder
                           ref_coeff_pvq, // reference vector
                           coeff_pvq,     // target original vector
                           dqcoeff_pvq,   // de-quantized vector
                           quant,         // quantizer
                           plane,         // image plane
                           tx_size,       // transform size in log_2 - 2, ex: 0 is for 4x4
                           0,            // key frame? 0 for now
                           pvq_info);

  // Encode residue of DC coeff, if required.
  if (!has_dc_skip || dqcoeff_pvq[0]) {
    generic_encode(&daala_enc.ec, &daala_enc.state.adapt.model_dc[plane],
     abs(dqcoeff_pvq[0]) - has_dc_skip, -1,
     &daala_enc.state.adapt.ex_dc[plane][tx_size][0], 2);
  }
  if (dqcoeff_pvq[0]) {
    od_ec_enc_bits(&daala_enc.ec, dqcoeff_pvq[0] < 0, 1);
    skip = 0;
  }
  // need to save quantized residue of DC coeff
  // so that final pvq bitstream writing can know whether DC is coded.
  pvq_info->dq_dc_residue = dqcoeff_pvq[0];

  dqcoeff_pvq[0] = dqcoeff_pvq[0] * pvq_dc_quant;
  dqcoeff_pvq[0] += ref_coeff_pvq[0];

  *rate = od_ec_enc_tell(&daala_enc.ec) - tell;

  // Safely initialize dqcoeff since some coeffs (band size > 128 coeffs)
  // are skipped by PVQ.
  od_init_skipped_coeffs(dqcoeff, ref_coeff, 0, 0, tx_blk_size, tx_blk_size);

  // Back to original coefficient order
  od_coding_order_to_raster(dqcoeff, tx_blk_size, dqcoeff_pvq, tx_blk_size);

  // Mark last nonzero coeff position.
  for (j = 0; j < tx_blk_size*tx_blk_size; j++)
    if (dqcoeff[j]) *eob = j + 1;

  pvq_info->eob = *eob;

  if (skip)
    assert(pvq_info->ac_dc_coded == 0);

  if (!skip)
    assert(pvq_info->ac_dc_coded > 0);

  return skip;
}

void store_pvq_enc_info(PVQ_INFO *pvq_info,
                        int *qg,
                        int *theta,
                        int *max_theta,
                        int *k,
                        od_coeff *y,
                        int *exg,
                        int *ext,
                        int nb_bands,
                        int *off,
                        int *size,
                        int skip_rest,
                        int skip_dir,
                        int bs) {       // log of the block size minus two
  int i;

  for (i=0; i < PVQ_MAX_PARTITIONS; i++) {
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
