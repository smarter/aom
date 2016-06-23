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

#ifndef AV1_ENCODER_ENCODEMB_H_
#define AV1_ENCODER_ENCODEMB_H_

#include "./aom_config.h"
#include "av1/encoder/block.h"

#ifdef __cplusplus
extern "C" {
#endif

struct encode_b_args {
  MACROBLOCK *x;
  struct optimize_ctx *ctx;
  int8_t *skip;
};
void av1_encode_sb(MACROBLOCK *x, BLOCK_SIZE bsize);
void av1_encode_sby_pass1(MACROBLOCK *x, BLOCK_SIZE bsize);
void av1_xform_quant_fp(MACROBLOCK *x, int plane, int block, int blk_row,
                        int blk_col, BLOCK_SIZE plane_bsize, TX_SIZE tx_size);
void av1_xform_quant_dc(MACROBLOCK *x, int plane, int block, int blk_row,
                        int blk_col, BLOCK_SIZE plane_bsize, TX_SIZE tx_size);
void av1_xform_quant(MACROBLOCK *x, int plane, int block, int blk_row,
                     int blk_col, BLOCK_SIZE plane_bsize, TX_SIZE tx_size);

void av1_subtract_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane);

void av1_encode_block_intra(int plane, int block, int blk_row, int blk_col,
                            BLOCK_SIZE plane_bsize, TX_SIZE tx_size, void *arg);

void av1_encode_intra_block_plane(MACROBLOCK *x, BLOCK_SIZE bsize, int plane);

void av1_fwd_txfm_4x4(const int16_t *src_diff, tran_low_t *coeff,
                      int diff_stride, TX_TYPE tx_type, int lossless);

void fwd_txfm_8x8(const int16_t *src_diff, tran_low_t *coeff,
                  int diff_stride, TX_TYPE tx_type);

void fwd_txfm_16x16(const int16_t *src_diff, tran_low_t *coeff,
                    int diff_stride, TX_TYPE tx_type);

void fwd_txfm_32x32(int rd_transform, const int16_t *src_diff,
                    tran_low_t *coeff, int diff_stride, TX_TYPE tx_type);

#if CONFIG_PVQ
int pvq_encode_helper(daala_enc_ctx *daala_enc,
    tran_low_t *const coeff, tran_low_t *ref_coeff,
    tran_low_t *const dqcoeff,
    uint16_t *eob, const int16_t *quant,
    int plane, int tx_size, int *rate, PVQ_INFO *pvq_info);

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
                        int bs);
#endif

#if CONFIG_AOM_HIGHBITDEPTH
void av1_highbd_fwd_txfm_4x4(const int16_t *src_diff, tran_low_t *coeff,
                             int diff_stride, TX_TYPE tx_type, int lossless);
#endif  // CONFIG_AOM_HIGHBITDEPTH

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AV1_ENCODER_ENCODEMB_H_
