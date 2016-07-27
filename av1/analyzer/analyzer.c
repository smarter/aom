#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "av1/decoder/decoder.h"
#include "aom/aom_decoder.h"
#include "aom/aomdx.h"
#include "../decoder/decoder.h"

#include "analyzer.h"
#include "../common/blockd.h"
#include "../common/onyxc_int.h"

aom_codec_err_t analyzer_record_predicted_block(struct AV1Decoder *pbi,
                                                uint32_t plane,
                                                uint32_t subsampling_x,
                                                uint32_t subsampling_y,
                                                uint32_t mi_col,
                                                uint32_t mi_row,
                                                uint32_t col,
                                                uint32_t row,
                                                unsigned char* src,
                                                uint32_t src_stride,
                                                uint32_t transform_size) {
  if (pbi->analyzer_data == NULL) {
    return AOM_CODEC_OK;
  }
  AnalyzerImage *image = &pbi->analyzer_data->predicted_image;
  AnalyzerImagePlane *imagePlane = &image->planes[plane];
  if (imagePlane->buffer == NULL) {
    return AOM_CODEC_OK;
  }
  int size = 4 << transform_size;
  int mi_size_y = MI_SIZE >> subsampling_y;
  int mi_size_x = MI_SIZE >> subsampling_x;
  size_t offset = mi_row * mi_size_y * imagePlane->stride + mi_col * mi_size_x +
               row * (mi_size_y >> 1) * imagePlane->stride + col * (mi_size_x >> 1);
  if (offset > imagePlane->size) {
    return AOM_CODEC_ERROR;
  }

  unsigned char* dst = imagePlane->buffer + offset;
  for (int i = 0; i < size; i++) {
    memcpy(dst, src, size);
    dst += imagePlane->stride;
    src += src_stride;
  }
  return AOM_CODEC_OK;
}

// Saves the decoder state.
aom_codec_err_t analyzer_record_frame(struct AV1Decoder *pbi) {
  AV1_COMMON *const cm = &pbi->common;
  const uint32_t mi_rows = cm->mi_rows;
  const uint32_t mi_cols = cm->mi_cols;
  uint32_t r, c;
  if (pbi->analyzer_data == NULL) {
    return AOM_CODEC_OK;
  }

  pbi->analyzer_data->mi_rows = mi_rows;
  pbi->analyzer_data->mi_cols = mi_cols;

  pbi->analyzer_data->tile_rows_log2 = cm->log2_tile_rows;
  pbi->analyzer_data->tile_cols_log2 = cm->log2_tile_cols;

  // Save mode info.
  AnalyzerMIBuffer mi_grid = pbi->analyzer_data->mi_grid;
  if (mi_grid.length > 0) {
    if (mi_rows * mi_cols > mi_grid.length) {
      return AOM_CODEC_ERROR;
    }
    for (r = 0; r < mi_rows; ++r) {
      for (c = 0; c < mi_cols; ++c) {
        const MB_MODE_INFO *mbmi =
          &cm->mi_grid_visible[r * cm->mi_stride + c]->mbmi;
        AnalyzerMI *mi = &mi_grid.buffer[r * mi_cols + c];
        // MVs
        mi->mv[0].col = mbmi->mv[0].as_mv.col;
        mi->mv[0].row = mbmi->mv[0].as_mv.row;
        mi->mv[1].col = mbmi->mv[1].as_mv.col;
        mi->mv[1].row = mbmi->mv[1].as_mv.row;
        // Reference Frames
        mi->reference_frame[0] = mbmi->ref_frame[0];
        mi->reference_frame[1] = mbmi->ref_frame[1];
        // Prediction Mode
        mi->mode = mbmi->mode;
        // Deringing Gain
        mi->dering_gain = mbmi->dering_gain;
        // Block size
        mi->block_size = mbmi->sb_type;
        // Skip flag
        mi->skip = mbmi->skip;
        // Filters
        mi->filter = mbmi->interp_filter;
        // Transform
        mi->transform_type = mbmi->tx_type;
        mi->transform_size = mbmi->tx_size;
      }
    }
  }
  return AOM_CODEC_OK;
}

aom_codec_err_t analyzer_record_mi_bits(struct AV1Decoder *pbi,
                                        uint32_t plane,
                                        uint32_t mi_col,
                                        uint32_t mi_row,
                                        uint32_t bits) {
  AV1_COMMON *const cm = &pbi->common;
  if (pbi->analyzer_data == NULL) {
    return AOM_CODEC_OK;
  }
  if (plane != 0) {
    return AOM_CODEC_OK;
  }
  pbi->analyzer_data->mi_grid.buffer[mi_row * cm->mi_cols + mi_col].bits += bits;
  return AOM_CODEC_OK;
}
