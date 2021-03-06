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

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <assert.h>

#include "./av1_rtcd.h"
#include "av1/common/warped_motion.h"

/* clang-format off */
static const int error_measure_lut[512] = {
  // pow 0.7
  16384, 16339, 16294, 16249, 16204, 16158, 16113, 16068,
  16022, 15977, 15932, 15886, 15840, 15795, 15749, 15703,
  15657, 15612, 15566, 15520, 15474, 15427, 15381, 15335,
  15289, 15242, 15196, 15149, 15103, 15056, 15010, 14963,
  14916, 14869, 14822, 14775, 14728, 14681, 14634, 14587,
  14539, 14492, 14445, 14397, 14350, 14302, 14254, 14206,
  14159, 14111, 14063, 14015, 13967, 13918, 13870, 13822,
  13773, 13725, 13676, 13628, 13579, 13530, 13481, 13432,
  13383, 13334, 13285, 13236, 13187, 13137, 13088, 13038,
  12988, 12939, 12889, 12839, 12789, 12739, 12689, 12639,
  12588, 12538, 12487, 12437, 12386, 12335, 12285, 12234,
  12183, 12132, 12080, 12029, 11978, 11926, 11875, 11823,
  11771, 11719, 11667, 11615, 11563, 11511, 11458, 11406,
  11353, 11301, 11248, 11195, 11142, 11089, 11036, 10982,
  10929, 10875, 10822, 10768, 10714, 10660, 10606, 10552,
  10497, 10443, 10388, 10333, 10279, 10224, 10168, 10113,
  10058, 10002,  9947,  9891,  9835,  9779,  9723,  9666,
  9610, 9553, 9497, 9440, 9383, 9326, 9268, 9211,
  9153, 9095, 9037, 8979, 8921, 8862, 8804, 8745,
  8686, 8627, 8568, 8508, 8449, 8389, 8329, 8269,
  8208, 8148, 8087, 8026, 7965, 7903, 7842, 7780,
  7718, 7656, 7593, 7531, 7468, 7405, 7341, 7278,
  7214, 7150, 7086, 7021, 6956, 6891, 6826, 6760,
  6695, 6628, 6562, 6495, 6428, 6361, 6293, 6225,
  6157, 6089, 6020, 5950, 5881, 5811, 5741, 5670,
  5599, 5527, 5456, 5383, 5311, 5237, 5164, 5090,
  5015, 4941, 4865, 4789, 4713, 4636, 4558, 4480,
  4401, 4322, 4242, 4162, 4080, 3998, 3916, 3832,
  3748, 3663, 3577, 3490, 3402, 3314, 3224, 3133,
  3041, 2948, 2854, 2758, 2661, 2562, 2461, 2359,
  2255, 2148, 2040, 1929, 1815, 1698, 1577, 1452,
  1323, 1187, 1045,  894,  731,  550,  339,    0,
  339,  550,  731,  894, 1045, 1187, 1323, 1452,
  1577, 1698, 1815, 1929, 2040, 2148, 2255, 2359,
  2461, 2562, 2661, 2758, 2854, 2948, 3041, 3133,
  3224, 3314, 3402, 3490, 3577, 3663, 3748, 3832,
  3916, 3998, 4080, 4162, 4242, 4322, 4401, 4480,
  4558, 4636, 4713, 4789, 4865, 4941, 5015, 5090,
  5164, 5237, 5311, 5383, 5456, 5527, 5599, 5670,
  5741, 5811, 5881, 5950, 6020, 6089, 6157, 6225,
  6293, 6361, 6428, 6495, 6562, 6628, 6695, 6760,
  6826, 6891, 6956, 7021, 7086, 7150, 7214, 7278,
  7341, 7405, 7468, 7531, 7593, 7656, 7718, 7780,
  7842, 7903, 7965, 8026, 8087, 8148, 8208, 8269,
  8329, 8389, 8449, 8508, 8568, 8627, 8686, 8745,
  8804, 8862, 8921, 8979, 9037, 9095, 9153, 9211,
  9268, 9326, 9383, 9440, 9497, 9553, 9610, 9666,
  9723,  9779,  9835,  9891,  9947, 10002, 10058, 10113,
  10168, 10224, 10279, 10333, 10388, 10443, 10497, 10552,
  10606, 10660, 10714, 10768, 10822, 10875, 10929, 10982,
  11036, 11089, 11142, 11195, 11248, 11301, 11353, 11406,
  11458, 11511, 11563, 11615, 11667, 11719, 11771, 11823,
  11875, 11926, 11978, 12029, 12080, 12132, 12183, 12234,
  12285, 12335, 12386, 12437, 12487, 12538, 12588, 12639,
  12689, 12739, 12789, 12839, 12889, 12939, 12988, 13038,
  13088, 13137, 13187, 13236, 13285, 13334, 13383, 13432,
  13481, 13530, 13579, 13628, 13676, 13725, 13773, 13822,
  13870, 13918, 13967, 14015, 14063, 14111, 14159, 14206,
  14254, 14302, 14350, 14397, 14445, 14492, 14539, 14587,
  14634, 14681, 14728, 14775, 14822, 14869, 14916, 14963,
  15010, 15056, 15103, 15149, 15196, 15242, 15289, 15335,
  15381, 15427, 15474, 15520, 15566, 15612, 15657, 15703,
  15749, 15795, 15840, 15886, 15932, 15977, 16022, 16068,
  16113, 16158, 16204, 16249, 16294, 16339, 16384, 16384,
};
/* clang-format on */

static ProjectPointsFunc get_project_points_type(TransformationType type) {
  switch (type) {
    case HOMOGRAPHY: return project_points_homography;
    case AFFINE: return project_points_affine;
    case ROTZOOM: return project_points_rotzoom;
    case TRANSLATION: return project_points_translation;
    default: assert(0); return NULL;
  }
}

void project_points_translation(int32_t *mat, int *points, int *proj,
                                const int n, const int stride_points,
                                const int stride_proj, const int subsampling_x,
                                const int subsampling_y) {
  int i;
  for (i = 0; i < n; ++i) {
    const int x = *(points++), y = *(points++);
    if (subsampling_x)
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          ((x * (1 << (WARPEDMODEL_PREC_BITS + 1))) + mat[0]),
          WARPEDDIFF_PREC_BITS + 1);
    else
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          ((x * (1 << WARPEDMODEL_PREC_BITS)) + mat[0]), WARPEDDIFF_PREC_BITS);
    if (subsampling_y)
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          ((y * (1 << (WARPEDMODEL_PREC_BITS + 1))) + mat[1]),
          WARPEDDIFF_PREC_BITS + 1);
    else
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          ((y * (1 << WARPEDMODEL_PREC_BITS))) + mat[1], WARPEDDIFF_PREC_BITS);
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

void project_points_rotzoom(int32_t *mat, int *points, int *proj, const int n,
                            const int stride_points, const int stride_proj,
                            const int subsampling_x, const int subsampling_y) {
  int i;
  for (i = 0; i < n; ++i) {
    const int x = *(points++), y = *(points++);
    if (subsampling_x)
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          mat[2] * 2 * x + mat[3] * 2 * y + mat[0] +
              (mat[2] + mat[3] - (1 << WARPEDMODEL_PREC_BITS)) / 2,
          WARPEDDIFF_PREC_BITS + 1);
    else
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(mat[2] * x + mat[3] * y + mat[0],
                                            WARPEDDIFF_PREC_BITS);
    if (subsampling_y)
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          -mat[3] * 2 * x + mat[2] * 2 * y + mat[1] +
              (-mat[3] + mat[2] - (1 << WARPEDMODEL_PREC_BITS)) / 2,
          WARPEDDIFF_PREC_BITS + 1);
    else
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(-mat[3] * x + mat[2] * y + mat[1],
                                            WARPEDDIFF_PREC_BITS);
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

void project_points_affine(int32_t *mat, int *points, int *proj, const int n,
                           const int stride_points, const int stride_proj,
                           const int subsampling_x, const int subsampling_y) {
  int i;
  for (i = 0; i < n; ++i) {
    const int x = *(points++), y = *(points++);
    if (subsampling_x)
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          mat[2] * 2 * x + mat[3] * 2 * y + mat[0] +
              (mat[2] + mat[3] - (1 << WARPEDMODEL_PREC_BITS)) / 2,
          WARPEDDIFF_PREC_BITS + 1);
    else
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(mat[2] * x + mat[3] * y + mat[0],
                                            WARPEDDIFF_PREC_BITS);
    if (subsampling_y)
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(
          mat[4] * 2 * x + mat[5] * 2 * y + mat[1] +
              (mat[4] + mat[5] - (1 << WARPEDMODEL_PREC_BITS)) / 2,
          WARPEDDIFF_PREC_BITS + 1);
    else
      *(proj++) = ROUND_POWER_OF_TWO_SIGNED(mat[4] * x + mat[5] * y + mat[1],
                                            WARPEDDIFF_PREC_BITS);
    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

void project_points_hortrapezoid(int32_t *mat, int *points, int *proj,
                                 const int n, const int stride_points,
                                 const int stride_proj, const int subsampling_x,
                                 const int subsampling_y) {
  int i;
  int64_t x, y, Z;
  int64_t xp, yp;
  for (i = 0; i < n; ++i) {
    x = *(points++), y = *(points++);
    x = (subsampling_x ? 4 * x + 1 : 2 * x);
    y = (subsampling_y ? 4 * y + 1 : 2 * y);

    Z = (mat[7] * y + (1 << (WARPEDMODEL_ROW3HOMO_PREC_BITS + 1)));
    xp = (mat[2] * x + mat[3] * y + 2 * mat[0]) *
         (1 << (WARPEDPIXEL_PREC_BITS + WARPEDMODEL_ROW3HOMO_PREC_BITS -
                WARPEDMODEL_PREC_BITS));
    yp = (mat[5] * y + 2 * mat[1]) *
         (1 << (WARPEDPIXEL_PREC_BITS + WARPEDMODEL_ROW3HOMO_PREC_BITS -
                WARPEDMODEL_PREC_BITS));

    xp = xp > 0 ? (xp + Z / 2) / Z : (xp - Z / 2) / Z;
    yp = yp > 0 ? (yp + Z / 2) / Z : (yp - Z / 2) / Z;

    if (subsampling_x) xp = (xp - (1 << (WARPEDPIXEL_PREC_BITS - 1))) / 2;
    if (subsampling_y) yp = (yp - (1 << (WARPEDPIXEL_PREC_BITS - 1))) / 2;
    *(proj++) = xp;
    *(proj++) = yp;

    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

void project_points_vertrapezoid(int32_t *mat, int *points, int *proj,
                                 const int n, const int stride_points,
                                 const int stride_proj, const int subsampling_x,
                                 const int subsampling_y) {
  int i;
  int64_t x, y, Z;
  int64_t xp, yp;
  for (i = 0; i < n; ++i) {
    x = *(points++), y = *(points++);
    x = (subsampling_x ? 4 * x + 1 : 2 * x);
    y = (subsampling_y ? 4 * y + 1 : 2 * y);

    Z = (mat[6] * x + (1 << (WARPEDMODEL_ROW3HOMO_PREC_BITS + 1)));
    xp = (mat[2] * x + 2 * mat[0]) *
         (1 << (WARPEDPIXEL_PREC_BITS + WARPEDMODEL_ROW3HOMO_PREC_BITS -
                WARPEDMODEL_PREC_BITS));
    yp = (mat[4] * x + mat[5] * y + 2 * mat[1]) *
         (1 << (WARPEDPIXEL_PREC_BITS + WARPEDMODEL_ROW3HOMO_PREC_BITS -
                WARPEDMODEL_PREC_BITS));

    xp = xp > 0 ? (xp + Z / 2) / Z : (xp - Z / 2) / Z;
    yp = yp > 0 ? (yp + Z / 2) / Z : (yp - Z / 2) / Z;

    if (subsampling_x) xp = (xp - (1 << (WARPEDPIXEL_PREC_BITS - 1))) / 2;
    if (subsampling_y) yp = (yp - (1 << (WARPEDPIXEL_PREC_BITS - 1))) / 2;
    *(proj++) = xp;
    *(proj++) = yp;

    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

void project_points_homography(int32_t *mat, int *points, int *proj,
                               const int n, const int stride_points,
                               const int stride_proj, const int subsampling_x,
                               const int subsampling_y) {
  int i;
  int64_t x, y, Z;
  int64_t xp, yp;
  for (i = 0; i < n; ++i) {
    x = *(points++), y = *(points++);
    x = (subsampling_x ? 4 * x + 1 : 2 * x);
    y = (subsampling_y ? 4 * y + 1 : 2 * y);

    Z = (mat[6] * x + mat[7] * y + (1 << (WARPEDMODEL_ROW3HOMO_PREC_BITS + 1)));
    xp = (mat[2] * x + mat[3] * y + 2 * mat[0]) *
         (1 << (WARPEDPIXEL_PREC_BITS + WARPEDMODEL_ROW3HOMO_PREC_BITS -
                WARPEDMODEL_PREC_BITS));
    yp = (mat[4] * x + mat[5] * y + 2 * mat[1]) *
         (1 << (WARPEDPIXEL_PREC_BITS + WARPEDMODEL_ROW3HOMO_PREC_BITS -
                WARPEDMODEL_PREC_BITS));

    xp = xp > 0 ? (xp + Z / 2) / Z : (xp - Z / 2) / Z;
    yp = yp > 0 ? (yp + Z / 2) / Z : (yp - Z / 2) / Z;

    if (subsampling_x) xp = (xp - (1 << (WARPEDPIXEL_PREC_BITS - 1))) / 2;
    if (subsampling_y) yp = (yp - (1 << (WARPEDPIXEL_PREC_BITS - 1))) / 2;
    *(proj++) = xp;
    *(proj++) = yp;

    points += stride_points - 2;
    proj += stride_proj - 2;
  }
}

// 'points' are at original scale, output 'proj's are scaled up by
// 1 << WARPEDPIXEL_PREC_BITS
void project_points(WarpedMotionParams *wm_params, int *points, int *proj,
                    const int n, const int stride_points, const int stride_proj,
                    const int subsampling_x, const int subsampling_y) {
  switch (wm_params->wmtype) {
    case AFFINE:
      project_points_affine(wm_params->wmmat, points, proj, n, stride_points,
                            stride_proj, subsampling_x, subsampling_y);
      break;
    case ROTZOOM:
      project_points_rotzoom(wm_params->wmmat, points, proj, n, stride_points,
                             stride_proj, subsampling_x, subsampling_y);
      break;
    case HOMOGRAPHY:
      project_points_homography(wm_params->wmmat, points, proj, n,
                                stride_points, stride_proj, subsampling_x,
                                subsampling_y);
      break;
    default: assert(0 && "Invalid warped motion type!"); return;
  }
}

static const int16_t
    filter_ntap[WARPEDPIXEL_PREC_SHIFTS][WARPEDPIXEL_FILTER_TAPS] = {
      { 0, 0, 128, 0, 0, 0 },      { 0, -1, 128, 2, -1, 0 },
      { 1, -3, 127, 4, -1, 0 },    { 1, -4, 126, 6, -2, 1 },
      { 1, -5, 126, 8, -3, 1 },    { 1, -6, 125, 11, -4, 1 },
      { 1, -7, 124, 13, -4, 1 },   { 2, -8, 123, 15, -5, 1 },
      { 2, -9, 122, 18, -6, 1 },   { 2, -10, 121, 20, -6, 1 },
      { 2, -11, 120, 22, -7, 2 },  { 2, -12, 119, 25, -8, 2 },
      { 3, -13, 117, 27, -8, 2 },  { 3, -13, 116, 29, -9, 2 },
      { 3, -14, 114, 32, -10, 3 }, { 3, -15, 113, 35, -10, 2 },
      { 3, -15, 111, 37, -11, 3 }, { 3, -16, 109, 40, -11, 3 },
      { 3, -16, 108, 42, -12, 3 }, { 4, -17, 106, 45, -13, 3 },
      { 4, -17, 104, 47, -13, 3 }, { 4, -17, 102, 50, -14, 3 },
      { 4, -17, 100, 52, -14, 3 }, { 4, -18, 98, 55, -15, 4 },
      { 4, -18, 96, 58, -15, 3 },  { 4, -18, 94, 60, -16, 4 },
      { 4, -18, 91, 63, -16, 4 },  { 4, -18, 89, 65, -16, 4 },
      { 4, -18, 87, 68, -17, 4 },  { 4, -18, 85, 70, -17, 4 },
      { 4, -18, 82, 73, -17, 4 },  { 4, -18, 80, 75, -17, 4 },
      { 4, -18, 78, 78, -18, 4 },  { 4, -17, 75, 80, -18, 4 },
      { 4, -17, 73, 82, -18, 4 },  { 4, -17, 70, 85, -18, 4 },
      { 4, -17, 68, 87, -18, 4 },  { 4, -16, 65, 89, -18, 4 },
      { 4, -16, 63, 91, -18, 4 },  { 4, -16, 60, 94, -18, 4 },
      { 3, -15, 58, 96, -18, 4 },  { 4, -15, 55, 98, -18, 4 },
      { 3, -14, 52, 100, -17, 4 }, { 3, -14, 50, 102, -17, 4 },
      { 3, -13, 47, 104, -17, 4 }, { 3, -13, 45, 106, -17, 4 },
      { 3, -12, 42, 108, -16, 3 }, { 3, -11, 40, 109, -16, 3 },
      { 3, -11, 37, 111, -15, 3 }, { 2, -10, 35, 113, -15, 3 },
      { 3, -10, 32, 114, -14, 3 }, { 2, -9, 29, 116, -13, 3 },
      { 2, -8, 27, 117, -13, 3 },  { 2, -8, 25, 119, -12, 2 },
      { 2, -7, 22, 120, -11, 2 },  { 1, -6, 20, 121, -10, 2 },
      { 1, -6, 18, 122, -9, 2 },   { 1, -5, 15, 123, -8, 2 },
      { 1, -4, 13, 124, -7, 1 },   { 1, -4, 11, 125, -6, 1 },
      { 1, -3, 8, 126, -5, 1 },    { 1, -2, 6, 126, -4, 1 },
      { 0, -1, 4, 127, -3, 1 },    { 0, -1, 2, 128, -1, 0 },
    };

static int32_t do_ntap_filter(int32_t *p, int x) {
  int i;
  int32_t sum = 0;
  for (i = 0; i < WARPEDPIXEL_FILTER_TAPS; ++i) {
    sum += p[i - WARPEDPIXEL_FILTER_TAPS / 2 + 1] * filter_ntap[x][i];
  }
  return sum;
}

static int32_t do_cubic_filter(int32_t *p, int x) {
  if (x == 0) {
    return p[0] * (1 << WARPEDPIXEL_FILTER_BITS);
  } else if (x == (1 << WARPEDPIXEL_PREC_BITS)) {
    return p[1] * (1 << WARPEDPIXEL_FILTER_BITS);
  } else {
    const int64_t v1 = (int64_t)x * x * x * (3 * (p[0] - p[1]) + p[2] - p[-1]);
    const int64_t v2 = x * x * (2 * p[-1] - 5 * p[0] + 4 * p[1] - p[2]);
    const int64_t v3 = x * (p[1] - p[-1]);
    const int64_t v4 = 2 * p[0];
    return (int32_t)ROUND_POWER_OF_TWO_SIGNED(
        (v4 * (1 << (3 * WARPEDPIXEL_PREC_BITS))) +
            (v3 * (1 << (2 * WARPEDPIXEL_PREC_BITS))) +
            (v2 * (1 << WARPEDPIXEL_PREC_BITS)) + v1,
        3 * WARPEDPIXEL_PREC_BITS + 1 - WARPEDPIXEL_FILTER_BITS);
  }
}

static INLINE void get_subcolumn(int taps, uint8_t *ref, int32_t *col,
                                 int stride, int x, int y_start) {
  int i;
  for (i = 0; i < taps; ++i) {
    col[i] = ref[(i + y_start) * stride + x];
  }
}

static uint8_t bi_ntap_filter(uint8_t *ref, int x, int y, int stride) {
  int32_t val, arr[WARPEDPIXEL_FILTER_TAPS];
  int k;
  int i = (int)x >> WARPEDPIXEL_PREC_BITS;
  int j = (int)y >> WARPEDPIXEL_PREC_BITS;
  for (k = 0; k < WARPEDPIXEL_FILTER_TAPS; ++k) {
    int32_t arr_temp[WARPEDPIXEL_FILTER_TAPS];
    get_subcolumn(WARPEDPIXEL_FILTER_TAPS, ref, arr_temp, stride,
                  i + k + 1 - WARPEDPIXEL_FILTER_TAPS / 2,
                  j + 1 - WARPEDPIXEL_FILTER_TAPS / 2);
    arr[k] = do_ntap_filter(arr_temp + WARPEDPIXEL_FILTER_TAPS / 2 - 1,
                            y - (j * (1 << WARPEDPIXEL_PREC_BITS)));
  }
  val = do_ntap_filter(arr + WARPEDPIXEL_FILTER_TAPS / 2 - 1,
                       x - (i * (1 << WARPEDPIXEL_PREC_BITS)));
  val = ROUND_POWER_OF_TWO_SIGNED(val, WARPEDPIXEL_FILTER_BITS * 2);
  return (uint8_t)clip_pixel(val);
}

static uint8_t bi_cubic_filter(uint8_t *ref, int x, int y, int stride) {
  int32_t val, arr[4];
  int k;
  int i = (int)x >> WARPEDPIXEL_PREC_BITS;
  int j = (int)y >> WARPEDPIXEL_PREC_BITS;
  for (k = 0; k < 4; ++k) {
    int32_t arr_temp[4];
    get_subcolumn(4, ref, arr_temp, stride, i + k - 1, j - 1);
    arr[k] =
        do_cubic_filter(arr_temp + 1, y - (j * (1 << WARPEDPIXEL_PREC_BITS)));
  }
  val = do_cubic_filter(arr + 1, x - (i * (1 << WARPEDPIXEL_PREC_BITS)));
  val = ROUND_POWER_OF_TWO_SIGNED(val, WARPEDPIXEL_FILTER_BITS * 2);
  return (uint8_t)clip_pixel(val);
}

static uint8_t bi_linear_filter(uint8_t *ref, int x, int y, int stride) {
  const int ix = x >> WARPEDPIXEL_PREC_BITS;
  const int iy = y >> WARPEDPIXEL_PREC_BITS;
  const int sx = x - (ix * (1 << WARPEDPIXEL_PREC_BITS));
  const int sy = y - (iy * (1 << WARPEDPIXEL_PREC_BITS));
  int32_t val;
  val = ROUND_POWER_OF_TWO_SIGNED(
      ref[iy * stride + ix] * (WARPEDPIXEL_PREC_SHIFTS - sy) *
              (WARPEDPIXEL_PREC_SHIFTS - sx) +
          ref[iy * stride + ix + 1] * (WARPEDPIXEL_PREC_SHIFTS - sy) * sx +
          ref[(iy + 1) * stride + ix] * sy * (WARPEDPIXEL_PREC_SHIFTS - sx) +
          ref[(iy + 1) * stride + ix + 1] * sy * sx,
      WARPEDPIXEL_PREC_BITS * 2);
  return (uint8_t)clip_pixel(val);
}

static uint8_t warp_interpolate(uint8_t *ref, int x, int y, int width,
                                int height, int stride) {
  int ix = x >> WARPEDPIXEL_PREC_BITS;
  int iy = y >> WARPEDPIXEL_PREC_BITS;
  int sx = x - (ix * (1 << WARPEDPIXEL_PREC_BITS));
  int sy = y - (iy * (1 << WARPEDPIXEL_PREC_BITS));
  int32_t v;

  if (ix < 0 && iy < 0)
    return ref[0];
  else if (ix < 0 && iy >= height - 1)
    return ref[(height - 1) * stride];
  else if (ix >= width - 1 && iy < 0)
    return ref[width - 1];
  else if (ix >= width - 1 && iy >= height - 1)
    return ref[(height - 1) * stride + (width - 1)];
  else if (ix < 0) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[iy * stride] * (WARPEDPIXEL_PREC_SHIFTS - sy) +
            ref[(iy + 1) * stride] * sy,
        WARPEDPIXEL_PREC_BITS);
    return clip_pixel(v);
  } else if (iy < 0) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[ix] * (WARPEDPIXEL_PREC_SHIFTS - sx) + ref[ix + 1] * sx,
        WARPEDPIXEL_PREC_BITS);
    return clip_pixel(v);
  } else if (ix >= width - 1) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[iy * stride + width - 1] * (WARPEDPIXEL_PREC_SHIFTS - sy) +
            ref[(iy + 1) * stride + width - 1] * sy,
        WARPEDPIXEL_PREC_BITS);
    return clip_pixel(v);
  } else if (iy >= height - 1) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[(height - 1) * stride + ix] * (WARPEDPIXEL_PREC_SHIFTS - sx) +
            ref[(height - 1) * stride + ix + 1] * sx,
        WARPEDPIXEL_PREC_BITS);
    return clip_pixel(v);
  } else if (ix >= WARPEDPIXEL_FILTER_TAPS / 2 - 1 &&
             iy >= WARPEDPIXEL_FILTER_TAPS / 2 - 1 &&
             ix < width - WARPEDPIXEL_FILTER_TAPS / 2 &&
             iy < height - WARPEDPIXEL_FILTER_TAPS / 2) {
    return bi_ntap_filter(ref, x, y, stride);
  } else if (ix >= 1 && iy >= 1 && ix < width - 2 && iy < height - 2) {
    return bi_cubic_filter(ref, x, y, stride);
  } else {
    return bi_linear_filter(ref, x, y, stride);
  }
}

// For warping, we really use a 6-tap filter, but we do blocks of 8 pixels
// at a time. The zoom/rotation/shear in the model are applied to the
// "fractional" position of each pixel, which therefore varies within
// [-1, 2) * WARPEDPIXEL_PREC_SHIFTS.
// We need an extra 2 taps to fit this in, for a total of 8 taps.
/* clang-format off */
const int16_t warped_filter[WARPEDPIXEL_PREC_SHIFTS * 3 + 1][8] = {
  // [-1, 0)
  { 0,   0, 128,   0,   0, 0, 0, 0 }, { 0, - 1, 128,   2, - 1, 0, 0, 0 },
  { 1, - 3, 127,   4, - 1, 0, 0, 0 }, { 1, - 4, 126,   6, - 2, 1, 0, 0 },
  { 1, - 5, 126,   8, - 3, 1, 0, 0 }, { 1, - 6, 125,  11, - 4, 1, 0, 0 },
  { 1, - 7, 124,  13, - 4, 1, 0, 0 }, { 2, - 8, 123,  15, - 5, 1, 0, 0 },
  { 2, - 9, 122,  18, - 6, 1, 0, 0 }, { 2, -10, 121,  20, - 6, 1, 0, 0 },
  { 2, -11, 120,  22, - 7, 2, 0, 0 }, { 2, -12, 119,  25, - 8, 2, 0, 0 },
  { 3, -13, 117,  27, - 8, 2, 0, 0 }, { 3, -13, 116,  29, - 9, 2, 0, 0 },
  { 3, -14, 114,  32, -10, 3, 0, 0 }, { 3, -15, 113,  35, -10, 2, 0, 0 },
  { 3, -15, 111,  37, -11, 3, 0, 0 }, { 3, -16, 109,  40, -11, 3, 0, 0 },
  { 3, -16, 108,  42, -12, 3, 0, 0 }, { 4, -17, 106,  45, -13, 3, 0, 0 },
  { 4, -17, 104,  47, -13, 3, 0, 0 }, { 4, -17, 102,  50, -14, 3, 0, 0 },
  { 4, -17, 100,  52, -14, 3, 0, 0 }, { 4, -18,  98,  55, -15, 4, 0, 0 },
  { 4, -18,  96,  58, -15, 3, 0, 0 }, { 4, -18,  94,  60, -16, 4, 0, 0 },
  { 4, -18,  91,  63, -16, 4, 0, 0 }, { 4, -18,  89,  65, -16, 4, 0, 0 },
  { 4, -18,  87,  68, -17, 4, 0, 0 }, { 4, -18,  85,  70, -17, 4, 0, 0 },
  { 4, -18,  82,  73, -17, 4, 0, 0 }, { 4, -18,  80,  75, -17, 4, 0, 0 },
  { 4, -18,  78,  78, -18, 4, 0, 0 }, { 4, -17,  75,  80, -18, 4, 0, 0 },
  { 4, -17,  73,  82, -18, 4, 0, 0 }, { 4, -17,  70,  85, -18, 4, 0, 0 },
  { 4, -17,  68,  87, -18, 4, 0, 0 }, { 4, -16,  65,  89, -18, 4, 0, 0 },
  { 4, -16,  63,  91, -18, 4, 0, 0 }, { 4, -16,  60,  94, -18, 4, 0, 0 },
  { 3, -15,  58,  96, -18, 4, 0, 0 }, { 4, -15,  55,  98, -18, 4, 0, 0 },
  { 3, -14,  52, 100, -17, 4, 0, 0 }, { 3, -14,  50, 102, -17, 4, 0, 0 },
  { 3, -13,  47, 104, -17, 4, 0, 0 }, { 3, -13,  45, 106, -17, 4, 0, 0 },
  { 3, -12,  42, 108, -16, 3, 0, 0 }, { 3, -11,  40, 109, -16, 3, 0, 0 },
  { 3, -11,  37, 111, -15, 3, 0, 0 }, { 2, -10,  35, 113, -15, 3, 0, 0 },
  { 3, -10,  32, 114, -14, 3, 0, 0 }, { 2, - 9,  29, 116, -13, 3, 0, 0 },
  { 2, - 8,  27, 117, -13, 3, 0, 0 }, { 2, - 8,  25, 119, -12, 2, 0, 0 },
  { 2, - 7,  22, 120, -11, 2, 0, 0 }, { 1, - 6,  20, 121, -10, 2, 0, 0 },
  { 1, - 6,  18, 122, - 9, 2, 0, 0 }, { 1, - 5,  15, 123, - 8, 2, 0, 0 },
  { 1, - 4,  13, 124, - 7, 1, 0, 0 }, { 1, - 4,  11, 125, - 6, 1, 0, 0 },
  { 1, - 3,   8, 126, - 5, 1, 0, 0 }, { 1, - 2,   6, 126, - 4, 1, 0, 0 },
  { 0, - 1,   4, 127, - 3, 1, 0, 0 }, { 0, - 1,   2, 128, - 1, 0, 0, 0 },

  // [0, 1)
  { 0,  0,   0, 128,   0,   0,  0,  0}, { 0,  1,  -2, 128,   2,  -1,  0,  0},
  { 0,  1,  -3, 127,   4,  -2,  1,  0}, { 0,  1,  -5, 127,   6,  -2,  1,  0},
  { 0,  2,  -6, 126,   8,  -3,  1,  0}, {-1,  2,  -7, 126,  11,  -4,  2, -1},
  {-1,  3,  -8, 125,  13,  -5,  2, -1}, {-1,  3, -10, 124,  16,  -6,  3, -1},
  {-1,  4, -11, 123,  18,  -7,  3, -1}, {-1,  4, -12, 122,  20,  -7,  3, -1},
  {-1,  4, -13, 121,  23,  -8,  3, -1}, {-2,  5, -14, 120,  25,  -9,  4, -1},
  {-1,  5, -15, 119,  27, -10,  4, -1}, {-1,  5, -16, 118,  30, -11,  4, -1},
  {-2,  6, -17, 116,  33, -12,  5, -1}, {-2,  6, -17, 114,  35, -12,  5, -1},
  {-2,  6, -18, 113,  38, -13,  5, -1}, {-2,  7, -19, 111,  41, -14,  6, -2},
  {-2,  7, -19, 110,  43, -15,  6, -2}, {-2,  7, -20, 108,  46, -15,  6, -2},
  {-2,  7, -20, 106,  49, -16,  6, -2}, {-2,  7, -21, 104,  51, -16,  7, -2},
  {-2,  7, -21, 102,  54, -17,  7, -2}, {-2,  8, -21, 100,  56, -18,  7, -2},
  {-2,  8, -22,  98,  59, -18,  7, -2}, {-2,  8, -22,  96,  62, -19,  7, -2},
  {-2,  8, -22,  94,  64, -19,  7, -2}, {-2,  8, -22,  91,  67, -20,  8, -2},
  {-2,  8, -22,  89,  69, -20,  8, -2}, {-2,  8, -22,  87,  72, -21,  8, -2},
  {-2,  8, -21,  84,  74, -21,  8, -2}, {-2,  8, -22,  82,  77, -21,  8, -2},
  {-2,  8, -21,  79,  79, -21,  8, -2}, {-2,  8, -21,  77,  82, -22,  8, -2},
  {-2,  8, -21,  74,  84, -21,  8, -2}, {-2,  8, -21,  72,  87, -22,  8, -2},
  {-2,  8, -20,  69,  89, -22,  8, -2}, {-2,  8, -20,  67,  91, -22,  8, -2},
  {-2,  7, -19,  64,  94, -22,  8, -2}, {-2,  7, -19,  62,  96, -22,  8, -2},
  {-2,  7, -18,  59,  98, -22,  8, -2}, {-2,  7, -18,  56, 100, -21,  8, -2},
  {-2,  7, -17,  54, 102, -21,  7, -2}, {-2,  7, -16,  51, 104, -21,  7, -2},
  {-2,  6, -16,  49, 106, -20,  7, -2}, {-2,  6, -15,  46, 108, -20,  7, -2},
  {-2,  6, -15,  43, 110, -19,  7, -2}, {-2,  6, -14,  41, 111, -19,  7, -2},
  {-1,  5, -13,  38, 113, -18,  6, -2}, {-1,  5, -12,  35, 114, -17,  6, -2},
  {-1,  5, -12,  33, 116, -17,  6, -2}, {-1,  4, -11,  30, 118, -16,  5, -1},
  {-1,  4, -10,  27, 119, -15,  5, -1}, {-1,  4,  -9,  25, 120, -14,  5, -2},
  {-1,  3,  -8,  23, 121, -13,  4, -1}, {-1,  3,  -7,  20, 122, -12,  4, -1},
  {-1,  3,  -7,  18, 123, -11,  4, -1}, {-1,  3,  -6,  16, 124, -10,  3, -1},
  {-1,  2,  -5,  13, 125,  -8,  3, -1}, {-1,  2,  -4,  11, 126,  -7,  2, -1},
  { 0,  1,  -3,   8, 126,  -6,  2,  0}, { 0,  1,  -2,   6, 127,  -5,  1,  0},
  { 0,  1,  -2,   4, 127,  -3,  1,  0}, { 0,  0,  -1,   2, 128,  -2,  1,  0},

  // [1, 2)
  { 0, 0, 0,   0, 128,   0,   0, 0 }, { 0, 0, 0, - 1, 128,   2, - 1, 0 },
  { 0, 0, 1, - 3, 127,   4, - 1, 0 }, { 0, 0, 1, - 4, 126,   6, - 2, 1 },
  { 0, 0, 1, - 5, 126,   8, - 3, 1 }, { 0, 0, 1, - 6, 125,  11, - 4, 1 },
  { 0, 0, 1, - 7, 124,  13, - 4, 1 }, { 0, 0, 2, - 8, 123,  15, - 5, 1 },
  { 0, 0, 2, - 9, 122,  18, - 6, 1 }, { 0, 0, 2, -10, 121,  20, - 6, 1 },
  { 0, 0, 2, -11, 120,  22, - 7, 2 }, { 0, 0, 2, -12, 119,  25, - 8, 2 },
  { 0, 0, 3, -13, 117,  27, - 8, 2 }, { 0, 0, 3, -13, 116,  29, - 9, 2 },
  { 0, 0, 3, -14, 114,  32, -10, 3 }, { 0, 0, 3, -15, 113,  35, -10, 2 },
  { 0, 0, 3, -15, 111,  37, -11, 3 }, { 0, 0, 3, -16, 109,  40, -11, 3 },
  { 0, 0, 3, -16, 108,  42, -12, 3 }, { 0, 0, 4, -17, 106,  45, -13, 3 },
  { 0, 0, 4, -17, 104,  47, -13, 3 }, { 0, 0, 4, -17, 102,  50, -14, 3 },
  { 0, 0, 4, -17, 100,  52, -14, 3 }, { 0, 0, 4, -18,  98,  55, -15, 4 },
  { 0, 0, 4, -18,  96,  58, -15, 3 }, { 0, 0, 4, -18,  94,  60, -16, 4 },
  { 0, 0, 4, -18,  91,  63, -16, 4 }, { 0, 0, 4, -18,  89,  65, -16, 4 },
  { 0, 0, 4, -18,  87,  68, -17, 4 }, { 0, 0, 4, -18,  85,  70, -17, 4 },
  { 0, 0, 4, -18,  82,  73, -17, 4 }, { 0, 0, 4, -18,  80,  75, -17, 4 },
  { 0, 0, 4, -18,  78,  78, -18, 4 }, { 0, 0, 4, -17,  75,  80, -18, 4 },
  { 0, 0, 4, -17,  73,  82, -18, 4 }, { 0, 0, 4, -17,  70,  85, -18, 4 },
  { 0, 0, 4, -17,  68,  87, -18, 4 }, { 0, 0, 4, -16,  65,  89, -18, 4 },
  { 0, 0, 4, -16,  63,  91, -18, 4 }, { 0, 0, 4, -16,  60,  94, -18, 4 },
  { 0, 0, 3, -15,  58,  96, -18, 4 }, { 0, 0, 4, -15,  55,  98, -18, 4 },
  { 0, 0, 3, -14,  52, 100, -17, 4 }, { 0, 0, 3, -14,  50, 102, -17, 4 },
  { 0, 0, 3, -13,  47, 104, -17, 4 }, { 0, 0, 3, -13,  45, 106, -17, 4 },
  { 0, 0, 3, -12,  42, 108, -16, 3 }, { 0, 0, 3, -11,  40, 109, -16, 3 },
  { 0, 0, 3, -11,  37, 111, -15, 3 }, { 0, 0, 2, -10,  35, 113, -15, 3 },
  { 0, 0, 3, -10,  32, 114, -14, 3 }, { 0, 0, 2, - 9,  29, 116, -13, 3 },
  { 0, 0, 2, - 8,  27, 117, -13, 3 }, { 0, 0, 2, - 8,  25, 119, -12, 2 },
  { 0, 0, 2, - 7,  22, 120, -11, 2 }, { 0, 0, 1, - 6,  20, 121, -10, 2 },
  { 0, 0, 1, - 6,  18, 122, - 9, 2 }, { 0, 0, 1, - 5,  15, 123, - 8, 2 },
  { 0, 0, 1, - 4,  13, 124, - 7, 1 }, { 0, 0, 1, - 4,  11, 125, - 6, 1 },
  { 0, 0, 1, - 3,   8, 126, - 5, 1 }, { 0, 0, 1, - 2,   6, 126, - 4, 1 },
  { 0, 0, 0, - 1,   4, 127, - 3, 1 }, { 0, 0, 0, - 1,   2, 128, - 1, 0 },

  // dummy
  { 0, 0, 0, 0,   0, 128, 0, 0 },
};
/* clang-format on */

#define DIV_LUT_PREC_BITS 14
#define DIV_LUT_BITS 8
#define DIV_LUT_NUM (1 << DIV_LUT_BITS)

static const uint16_t div_lut[DIV_LUT_NUM + 1] = {
  16384, 16320, 16257, 16194, 16132, 16070, 16009, 15948, 15888, 15828, 15768,
  15709, 15650, 15592, 15534, 15477, 15420, 15364, 15308, 15252, 15197, 15142,
  15087, 15033, 14980, 14926, 14873, 14821, 14769, 14717, 14665, 14614, 14564,
  14513, 14463, 14413, 14364, 14315, 14266, 14218, 14170, 14122, 14075, 14028,
  13981, 13935, 13888, 13843, 13797, 13752, 13707, 13662, 13618, 13574, 13530,
  13487, 13443, 13400, 13358, 13315, 13273, 13231, 13190, 13148, 13107, 13066,
  13026, 12985, 12945, 12906, 12866, 12827, 12788, 12749, 12710, 12672, 12633,
  12596, 12558, 12520, 12483, 12446, 12409, 12373, 12336, 12300, 12264, 12228,
  12193, 12157, 12122, 12087, 12053, 12018, 11984, 11950, 11916, 11882, 11848,
  11815, 11782, 11749, 11716, 11683, 11651, 11619, 11586, 11555, 11523, 11491,
  11460, 11429, 11398, 11367, 11336, 11305, 11275, 11245, 11215, 11185, 11155,
  11125, 11096, 11067, 11038, 11009, 10980, 10951, 10923, 10894, 10866, 10838,
  10810, 10782, 10755, 10727, 10700, 10673, 10645, 10618, 10592, 10565, 10538,
  10512, 10486, 10460, 10434, 10408, 10382, 10356, 10331, 10305, 10280, 10255,
  10230, 10205, 10180, 10156, 10131, 10107, 10082, 10058, 10034, 10010, 9986,
  9963,  9939,  9916,  9892,  9869,  9846,  9823,  9800,  9777,  9754,  9732,
  9709,  9687,  9664,  9642,  9620,  9598,  9576,  9554,  9533,  9511,  9489,
  9468,  9447,  9425,  9404,  9383,  9362,  9341,  9321,  9300,  9279,  9259,
  9239,  9218,  9198,  9178,  9158,  9138,  9118,  9098,  9079,  9059,  9039,
  9020,  9001,  8981,  8962,  8943,  8924,  8905,  8886,  8867,  8849,  8830,
  8812,  8793,  8775,  8756,  8738,  8720,  8702,  8684,  8666,  8648,  8630,
  8613,  8595,  8577,  8560,  8542,  8525,  8508,  8490,  8473,  8456,  8439,
  8422,  8405,  8389,  8372,  8355,  8339,  8322,  8306,  8289,  8273,  8257,
  8240,  8224,  8208,  8192,
};

static inline int16_t saturate_int16(int32_t v) {
  if (v > 32767)
    return 32767;
  else if (v < -32768)
    return -32768;
  return v;
}

#if CONFIG_WARPED_MOTION
// Decomposes a divisor D such that 1/D = y/2^shift, where y is returned
// at precision of DIV_LUT_PREC_BITS along with the shift.
static int16_t resolve_divisor_64(uint64_t D, int16_t *shift) {
  int64_t e, f;
  *shift = (D >> 32) ? get_msb(D >> 32) + 32 : get_msb(D);
  // e is obtained from D after resetting the most significant 1 bit.
  e = D - ((uint64_t)1 << *shift);
  // Get the most significant DIV_LUT_BITS (8) bits of e into f
  if (*shift > DIV_LUT_BITS)
    f = ROUND_POWER_OF_TWO_64(e, *shift - DIV_LUT_BITS);
  else
    f = e << (DIV_LUT_BITS - *shift);
  assert(f <= DIV_LUT_NUM);
  *shift += DIV_LUT_PREC_BITS;
  // Use f as lookup into the precomputed table of multipliers
  return div_lut[f];
}
#endif  // CONFIG_WARPED_MOTION

static int16_t resolve_divisor_32(uint32_t D, int16_t *shift) {
  int32_t e, f;
  *shift = get_msb(D);
  // e is obtained from D after resetting the most significant 1 bit.
  e = D - ((uint32_t)1 << *shift);
  // Get the most significant DIV_LUT_BITS (8) bits of e into f
  if (*shift > DIV_LUT_BITS)
    f = ROUND_POWER_OF_TWO(e, *shift - DIV_LUT_BITS);
  else
    f = e << (DIV_LUT_BITS - *shift);
  assert(f <= DIV_LUT_NUM);
  *shift += DIV_LUT_PREC_BITS;
  // Use f as lookup into the precomputed table of multipliers
  return div_lut[f];
}

static int is_affine_valid(WarpedMotionParams *wm) {
  const int32_t *mat = wm->wmmat;
  return (mat[2] > 0);
}

static int is_affine_shear_allowed(int32_t alpha, int32_t beta, int32_t gamma,
                                   int32_t delta) {
  if ((4 * abs(alpha) + 7 * abs(beta) > (1 << WARPEDMODEL_PREC_BITS)) ||
      (4 * abs(gamma) + 4 * abs(delta) > (1 << WARPEDMODEL_PREC_BITS)))
    return 0;
  else
    return 1;
}

// Returns 1 on success or 0 on an invalid affine set
int get_shear_params(WarpedMotionParams *wm) {
  const int32_t *mat = wm->wmmat;
  if (!is_affine_valid(wm)) return 0;
  wm->alpha = mat[2] - (1 << WARPEDMODEL_PREC_BITS);
  wm->beta = mat[3];
  int16_t shift;
  int16_t y = resolve_divisor_32(abs(mat[2]), &shift) * (mat[2] < 0 ? -1 : 1);
  int64_t v;
  v = ((int64_t)mat[4] << WARPEDMODEL_PREC_BITS) * y;
  wm->gamma = ROUND_POWER_OF_TWO_SIGNED_64(v, shift);
  v = ((int64_t)mat[3] * mat[4]) * y;
  wm->delta = mat[5] - ROUND_POWER_OF_TWO_SIGNED_64(v, shift) -
              (1 << WARPEDMODEL_PREC_BITS);
  if (!is_affine_shear_allowed(wm->alpha, wm->beta, wm->gamma, wm->delta))
    return 0;
  return 1;
}

#if CONFIG_AOM_HIGHBITDEPTH
static INLINE void highbd_get_subcolumn(int taps, uint16_t *ref, int32_t *col,
                                        int stride, int x, int y_start) {
  int i;
  for (i = 0; i < taps; ++i) {
    col[i] = ref[(i + y_start) * stride + x];
  }
}

static uint16_t highbd_bi_ntap_filter(uint16_t *ref, int x, int y, int stride,
                                      int bd) {
  int32_t val, arr[WARPEDPIXEL_FILTER_TAPS];
  int k;
  int i = (int)x >> WARPEDPIXEL_PREC_BITS;
  int j = (int)y >> WARPEDPIXEL_PREC_BITS;
  for (k = 0; k < WARPEDPIXEL_FILTER_TAPS; ++k) {
    int32_t arr_temp[WARPEDPIXEL_FILTER_TAPS];
    highbd_get_subcolumn(WARPEDPIXEL_FILTER_TAPS, ref, arr_temp, stride,
                         i + k + 1 - WARPEDPIXEL_FILTER_TAPS / 2,
                         j + 1 - WARPEDPIXEL_FILTER_TAPS / 2);
    arr[k] = do_ntap_filter(arr_temp + WARPEDPIXEL_FILTER_TAPS / 2 - 1,
                            y - (j * (1 << WARPEDPIXEL_PREC_BITS)));
  }
  val = do_ntap_filter(arr + WARPEDPIXEL_FILTER_TAPS / 2 - 1,
                       x - (i * (1 << WARPEDPIXEL_PREC_BITS)));
  val = ROUND_POWER_OF_TWO_SIGNED(val, WARPEDPIXEL_FILTER_BITS * 2);
  return (uint16_t)clip_pixel_highbd(val, bd);
}

static uint16_t highbd_bi_cubic_filter(uint16_t *ref, int x, int y, int stride,
                                       int bd) {
  int32_t val, arr[4];
  int k;
  int i = (int)x >> WARPEDPIXEL_PREC_BITS;
  int j = (int)y >> WARPEDPIXEL_PREC_BITS;
  for (k = 0; k < 4; ++k) {
    int32_t arr_temp[4];
    highbd_get_subcolumn(4, ref, arr_temp, stride, i + k - 1, j - 1);
    arr[k] =
        do_cubic_filter(arr_temp + 1, y - (j * (1 << WARPEDPIXEL_PREC_BITS)));
  }
  val = do_cubic_filter(arr + 1, x - (i * (1 << WARPEDPIXEL_PREC_BITS)));
  val = ROUND_POWER_OF_TWO_SIGNED(val, WARPEDPIXEL_FILTER_BITS * 2);
  return (uint16_t)clip_pixel_highbd(val, bd);
}

static uint16_t highbd_bi_linear_filter(uint16_t *ref, int x, int y, int stride,
                                        int bd) {
  const int ix = x >> WARPEDPIXEL_PREC_BITS;
  const int iy = y >> WARPEDPIXEL_PREC_BITS;
  const int sx = x - (ix * (1 << WARPEDPIXEL_PREC_BITS));
  const int sy = y - (iy * (1 << WARPEDPIXEL_PREC_BITS));
  int32_t val;
  val = ROUND_POWER_OF_TWO_SIGNED(
      ref[iy * stride + ix] * (WARPEDPIXEL_PREC_SHIFTS - sy) *
              (WARPEDPIXEL_PREC_SHIFTS - sx) +
          ref[iy * stride + ix + 1] * (WARPEDPIXEL_PREC_SHIFTS - sy) * sx +
          ref[(iy + 1) * stride + ix] * sy * (WARPEDPIXEL_PREC_SHIFTS - sx) +
          ref[(iy + 1) * stride + ix + 1] * sy * sx,
      WARPEDPIXEL_PREC_BITS * 2);
  return (uint16_t)clip_pixel_highbd(val, bd);
}

static uint16_t highbd_warp_interpolate(uint16_t *ref, int x, int y, int width,
                                        int height, int stride, int bd) {
  int ix = x >> WARPEDPIXEL_PREC_BITS;
  int iy = y >> WARPEDPIXEL_PREC_BITS;
  int sx = x - (ix * (1 << WARPEDPIXEL_PREC_BITS));
  int sy = y - (iy * (1 << WARPEDPIXEL_PREC_BITS));
  int32_t v;

  if (ix < 0 && iy < 0)
    return ref[0];
  else if (ix < 0 && iy > height - 1)
    return ref[(height - 1) * stride];
  else if (ix > width - 1 && iy < 0)
    return ref[width - 1];
  else if (ix > width - 1 && iy > height - 1)
    return ref[(height - 1) * stride + (width - 1)];
  else if (ix < 0) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[iy * stride] * (WARPEDPIXEL_PREC_SHIFTS - sy) +
            ref[(iy + 1) * stride] * sy,
        WARPEDPIXEL_PREC_BITS);
    return clip_pixel_highbd(v, bd);
  } else if (iy < 0) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[ix] * (WARPEDPIXEL_PREC_SHIFTS - sx) + ref[ix + 1] * sx,
        WARPEDPIXEL_PREC_BITS);
    return clip_pixel_highbd(v, bd);
  } else if (ix > width - 1) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[iy * stride + width - 1] * (WARPEDPIXEL_PREC_SHIFTS - sy) +
            ref[(iy + 1) * stride + width - 1] * sy,
        WARPEDPIXEL_PREC_BITS);
    return clip_pixel_highbd(v, bd);
  } else if (iy > height - 1) {
    v = ROUND_POWER_OF_TWO_SIGNED(
        ref[(height - 1) * stride + ix] * (WARPEDPIXEL_PREC_SHIFTS - sx) +
            ref[(height - 1) * stride + ix + 1] * sx,
        WARPEDPIXEL_PREC_BITS);
    return clip_pixel_highbd(v, bd);
  } else if (ix >= WARPEDPIXEL_FILTER_TAPS / 2 - 1 &&
             iy >= WARPEDPIXEL_FILTER_TAPS / 2 - 1 &&
             ix < width - WARPEDPIXEL_FILTER_TAPS / 2 &&
             iy < height - WARPEDPIXEL_FILTER_TAPS / 2) {
    return highbd_bi_ntap_filter(ref, x, y, stride, bd);
  } else if (ix >= 1 && iy >= 1 && ix < width - 2 && iy < height - 2) {
    return highbd_bi_cubic_filter(ref, x, y, stride, bd);
  } else {
    return highbd_bi_linear_filter(ref, x, y, stride, bd);
  }
}

static INLINE int highbd_error_measure(int err, int bd) {
  const int b = bd - 8;
  const int bmask = (1 << b) - 1;
  const int v = (1 << b);
  int e1, e2;
  err = abs(err);
  e1 = err >> b;
  e2 = err & bmask;
  return error_measure_lut[255 + e1] * (v - e2) +
         error_measure_lut[256 + e1] * e2;
}

static void highbd_warp_plane_old(WarpedMotionParams *wm, uint8_t *ref8,
                                  int width, int height, int stride,
                                  uint8_t *pred8, int p_col, int p_row,
                                  int p_width, int p_height, int p_stride,
                                  int subsampling_x, int subsampling_y,
                                  int x_scale, int y_scale, int bd,
                                  int ref_frm) {
  int i, j;
  ProjectPointsFunc projectpoints = get_project_points_type(wm->wmtype);
  uint16_t *pred = CONVERT_TO_SHORTPTR(pred8);
  uint16_t *ref = CONVERT_TO_SHORTPTR(ref8);
  if (projectpoints == NULL) return;
  for (i = p_row; i < p_row + p_height; ++i) {
    for (j = p_col; j < p_col + p_width; ++j) {
      int in[2], out[2];
      in[0] = j;
      in[1] = i;
      projectpoints(wm->wmmat, in, out, 1, 2, 2, subsampling_x, subsampling_y);
      out[0] = ROUND_POWER_OF_TWO_SIGNED(out[0] * x_scale, 4);
      out[1] = ROUND_POWER_OF_TWO_SIGNED(out[1] * y_scale, 4);
      if (ref_frm)
        pred[(j - p_col) + (i - p_row) * p_stride] = ROUND_POWER_OF_TWO(
            pred[(j - p_col) + (i - p_row) * p_stride] +
                highbd_warp_interpolate(ref, out[0], out[1], width, height,
                                        stride, bd),
            1);
      else
        pred[(j - p_col) + (i - p_row) * p_stride] = highbd_warp_interpolate(
            ref, out[0], out[1], width, height, stride, bd);
    }
  }
}

// Note: For an explanation of the warp algorithm, see the comment
// above warp_plane()
//
// Note also: The "worst case" in terms of modulus of the data stored into 'tmp'
// (ie, the result of 'sum' in the horizontal filter) occurs when:
// coeffs = { -2,   8, -22,  87,  72, -21,   8, -2}, and
// ref =    {  0, 255,   0, 255, 255,   0, 255,  0}
// Before rounding, this gives sum = 716625. After rounding,
// HORSHEAR_REDUCE_PREC_BITS = 4 => sum = 44789 > 2^15
// HORSHEAR_REDUCE_PREC_BITS = 5 => sum = 22395 < 2^15
//
// So, as long as HORSHEAR_REDUCE_PREC_BITS >= 5, we can safely use a 16-bit
// intermediate array.
void av1_highbd_warp_affine_c(int32_t *mat, uint16_t *ref, int width,
                              int height, int stride, uint16_t *pred, int p_col,
                              int p_row, int p_width, int p_height,
                              int p_stride, int subsampling_x,
                              int subsampling_y, int bd, int ref_frm,
                              int32_t alpha, int32_t beta, int32_t gamma,
                              int32_t delta) {
#if HORSHEAR_REDUCE_PREC_BITS >= 5
  int16_t tmp[15 * 8];
#else
  int32_t tmp[15 * 8];
#endif
  int i, j, k, l, m;

  /* Note: For this code to work, the left/right frame borders need to be
     extended by at least 13 pixels each. By the time we get here, other
     code will have set up this border, but we allow an explicit check
     for debugging purposes.
  */
  /*for (i = 0; i < height; ++i) {
    for (j = 0; j < 13; ++j) {
      assert(ref[i * stride - 13 + j] == ref[i * stride]);
      assert(ref[i * stride + width + j] == ref[i * stride + (width - 1)]);
    }
  }*/

  for (i = p_row; i < p_row + p_height; i += 8) {
    for (j = p_col; j < p_col + p_width; j += 8) {
      int32_t x4, y4, ix4, sx4, iy4, sy4;
      if (subsampling_x)
        x4 = ROUND_POWER_OF_TWO_SIGNED(
            mat[2] * 2 * (j + 4) + mat[3] * 2 * (i + 4) + mat[0] +
                (mat[2] + mat[3] - (1 << WARPEDMODEL_PREC_BITS)) / 2,
            1);
      else
        x4 = mat[2] * (j + 4) + mat[3] * (i + 4) + mat[0];

      if (subsampling_y)
        y4 = ROUND_POWER_OF_TWO_SIGNED(
            mat[4] * 2 * (j + 4) + mat[5] * 2 * (i + 4) + mat[1] +
                (mat[4] + mat[5] - (1 << WARPEDMODEL_PREC_BITS)) / 2,
            1);
      else
        y4 = mat[4] * (j + 4) + mat[5] * (i + 4) + mat[1];

      ix4 = x4 >> WARPEDMODEL_PREC_BITS;
      sx4 = x4 & ((1 << WARPEDMODEL_PREC_BITS) - 1);
      iy4 = y4 >> WARPEDMODEL_PREC_BITS;
      sy4 = y4 & ((1 << WARPEDMODEL_PREC_BITS) - 1);

      // Horizontal filter
      for (k = -7; k < 8; ++k) {
        int iy = iy4 + k;
        if (iy < 0)
          iy = 0;
        else if (iy > height - 1)
          iy = height - 1;

        if (ix4 <= -7) {
          for (l = 0; l < 8; ++l) {
            tmp[(k + 7) * 8 + l] =
                ref[iy * stride] *
                (1 << (WARPEDPIXEL_FILTER_BITS - HORSHEAR_REDUCE_PREC_BITS));
          }
        } else if (ix4 >= width + 6) {
          for (l = 0; l < 8; ++l) {
            tmp[(k + 7) * 8 + l] =
                ref[iy * stride + (width - 1)] *
                (1 << (WARPEDPIXEL_FILTER_BITS - HORSHEAR_REDUCE_PREC_BITS));
          }
        } else {
          int sx = sx4 + alpha * (-4) + beta * k;

          for (l = -4; l < 4; ++l) {
            int ix = ix4 + l - 3;
            const int offs = ROUND_POWER_OF_TWO(sx, WARPEDDIFF_PREC_BITS) +
                             WARPEDPIXEL_PREC_SHIFTS;
            const int16_t *coeffs = warped_filter[offs];
            int32_t sum = 0;
            // assert(offs >= 0 && offs <= WARPEDPIXEL_PREC_SHIFTS * 3);
            for (m = 0; m < 8; ++m) {
              sum += ref[iy * stride + ix + m] * coeffs[m];
            }
            sum = ROUND_POWER_OF_TWO(sum, HORSHEAR_REDUCE_PREC_BITS);
#if HORSHEAR_REDUCE_PREC_BITS >= 5
            tmp[(k + 7) * 8 + (l + 4)] = saturate_int16(sum);
#else
            tmp[(k + 7) * 8 + (l + 4)] = sum;
#endif
            sx += alpha;
          }
        }
      }

      // Vertical filter
      for (k = -4; k < AOMMIN(4, p_row + p_height - i - 4); ++k) {
        int sy = sy4 + gamma * (-4) + delta * k;
        for (l = -4; l < 4; ++l) {
          uint16_t *p =
              &pred[(i - p_row + k + 4) * p_stride + (j - p_col + l + 4)];
          const int offs = ROUND_POWER_OF_TWO(sy, WARPEDDIFF_PREC_BITS) +
                           WARPEDPIXEL_PREC_SHIFTS;
          const int16_t *coeffs = warped_filter[offs];
          int32_t sum = 0;
          // assert(offs >= 0 && offs <= WARPEDPIXEL_PREC_SHIFTS * 3);
          for (m = 0; m < 8; ++m) {
            sum += tmp[(k + m + 4) * 8 + (l + 4)] * coeffs[m];
          }
          sum = clip_pixel_highbd(
              ROUND_POWER_OF_TWO(sum, VERSHEAR_REDUCE_PREC_BITS), bd);
          if (ref_frm)
            *p = ROUND_POWER_OF_TWO(*p + sum, 1);
          else
            *p = sum;
          sy += gamma;
        }
      }
    }
  }
}

static void highbd_warp_plane(WarpedMotionParams *wm, uint8_t *ref8, int width,
                              int height, int stride, uint8_t *pred8, int p_col,
                              int p_row, int p_width, int p_height,
                              int p_stride, int subsampling_x,
                              int subsampling_y, int x_scale, int y_scale,
                              int bd, int ref_frm) {
  if (wm->wmtype == ROTZOOM) {
    wm->wmmat[5] = wm->wmmat[2];
    wm->wmmat[4] = -wm->wmmat[3];
  }
  if ((wm->wmtype == ROTZOOM || wm->wmtype == AFFINE) && x_scale == 16 &&
      y_scale == 16) {
    int32_t *mat = wm->wmmat;
    const int32_t alpha = wm->alpha;
    const int32_t beta = wm->beta;
    const int32_t gamma = wm->gamma;
    const int32_t delta = wm->delta;

    uint16_t *ref = CONVERT_TO_SHORTPTR(ref8);
    uint16_t *pred = CONVERT_TO_SHORTPTR(pred8);
    av1_highbd_warp_affine(mat, ref, width, height, stride, pred, p_col, p_row,
                           p_width, p_height, p_stride, subsampling_x,
                           subsampling_y, bd, ref_frm, alpha, beta, gamma,
                           delta);
  } else {
    highbd_warp_plane_old(wm, ref8, width, height, stride, pred8, p_col, p_row,
                          p_width, p_height, p_stride, subsampling_x,
                          subsampling_y, x_scale, y_scale, bd, ref_frm);
  }
}

static double highbd_warp_erroradv(WarpedMotionParams *wm, uint8_t *ref8,
                                   int width, int height, int stride,
                                   uint8_t *dst8, int p_col, int p_row,
                                   int p_width, int p_height, int p_stride,
                                   int subsampling_x, int subsampling_y,
                                   int x_scale, int y_scale, int bd) {
  int gm_err = 0, no_gm_err = 0;
  int64_t gm_sumerr = 0, no_gm_sumerr = 0;
  int i, j;
  uint16_t *tmp = aom_malloc(p_width * p_height * sizeof(*tmp));
  uint16_t *dst = CONVERT_TO_SHORTPTR(dst8);
  uint16_t *ref = CONVERT_TO_SHORTPTR(ref8);
  highbd_warp_plane(wm, ref8, width, height, stride, CONVERT_TO_BYTEPTR(tmp),
                    p_col, p_row, p_width, p_height, p_width, subsampling_x,
                    subsampling_y, x_scale, y_scale, bd, 0);
  for (i = 0; i < p_height; ++i) {
    for (j = 0; j < p_width; ++j) {
      gm_err = dst[j + i * p_stride] - tmp[j + i * p_width];
      no_gm_err =
          dst[j + i * p_stride] - ref[(j + p_col) + (i + p_row) * stride];
      gm_sumerr += highbd_error_measure(gm_err, bd);
      no_gm_sumerr += highbd_error_measure(no_gm_err, bd);
    }
  }
  aom_free(tmp);
  return (double)gm_sumerr / no_gm_sumerr;
}
#endif  // CONFIG_AOM_HIGHBITDEPTH

static INLINE int error_measure(int err) {
  return error_measure_lut[255 + err];
}

static void warp_plane_old(WarpedMotionParams *wm, uint8_t *ref, int width,
                           int height, int stride, uint8_t *pred, int p_col,
                           int p_row, int p_width, int p_height, int p_stride,
                           int subsampling_x, int subsampling_y, int x_scale,
                           int y_scale, int ref_frm) {
  int i, j;
  ProjectPointsFunc projectpoints = get_project_points_type(wm->wmtype);
  if (projectpoints == NULL) return;
  for (i = p_row; i < p_row + p_height; ++i) {
    for (j = p_col; j < p_col + p_width; ++j) {
      int in[2], out[2];
      in[0] = j;
      in[1] = i;
      projectpoints(wm->wmmat, in, out, 1, 2, 2, subsampling_x, subsampling_y);
      out[0] = ROUND_POWER_OF_TWO_SIGNED(out[0] * x_scale, 4);
      out[1] = ROUND_POWER_OF_TWO_SIGNED(out[1] * y_scale, 4);
      if (ref_frm)
        pred[(j - p_col) + (i - p_row) * p_stride] = ROUND_POWER_OF_TWO(
            pred[(j - p_col) + (i - p_row) * p_stride] +
                warp_interpolate(ref, out[0], out[1], width, height, stride),
            1);
      else
        pred[(j - p_col) + (i - p_row) * p_stride] =
            warp_interpolate(ref, out[0], out[1], width, height, stride);
    }
  }
}

/* The warp filter for ROTZOOM and AFFINE models works as follows:
   * Split the input into 8x8 blocks
   * For each block, project the point (4, 4) within the block, to get the
     overall block position. Split into integer and fractional coordinates,
     maintaining full WARPEDMODEL precision
   * Filter horizontally: Generate 15 rows of 8 pixels each. Each pixel gets a
     variable horizontal offset. This means that, while the rows of the
     intermediate buffer align with the rows of the *reference* image, the
     columns align with the columns of the *destination* image.
   * Filter vertically: Generate the output block (up to 8x8 pixels, but if the
     destination is too small we crop the output at this stage). Each pixel has
     a variable vertical offset, so that the resulting rows are aligned with
     the rows of the destination image.

   To accomplish these alignments, we factor the warp matrix as a
   product of two shear / asymmetric zoom matrices:
   / a b \  = /   1       0    \ * / 1+alpha  beta \
   \ c d /    \ gamma  1+delta /   \    0      1   /
   where a, b, c, d are wmmat[2], wmmat[3], wmmat[4], wmmat[5] respectively.
   The second shear (with alpha and beta) is applied by the horizontal filter,
   then the first shear (with gamma and delta) is applied by the vertical
   filter.

   The only limitation is that, to fit this in a fixed 8-tap filter size,
   the fractional pixel offsets must be at most +-1. Since the horizontal filter
   generates 15 rows of 8 columns, and the initial point we project is at (4, 4)
   within the block, the parameters must satisfy
   4 * |alpha| + 7 * |beta| <= 1   and   4 * |gamma| + 7 * |delta| <= 1
   for this filter to be applicable.

   Note: warp_affine() assumes that the caller has done all of the relevant
   checks, ie. that we have a ROTZOOM or AFFINE model, that wm[4] and wm[5]
   are set appropriately (if using a ROTZOOM model), and that alpha, beta,
   gamma, delta are all in range.

   TODO(david.barker): Maybe support scaled references?
*/
void av1_warp_affine_c(int32_t *mat, uint8_t *ref, int width, int height,
                       int stride, uint8_t *pred, int p_col, int p_row,
                       int p_width, int p_height, int p_stride,
                       int subsampling_x, int subsampling_y, int ref_frm,
                       int32_t alpha, int32_t beta, int32_t gamma,
                       int32_t delta) {
  int16_t tmp[15 * 8];
  int i, j, k, l, m;

  /* Note: For this code to work, the left/right frame borders need to be
     extended by at least 13 pixels each. By the time we get here, other
     code will have set up this border, but we allow an explicit check
     for debugging purposes.
  */
  /*for (i = 0; i < height; ++i) {
    for (j = 0; j < 13; ++j) {
      assert(ref[i * stride - 13 + j] == ref[i * stride]);
      assert(ref[i * stride + width + j] == ref[i * stride + (width - 1)]);
    }
  }*/

  for (i = p_row; i < p_row + p_height; i += 8) {
    for (j = p_col; j < p_col + p_width; j += 8) {
      int32_t x4, y4, ix4, sx4, iy4, sy4;
      if (subsampling_x)
        x4 = ROUND_POWER_OF_TWO_SIGNED(
            mat[2] * 2 * (j + 4) + mat[3] * 2 * (i + 4) + mat[0] +
                (mat[2] + mat[3] - (1 << WARPEDMODEL_PREC_BITS)) / 2,
            1);
      else
        x4 = mat[2] * (j + 4) + mat[3] * (i + 4) + mat[0];

      if (subsampling_y)
        y4 = ROUND_POWER_OF_TWO_SIGNED(
            mat[4] * 2 * (j + 4) + mat[5] * 2 * (i + 4) + mat[1] +
                (mat[4] + mat[5] - (1 << WARPEDMODEL_PREC_BITS)) / 2,
            1);
      else
        y4 = mat[4] * (j + 4) + mat[5] * (i + 4) + mat[1];

      ix4 = x4 >> WARPEDMODEL_PREC_BITS;
      sx4 = x4 & ((1 << WARPEDMODEL_PREC_BITS) - 1);
      iy4 = y4 >> WARPEDMODEL_PREC_BITS;
      sy4 = y4 & ((1 << WARPEDMODEL_PREC_BITS) - 1);

      // Horizontal filter
      for (k = -7; k < 8; ++k) {
        int iy = iy4 + k;
        if (iy < 0)
          iy = 0;
        else if (iy > height - 1)
          iy = height - 1;

        if (ix4 <= -7) {
          // In this case, the rightmost pixel sampled is in column
          // ix4 + 3 + 7 - 3 = ix4 + 7 <= 0, ie. the entire block
          // will sample only from the leftmost column
          // (once border extension is taken into account)
          for (l = 0; l < 8; ++l) {
            tmp[(k + 7) * 8 + l] =
                ref[iy * stride] *
                (1 << (WARPEDPIXEL_FILTER_BITS - HORSHEAR_REDUCE_PREC_BITS));
          }
        } else if (ix4 >= width + 6) {
          // In this case, the leftmost pixel sampled is in column
          // ix4 - 4 + 0 - 3 = ix4 - 7 >= width - 1, ie. the entire block
          // will sample only from the rightmost column
          // (once border extension is taken into account)
          for (l = 0; l < 8; ++l) {
            tmp[(k + 7) * 8 + l] =
                ref[iy * stride + (width - 1)] *
                (1 << (WARPEDPIXEL_FILTER_BITS - HORSHEAR_REDUCE_PREC_BITS));
          }
        } else {
          // If we get here, then
          // the leftmost pixel sampled is
          // ix4 - 4 + 0 - 3 = ix4 - 7 >= -13
          // and the rightmost pixel sampled is at most
          // ix4 + 3 + 7 - 3 = ix4 + 7 <= width + 12
          // So, assuming that border extension has been done, we
          // don't need to explicitly clamp values.
          int sx = sx4 + alpha * (-4) + beta * k;

          for (l = -4; l < 4; ++l) {
            int ix = ix4 + l - 3;
            // At this point, sx = sx4 + alpha * l + beta * k
            const int offs = ROUND_POWER_OF_TWO(sx, WARPEDDIFF_PREC_BITS) +
                             WARPEDPIXEL_PREC_SHIFTS;
            const int16_t *coeffs = warped_filter[offs];
            int32_t sum = 0;
            // assert(offs >= 0 && offs <= WARPEDPIXEL_PREC_SHIFTS * 3);
            for (m = 0; m < 8; ++m) {
              sum += ref[iy * stride + ix + m] * coeffs[m];
            }
            sum = ROUND_POWER_OF_TWO(sum, HORSHEAR_REDUCE_PREC_BITS);
            tmp[(k + 7) * 8 + (l + 4)] = saturate_int16(sum);
            sx += alpha;
          }
        }
      }

      // Vertical filter
      for (k = -4; k < AOMMIN(4, p_row + p_height - i - 4); ++k) {
        int sy = sy4 + gamma * (-4) + delta * k;
        for (l = -4; l < 4; ++l) {
          uint8_t *p =
              &pred[(i - p_row + k + 4) * p_stride + (j - p_col + l + 4)];
          // At this point, sy = sy4 + gamma * l + delta * k
          const int offs = ROUND_POWER_OF_TWO(sy, WARPEDDIFF_PREC_BITS) +
                           WARPEDPIXEL_PREC_SHIFTS;
          const int16_t *coeffs = warped_filter[offs];
          int32_t sum = 0;
          // assert(offs >= 0 && offs <= WARPEDPIXEL_PREC_SHIFTS * 3);
          for (m = 0; m < 8; ++m) {
            sum += tmp[(k + m + 4) * 8 + (l + 4)] * coeffs[m];
          }
          sum = clip_pixel(ROUND_POWER_OF_TWO(sum, VERSHEAR_REDUCE_PREC_BITS));
          if (ref_frm)
            *p = ROUND_POWER_OF_TWO(*p + sum, 1);
          else
            *p = sum;
          sy += gamma;
        }
      }
    }
  }
}

static void warp_plane(WarpedMotionParams *wm, uint8_t *ref, int width,
                       int height, int stride, uint8_t *pred, int p_col,
                       int p_row, int p_width, int p_height, int p_stride,
                       int subsampling_x, int subsampling_y, int x_scale,
                       int y_scale, int ref_frm) {
  if (wm->wmtype == ROTZOOM) {
    wm->wmmat[5] = wm->wmmat[2];
    wm->wmmat[4] = -wm->wmmat[3];
  }
  if ((wm->wmtype == ROTZOOM || wm->wmtype == AFFINE) && x_scale == 16 &&
      y_scale == 16) {
    int32_t *mat = wm->wmmat;
    const int32_t alpha = wm->alpha;
    const int32_t beta = wm->beta;
    const int32_t gamma = wm->gamma;
    const int32_t delta = wm->delta;

    av1_warp_affine(mat, ref, width, height, stride, pred, p_col, p_row,
                    p_width, p_height, p_stride, subsampling_x, subsampling_y,
                    ref_frm, alpha, beta, gamma, delta);
  } else {
    warp_plane_old(wm, ref, width, height, stride, pred, p_col, p_row, p_width,
                   p_height, p_stride, subsampling_x, subsampling_y, x_scale,
                   y_scale, ref_frm);
  }
}

static double warp_erroradv(WarpedMotionParams *wm, uint8_t *ref, int width,
                            int height, int stride, uint8_t *dst, int p_col,
                            int p_row, int p_width, int p_height, int p_stride,
                            int subsampling_x, int subsampling_y, int x_scale,
                            int y_scale) {
  int gm_err = 0, no_gm_err = 0;
  int gm_sumerr = 0, no_gm_sumerr = 0;
  int i, j;
  uint8_t *tmp = aom_malloc(p_width * p_height);
  warp_plane(wm, ref, width, height, stride, tmp, p_col, p_row, p_width,
             p_height, p_width, subsampling_x, subsampling_y, x_scale, y_scale,
             0);

  for (i = 0; i < p_height; ++i) {
    for (j = 0; j < p_width; ++j) {
      gm_err = dst[j + i * p_stride] - tmp[j + i * p_width];
      no_gm_err =
          dst[j + i * p_stride] - ref[(j + p_col) + (i + p_row) * stride];
      gm_sumerr += error_measure(gm_err);
      no_gm_sumerr += error_measure(no_gm_err);
    }
  }

  aom_free(tmp);
  return (double)gm_sumerr / no_gm_sumerr;
}

double av1_warp_erroradv(WarpedMotionParams *wm,
#if CONFIG_AOM_HIGHBITDEPTH
                         int use_hbd, int bd,
#endif  // CONFIG_AOM_HIGHBITDEPTH
                         uint8_t *ref, int width, int height, int stride,
                         uint8_t *dst, int p_col, int p_row, int p_width,
                         int p_height, int p_stride, int subsampling_x,
                         int subsampling_y, int x_scale, int y_scale) {
  if (wm->wmtype <= AFFINE)
    if (!get_shear_params(wm)) return 1;
#if CONFIG_AOM_HIGHBITDEPTH
  if (use_hbd)
    return highbd_warp_erroradv(
        wm, ref, width, height, stride, dst, p_col, p_row, p_width, p_height,
        p_stride, subsampling_x, subsampling_y, x_scale, y_scale, bd);
#endif  // CONFIG_AOM_HIGHBITDEPTH
  return warp_erroradv(wm, ref, width, height, stride, dst, p_col, p_row,
                       p_width, p_height, p_stride, subsampling_x,
                       subsampling_y, x_scale, y_scale);
}

void av1_warp_plane(WarpedMotionParams *wm,
#if CONFIG_AOM_HIGHBITDEPTH
                    int use_hbd, int bd,
#endif  // CONFIG_AOM_HIGHBITDEPTH
                    uint8_t *ref, int width, int height, int stride,
                    uint8_t *pred, int p_col, int p_row, int p_width,
                    int p_height, int p_stride, int subsampling_x,
                    int subsampling_y, int x_scale, int y_scale, int ref_frm) {
#if CONFIG_AOM_HIGHBITDEPTH
  if (use_hbd)
    highbd_warp_plane(wm, ref, width, height, stride, pred, p_col, p_row,
                      p_width, p_height, p_stride, subsampling_x, subsampling_y,
                      x_scale, y_scale, bd, ref_frm);
  else
#endif  // CONFIG_AOM_HIGHBITDEPTH
    warp_plane(wm, ref, width, height, stride, pred, p_col, p_row, p_width,
               p_height, p_stride, subsampling_x, subsampling_y, x_scale,
               y_scale, ref_frm);
}

#if CONFIG_WARPED_MOTION
#define LEAST_SQUARES_ORDER 2

#define LS_MV_MAX 512  // max mv in 1/8-pel
#define LS_STEP 2

#define LS_SUM(a) ((a)*4 + LS_STEP * 2)
#define LS_SQUARE(a) \
  (((a) * (a)*4 + (a)*4 * LS_STEP + LS_STEP * LS_STEP * 2) >> 2)
#define LS_PRODUCT1(a, b) \
  (((a) * (b)*4 + ((a) + (b)) * 2 * LS_STEP + LS_STEP * LS_STEP) >> 2)
#define LS_PRODUCT2(a, b) \
  (((a) * (b)*4 + ((a) + (b)) * 2 * LS_STEP + LS_STEP * LS_STEP * 2) >> 2)

#if LEAST_SQUARES_ORDER == 2
static int find_affine_int(const int np, int *pts1, int *pts2, BLOCK_SIZE bsize,
                           int mvy, int mvx, WarpedMotionParams *wm, int mi_row,
                           int mi_col) {
  int32_t A[2][2] = { { 0, 0 }, { 0, 0 } };
  int32_t Bx[2] = { 0, 0 };
  int32_t By[2] = { 0, 0 };
  int i, n = 0;

  const int bw = block_size_wide[bsize];
  const int bh = block_size_high[bsize];
  const int suy = (mi_row * MI_SIZE + AOMMAX(bh, MI_SIZE) / 2 - 1) * 8;
  const int sux = (mi_col * MI_SIZE + AOMMAX(bw, MI_SIZE) / 2 - 1) * 8;
  const int duy = suy + mvy;
  const int dux = sux + mvx;

  // Assume the center pixel of the block has exactly the same motion vector
  // as transmitted for the block. First shift the origin of the source
  // points to the block center, and the origin of the destination points to
  // the block center added to the motion vector transmitted.
  // Let (xi, yi) denote the source points and (xi', yi') denote destination
  // points after origin shfifting, for i = 0, 1, 2, .... n-1.
  // Then if  P = [x0, y0,
  //               x1, y1
  //               x2, y1,
  //                ....
  //              ]
  //          q = [x0', x1', x2', ... ]'
  //          r = [y0', y1', y2', ... ]'
  // the least squares problems that need to be solved are:
  //          [h1, h2]' = inv(P'P)P'q and
  //          [h3, h4]' = inv(P'P)P'r
  // where the affine transformation is given by:
  //          x' = h1.x + h2.y
  //          y' = h3.x + h4.y
  //
  // The loop below computes: A = P'P, Bx = P'q, By = P'r
  // We need to just compute inv(A).Bx and inv(A).By for the solutions.
  int sx, sy, dx, dy;
  // Contribution from neighbor block
  for (i = 0; i < np && n < LEAST_SQUARES_SAMPLES_MAX; i++) {
    dx = pts2[i * 2] - dux;
    dy = pts2[i * 2 + 1] - duy;
    sx = pts1[i * 2] - sux;
    sy = pts1[i * 2 + 1] - suy;
    if (abs(sx - dx) < LS_MV_MAX && abs(sy - dy) < LS_MV_MAX) {
      A[0][0] += LS_SQUARE(sx);
      A[0][1] += LS_PRODUCT1(sx, sy);
      A[1][1] += LS_SQUARE(sy);
      Bx[0] += LS_PRODUCT2(sx, dx);
      Bx[1] += LS_PRODUCT1(sy, dx);
      By[0] += LS_PRODUCT1(sx, dy);
      By[1] += LS_PRODUCT2(sy, dy);
      n++;
    }
  }
  int64_t Px[2], Py[2];
  int64_t iDet, Det, v;
  int16_t shift;

  // These divided by the Det, are the least squares solutions
  Px[0] = (int64_t)A[1][1] * Bx[0] - (int64_t)A[0][1] * Bx[1];
  Px[1] = -(int64_t)A[0][1] * Bx[0] + (int64_t)A[0][0] * Bx[1];
  Py[0] = (int64_t)A[1][1] * By[0] - (int64_t)A[0][1] * By[1];
  Py[1] = -(int64_t)A[0][1] * By[0] + (int64_t)A[0][0] * By[1];

  // Compute Determinant of A
  Det = (int64_t)A[0][0] * A[1][1] - (int64_t)A[0][1] * A[0][1];
  if (Det == 0) return 1;
  iDet = resolve_divisor_64(labs(Det), &shift) * (Det < 0 ? -1 : 1);
  shift -= WARPEDMODEL_PREC_BITS;
  if (shift < 0) {
    iDet <<= (-shift);
    shift = 0;
  }

  v = Px[0] * iDet;
  wm->wmmat[2] = ROUND_POWER_OF_TWO_SIGNED_64(v, shift);
  v = Px[1] * iDet;
  wm->wmmat[3] = ROUND_POWER_OF_TWO_SIGNED_64(v, shift);
  v = (dux << WARPEDMODEL_PREC_BITS) - sux * wm->wmmat[2] - suy * wm->wmmat[3];
  wm->wmmat[0] = ROUND_POWER_OF_TWO_SIGNED(v, 3);

  v = Py[0] * iDet;
  wm->wmmat[4] = ROUND_POWER_OF_TWO_SIGNED_64(v, shift);
  v = Py[1] * iDet;
  wm->wmmat[5] = ROUND_POWER_OF_TWO_SIGNED_64(v, shift);
  v = (duy << WARPEDMODEL_PREC_BITS) - sux * wm->wmmat[4] - suy * wm->wmmat[5];
  wm->wmmat[1] = ROUND_POWER_OF_TWO_SIGNED(v, 3);

  wm->wmmat[6] = wm->wmmat[7] = 0;

  // Clamp values
  wm->wmmat[0] = clamp(wm->wmmat[0], -WARPEDMODEL_TRANS_CLAMP,
                       WARPEDMODEL_TRANS_CLAMP - 1);
  wm->wmmat[1] = clamp(wm->wmmat[1], -WARPEDMODEL_TRANS_CLAMP,
                       WARPEDMODEL_TRANS_CLAMP - 1);
  wm->wmmat[2] = clamp(wm->wmmat[2], -WARPEDMODEL_DIAGAFFINE_CLAMP,
                       WARPEDMODEL_DIAGAFFINE_CLAMP - 1);
  wm->wmmat[5] = clamp(wm->wmmat[5], -WARPEDMODEL_DIAGAFFINE_CLAMP,
                       WARPEDMODEL_DIAGAFFINE_CLAMP - 1);
  wm->wmmat[3] = clamp(wm->wmmat[3], -WARPEDMODEL_NONDIAGAFFINE_CLAMP,
                       WARPEDMODEL_NONDIAGAFFINE_CLAMP - 1);
  wm->wmmat[4] = clamp(wm->wmmat[4], -WARPEDMODEL_NONDIAGAFFINE_CLAMP,
                       WARPEDMODEL_NONDIAGAFFINE_CLAMP - 1);
  return 0;
}

#else

static int find_affine_int(const int np, int *pts1, int *pts2, BLOCK_SIZE bsize,
                           int mvy, int mvx, WarpedMotionParams *wm, int mi_row,
                           int mi_col) {
  int32_t A[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
  int32_t Bx[3] = { 0, 0, 0 };
  int32_t By[3] = { 0, 0, 0 };
  int i, n = 0, off;

  int64_t C00, C01, C02, C11, C12, C22;
  int64_t Px[3], Py[3];
  int64_t Det, v;
  const int bw = block_size_wide[bsize];
  const int bh = block_size_high[bsize];
  const int cy_offset = AOMMAX(bh, MI_SIZE) / 2 - 1;
  const int cx_offset = AOMMAX(bw, MI_SIZE) / 2 - 1;

  // Offsets to make the values in the arrays smaller
  const int ux = mi_col * MI_SIZE * 8, uy = mi_row * MI_SIZE * 8;
  // Let source points (xi, yi) map to destimation points (xi', yi'),
  //     for i = 0, 1, 2, .... n-1
  // Then if  P = [x0, y0, 1,
  //               x1, y1, 1
  //               x2, y2, 1,
  //                ....
  //              ]
  //          q = [x0', x1', x2', ... ]'
  //          r = [y0', y1', y2', ... ]'
  // the least squares problems that need to be solved are:
  //          [h1, h2, dx]' = inv(P'P)P'q and
  //          [h3, h4, dy]' = inv(P'P)P'r
  // where the affine transformation is given by:
  //          x' = h1.x + h2.y + dx
  //          y' = h3.x + h4.y + dy
  //
  // The loop below computes: A = P'P, Bx = P'q, By = P'r
  // We need to just compute inv(A).Bx and inv(A).By for the solutions.
  //
  int sx, sy, dx, dy;
  // Contribution from sample in current block
  sx = cx_offset * 8;
  sy = cy_offset * 8;
  dx = sx + mvx;
  dy = sy + mvy;
  if (abs(sx - dx) < LS_MV_MAX && abs(sy - dy) < LS_MV_MAX) {
    A[0][0] += LS_SQUARE(sx);
    A[0][1] += LS_PRODUCT1(sx, sy);
    A[0][2] += LS_SUM(sx);
    A[1][1] += LS_SQUARE(sy);
    A[1][2] += LS_SUM(sy);
    A[2][2] += 4;
    Bx[0] += LS_PRODUCT2(sx, dx);
    Bx[1] += LS_PRODUCT1(sy, dx);
    Bx[2] += LS_SUM(dx);
    By[0] += LS_PRODUCT1(sx, dy);
    By[1] += LS_PRODUCT2(sy, dy);
    By[2] += LS_SUM(dy);
    n++;
  }
  // Contribution from neighbor block
  for (i = 0; i < np && n < LEAST_SQUARES_SAMPLES_MAX; i++) {
    dx = pts2[i * 2] - ux;
    dy = pts2[i * 2 + 1] - uy;
    sx = pts1[i * 2] - ux;
    sy = pts1[i * 2 + 1] - uy;
    if (abs(sx - dx) < LS_MV_MAX && abs(sy - dy) < LS_MV_MAX) {
      A[0][0] += LS_SQUARE(sx);
      A[0][1] += LS_PRODUCT1(sx, sy);
      A[0][2] += LS_SUM(sx);
      A[1][1] += LS_SQUARE(sy);
      A[1][2] += LS_SUM(sy);
      A[2][2] += 4;
      Bx[0] += LS_PRODUCT2(sx, dx);
      Bx[1] += LS_PRODUCT1(sy, dx);
      Bx[2] += LS_SUM(dx);
      By[0] += LS_PRODUCT1(sx, dy);
      By[1] += LS_PRODUCT2(sy, dy);
      By[2] += LS_SUM(dy);
      n++;
    }
  }
  // Compute Cofactors of A
  C00 = (int64_t)A[1][1] * A[2][2] - (int64_t)A[1][2] * A[1][2];
  C01 = (int64_t)A[1][2] * A[0][2] - (int64_t)A[0][1] * A[2][2];
  C02 = (int64_t)A[0][1] * A[1][2] - (int64_t)A[0][2] * A[1][1];
  C11 = (int64_t)A[0][0] * A[2][2] - (int64_t)A[0][2] * A[0][2];
  C12 = (int64_t)A[0][1] * A[0][2] - (int64_t)A[0][0] * A[1][2];
  C22 = (int64_t)A[0][0] * A[1][1] - (int64_t)A[0][1] * A[0][1];

  // Scale by 1/16
  C00 = ROUND_POWER_OF_TWO_SIGNED(C00, 6);
  C01 = ROUND_POWER_OF_TWO_SIGNED(C01, 6);
  C02 = ROUND_POWER_OF_TWO_SIGNED(C02, 6);
  C11 = ROUND_POWER_OF_TWO_SIGNED(C11, 6);
  C12 = ROUND_POWER_OF_TWO_SIGNED(C12, 6);
  C22 = ROUND_POWER_OF_TWO_SIGNED(C22, 6);

  // Compute Determinant of A
  Det = C00 * A[0][0] + C01 * A[0][1] + C02 * A[0][2];
  if (Det == 0) return 1;

  // These divided by the Det, are the least squares solutions
  Px[0] = C00 * Bx[0] + C01 * Bx[1] + C02 * Bx[2];
  Px[1] = C01 * Bx[0] + C11 * Bx[1] + C12 * Bx[2];
  Px[2] = C02 * Bx[0] + C12 * Bx[1] + C22 * Bx[2];
  Py[0] = C00 * By[0] + C01 * By[1] + C02 * By[2];
  Py[1] = C01 * By[0] + C11 * By[1] + C12 * By[2];
  Py[2] = C02 * By[0] + C12 * By[1] + C22 * By[2];

  int16_t shift;
  int64_t iDet;
  iDet = resolve_divisor_64(labs(Det), &shift) * (Det < 0 ? -1 : 1);
  shift -= WARPEDMODEL_PREC_BITS;
  if (shift < 0) {
    iDet <<= (-shift);
    shift = 0;
  }

  v = Px[0] * iDet;
  wm->wmmat[2] = ROUND_POWER_OF_TWO_SIGNED_64(v, shift);
  v = Px[1] * iDet;
  wm->wmmat[3] = ROUND_POWER_OF_TWO_SIGNED_64(v, shift);
  v = Px[2] * iDet;
  wm->wmmat[0] = ROUND_POWER_OF_TWO_SIGNED_64(v, shift + 3);
  // Adjust x displacement for the offset
  off = (ux << WARPEDMODEL_PREC_BITS) - ux * wm->wmmat[2] - uy * wm->wmmat[3];
  wm->wmmat[0] += ROUND_POWER_OF_TWO_SIGNED(off, 3);

  v = Py[0] * iDet;
  wm->wmmat[4] = ROUND_POWER_OF_TWO_SIGNED_64(v, shift);
  v = Py[1] * iDet;
  wm->wmmat[5] = ROUND_POWER_OF_TWO_SIGNED_64(v, shift);
  v = Py[2] * iDet;
  wm->wmmat[1] = ROUND_POWER_OF_TWO_SIGNED_64(v, shift + 3);
  // Adjust y displacement for the offset
  off = (uy << WARPEDMODEL_PREC_BITS) - ux * wm->wmmat[4] - uy * wm->wmmat[5];
  wm->wmmat[1] += ROUND_POWER_OF_TWO_SIGNED(off, 3);
  wm->wmmat[6] = wm->wmmat[7] = 0;

  // Clamp values
  wm->wmmat[0] = clamp(wm->wmmat[0], -WARPEDMODEL_TRANS_CLAMP,
                       WARPEDMODEL_TRANS_CLAMP - 1);
  wm->wmmat[1] = clamp(wm->wmmat[1], -WARPEDMODEL_TRANS_CLAMP,
                       WARPEDMODEL_TRANS_CLAMP - 1);
  wm->wmmat[2] = clamp(wm->wmmat[2], -WARPEDMODEL_DIAGAFFINE_CLAMP,
                       WARPEDMODEL_DIAGAFFINE_CLAMP - 1);
  wm->wmmat[5] = clamp(wm->wmmat[5], -WARPEDMODEL_DIAGAFFINE_CLAMP,
                       WARPEDMODEL_DIAGAFFINE_CLAMP - 1);
  wm->wmmat[3] = clamp(wm->wmmat[3], -WARPEDMODEL_NONDIAGAFFINE_CLAMP,
                       WARPEDMODEL_NONDIAGAFFINE_CLAMP - 1);
  wm->wmmat[4] = clamp(wm->wmmat[4], -WARPEDMODEL_NONDIAGAFFINE_CLAMP,
                       WARPEDMODEL_NONDIAGAFFINE_CLAMP - 1);

  return 0;
}
#endif  // LEAST_SQUARES_ORDER == 2

int find_projection(const int np, int *pts1, int *pts2, BLOCK_SIZE bsize,
                    int mvy, int mvx, WarpedMotionParams *wm_params, int mi_row,
                    int mi_col) {
  int result = 1;
  switch (wm_params->wmtype) {
    case AFFINE:
      result = find_affine_int(np, pts1, pts2, bsize, mvy, mvx, wm_params,
                               mi_row, mi_col);
      break;
    default: assert(0 && "Invalid warped motion type!"); return 1;
  }
  if (result == 0) {
    if (wm_params->wmtype == ROTZOOM) {
      wm_params->wmmat[5] = wm_params->wmmat[2];
      wm_params->wmmat[4] = -wm_params->wmmat[3];
    }
    if (wm_params->wmtype == AFFINE || wm_params->wmtype == ROTZOOM) {
      // check compatibility with the fast warp filter
      if (!get_shear_params(wm_params)) return 1;
    }
  }

  return result;
}
#endif  // CONFIG_WARPED_MOTION
