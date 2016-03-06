/*
 *  Copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <assert.h>

#include "./bitwriter.h"
#if CONFIG_DAALA_EC
#include "entenc.h"
#include <string.h>
#endif

void vpx_start_encode(vpx_writer *br, uint8_t *source) {
  br->lowvalue = 0;
  br->range = 255;
  br->count = -24;
  br->buffer = source;
  br->pos = 0;
#if CONFIG_DAALA_EC
  od_ec_enc_init(&br->ec, 62025);
#else
  vpx_write_bit(br, 0);
#endif
}

void vpx_stop_encode(vpx_writer *br) {
#if CONFIG_DAALA_EC
  uint32_t daala_bytes;
  unsigned char *daala_data;
  daala_data = od_ec_enc_done(&br->ec, &daala_bytes);
  memcpy(br->buffer, daala_data, daala_bytes);
  br->pos = daala_bytes;
  od_ec_enc_clear(&br->ec);
  // Ensure there's no ambigous collision with any index marker bytes
  // this always wastes a byte which is not great
  br->buffer[br->pos++] = 0;
#else
  int i;

  for (i = 0; i < 32; i++) vpx_write_bit(br, 0);

  // Ensure there's no ambigous collision with any index marker bytes
  if ((br->buffer[br->pos - 1] & 0xe0) == 0xc0) br->buffer[br->pos++] = 0;
#endif
}
