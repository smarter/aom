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

// Simple Decoder
// ==============
//
// This is an example of a simple decoder loop. It takes an input file
// containing the compressed data (in IVF format), passes it through the
// decoder, and writes the decompressed frames to disk. Other decoder
// examples build upon this one.
//
// The details of the IVF format have been elided from this example for
// simplicity of presentation, as IVF files will not generally be used by
// your application. In general, an IVF file consists of a file header,
// followed by a variable number of frames. Each frame consists of a frame
// header followed by a variable length payload. The length of the payload
// is specified in the first four bytes of the frame header. The payload is
// the raw compressed data.
//
// Standard Includes
// -----------------
// For decoders, you only have to include `aom_decoder.h` and then any
// header files for the specific codecs you use. In this case, we're using
// aom.
//
// Initializing The Codec
// ----------------------
// The libaom decoder is initialized by the call to aom_codec_dec_init().
// Determining the codec interface to use is handled by AvxVideoReader and the
// functions prefixed with aom_video_reader_. Discussion of those functions is
// beyond the scope of this example, but the main gist is to open the input file
// and parse just enough of it to determine if it's a AVx file and which AVx
// codec is contained within the file.
// Note the NULL pointer passed to aom_codec_dec_init(). We do that in this
// example because we want the algorithm to determine the stream configuration
// (width/height) and allocate memory automatically.
//
// Decoding A Frame
// ----------------
// Once the frame has been read into memory, it is decoded using the
// `aom_codec_decode` function. The call takes a pointer to the data
// (`frame`) and the length of the data (`frame_size`). No application data
// is associated with the frame in this example, so the `user_priv`
// parameter is NULL. The `deadline` parameter is left at zero for this
// example. This parameter is generally only used when doing adaptive post
// processing.
//
// Codecs may produce a variable number of output frames for every call to
// `aom_codec_decode`. These frames are retrieved by the
// `aom_codec_get_frame` iterator function. The iterator variable `iter` is
// initialized to NULL each time `aom_codec_decode` is called.
// `aom_codec_get_frame` is called in a loop, returning a pointer to a
// decoded image or NULL to indicate the end of list.
//
// Processing The Decoded Data
// ---------------------------
// In this example, we simply write the encoded data to disk. It is
// important to honor the image's `stride` values.
//
// Cleanup
// -------
// The `aom_codec_destroy` call frees any memory allocated by the codec.
//
// Error Handling
// --------------
// This example does not special case any error return codes. If there was
// an error, a descriptive message is printed and the program exits. With
// few exceptions, aom_codec functions return an enumerated error status,
// with the value `0` indicating success.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#else
#define EMSCRIPTEN_KEEPALIVE
#endif

#include "aom/aom_decoder.h"
#include "aom/aomdx.h"

#include "../tools_common.h"
#include "../video_reader.h"
#include "./aom_config.h"
#include "av1/av1_dx_iface.c"
#include "../av1/common/onyxc_int.h"
#include "../av1/analyzer/analyzer.h"
#include "../video_common.h"


static const char *exec_name;

void usage_exit(void) {
  fprintf(stderr, "Usage: %s <infile> <outfile>\n", exec_name);
  exit(EXIT_FAILURE);
}

int frame_count = 0;
FILE *outfile = NULL;
aom_codec_ctx_t codec;
AvxVideoReader *reader = NULL;
const AvxInterface *decoder = NULL;
const AvxVideoInfo *info = NULL;
aom_image_t *img = NULL;
AV1_COMMON *cm = NULL;

int read_frame();

AnalyzerData analyzer_data;

void init_analyzer() {
  const int aligned_width = ALIGN_POWER_OF_TWO(info->frame_width, MI_SIZE_LOG2);
  const int aligned_height = ALIGN_POWER_OF_TWO(info->frame_height, MI_SIZE_LOG2);
  int mi_cols = aligned_width >> MI_SIZE_LOG2;
  int mi_rows = aligned_height >> MI_SIZE_LOG2;
  int mi_length = mi_cols * mi_rows;
  printf("init_analyzer: %d:%d (mi)\n", mi_cols, mi_rows);
  analyzer_data.mi_grid.buffer = aom_malloc(sizeof(AnalyzerMI) * mi_length);
  analyzer_data.mi_grid.length = mi_length;

  size_t size = aligned_width * aligned_height * 2;
  AnalyzerImagePlane *planes = analyzer_data.predicted_image.planes;

  planes[AOM_PLANE_Y].size = size;
  planes[AOM_PLANE_U].size = size;
  planes[AOM_PLANE_V].size = size;

  planes[AOM_PLANE_Y].stride = info->frame_width;
  planes[AOM_PLANE_U].stride = info->frame_width >> 1;
  planes[AOM_PLANE_V].stride = info->frame_width >> 1;

  planes[AOM_PLANE_Y].buffer = aom_malloc(size);
  planes[AOM_PLANE_U].buffer = aom_malloc(size);
  planes[AOM_PLANE_V].buffer = aom_malloc(size);
}

void dump_analyzer() {
  aom_codec_control(&codec, ANALYZER_SET_DATA, &analyzer_data);
  const int mi_rows = analyzer_data.mi_rows;
  const int mi_cols = analyzer_data.mi_cols;
  int r, c;
  for (r = 0; r < mi_rows; ++r) {
    for (c = 0; c < mi_cols; ++c) {
      // AnalyzerMI mi = analyzer_data.mi_grid.buffer[r * mi_cols + c];
      // printf("%d ", mi.mode);
      // printf("%3d:%-3d ", abs(mv.row), abs(mv.col));
      // printf("%d:%d ", mv.row, mv.col);
      // ...
    }
    // printf("\n");
  }
}

EMSCRIPTEN_KEEPALIVE
int open_file(char *file) {
  if (file == NULL) {
    file = "/tmp/input.ivf";
  }
  reader = aom_video_reader_open(file);
  if (!reader) die("Failed to open %s for reading.", file);

  info = aom_video_reader_get_info(reader);

  decoder = get_aom_decoder_by_fourcc(info->codec_fourcc);
  if (!decoder) die("Unknown input codec.");
  printf("Using %s\n", aom_codec_iface_name(decoder->codec_interface()));

  if (aom_codec_dec_init(&codec, decoder->codec_interface(), NULL, 0))
    die_codec(&codec, "Failed to initialize decoder.");

  init_analyzer();
  aom_codec_control(&codec, ANALYZER_SET_DATA, &analyzer_data);

  printf("Opened file %s okay.\n", file);
  return EXIT_SUCCESS;
}

void clear_mi_bits();

EMSCRIPTEN_KEEPALIVE
int read_frame() {
  if (!aom_video_reader_read_frame(reader)) {
    return EXIT_FAILURE;
  }
  img = NULL;
  aom_codec_iter_t iter = NULL;
  
  size_t frame_size = 0;
  const unsigned char *frame =
      aom_video_reader_get_frame(reader, &frame_size);
  clear_mi_bits();
  if (aom_codec_decode(&codec, frame, (unsigned int)frame_size, NULL, 0) != AOM_CODEC_OK) {
    die_codec(&codec, "Failed to decode frame.");
  }
  img = aom_codec_get_frame(&codec, &iter);
  if (img == NULL) {
    return EXIT_FAILURE;
  }
  ++frame_count;
  aom_codec_alg_priv_t* t = (aom_codec_alg_priv_t*)codec.priv;
  AVxWorker *const worker = &t->frame_workers[0];
  FrameWorkerData *const frame_worker_data = (FrameWorkerData *)worker->data1;
  cm = &frame_worker_data->pbi->common;

  return EXIT_SUCCESS;
}

EMSCRIPTEN_KEEPALIVE
int get_frame_count() {
  return frame_count;
}

EMSCRIPTEN_KEEPALIVE
int get_mi_cols_and_rows() {
  return analyzer_data.mi_cols << 16 | analyzer_data.mi_rows;
}

EMSCRIPTEN_KEEPALIVE
int get_tile_cols_and_rows_log2() {
  return analyzer_data.tile_cols_log2 << 16 | analyzer_data.tile_rows_log2;
}

typedef enum {
  GET_CODEC_BUILD_CONFIG,
} GetProperty;

EMSCRIPTEN_KEEPALIVE
const char *get_property(GetProperty v) {
  switch (v) {
    case GET_CODEC_BUILD_CONFIG:
      return aom_codec_build_config();
  }
}

typedef enum {
  GET_MI_MV,
  GET_MI_MODE,
  GET_MI_SKIP,
  GET_MI_REFERENCE_FRAME,
  GET_MI_BLOCK_SIZE,
  GET_MI_TRANSFORM_TYPE,
  GET_MI_TRANSFORM_SIZE,
  GET_MI_DERING_GAIN,
  GET_MI_BITS
} GetMIProperty;

EMSCRIPTEN_KEEPALIVE
int get_mi_property(GetMIProperty v, int c, int r, int i) {
  AnalyzerMI *mi =
    &analyzer_data.mi_grid.buffer[r * analyzer_data.mi_cols + c];
  switch (v) {
    case GET_MI_MV:
      return mi->mv[i].row << 16 | mi->mv[i].col;
    case GET_MI_MODE:
      return mi->mode;
    case GET_MI_SKIP:
      return mi->skip;
    case GET_MI_REFERENCE_FRAME:
      return mi->reference_frame[i];
    case GET_MI_BLOCK_SIZE:
      return mi->block_size;
    case GET_MI_TRANSFORM_TYPE:
      return mi->transform_type;
    case GET_MI_TRANSFORM_SIZE:
      return mi->transform_size;
    case GET_MI_DERING_GAIN:
      return mi->dering_gain;
    case GET_MI_BITS:
      return mi->bits;
  }
}

EMSCRIPTEN_KEEPALIVE
void clear_mi_bits() {
  const int mi_rows = analyzer_data.mi_rows;
  const int mi_cols = analyzer_data.mi_cols;
  int r, c;
  for (r = 0; r < mi_rows; ++r) {
    for (c = 0; c < mi_cols; ++c) {
      AnalyzerMI *mi =
        &analyzer_data.mi_grid.buffer[r * analyzer_data.mi_cols + c];
      mi->bits = 0;
    }
  }
}

EMSCRIPTEN_KEEPALIVE
unsigned char *get_plane(int plane) {
  return img->planes[plane];
}

EMSCRIPTEN_KEEPALIVE
int get_plane_stride(int plane) {
  return img->stride[plane];
}

EMSCRIPTEN_KEEPALIVE
int get_plane_width(int plane) {
  return aom_img_plane_width(img, plane);
}

EMSCRIPTEN_KEEPALIVE
int get_plane_height(int plane) {
  return aom_img_plane_height(img, plane);
}

EMSCRIPTEN_KEEPALIVE
int get_frame_width() {
  return info->frame_width;
}

EMSCRIPTEN_KEEPALIVE
int get_frame_height() {
  return info->frame_height;
}

EMSCRIPTEN_KEEPALIVE
unsigned char *get_predicted_plane_buffer(int plane) {
  return analyzer_data.predicted_image.planes[plane].buffer;
}

EMSCRIPTEN_KEEPALIVE
int get_predicted_plane_stride(int plane) {
  return analyzer_data.predicted_image.planes[plane].stride;
}

EMSCRIPTEN_KEEPALIVE
int main(int argc, char **argv) {
  // ...
  if (argc == 2) {
    open_file(argv[1]);
    while (!read_frame()) {
      printf("%d\n", frame_count);
      dump_analyzer();
    }
  }
}

EMSCRIPTEN_KEEPALIVE
void quit() {
  printf("Processed %d frames.\n", frame_count);
  if (aom_codec_destroy(&codec)) 
    die_codec(&codec, "Failed to destroy codec");

  // printf("Play: ffplay -f rawvideo -pix_fmt yuv420p -s %dx%d %s\n",
  //        info->frame_width, info->frame_height, argv[2]);

  aom_video_reader_close(reader);

  fclose(outfile);
}