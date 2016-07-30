#ifndef ANALYZER_H_
#define ANALYZER_H_

struct AV1Decoder;

typedef enum {
  ANALYZER_OK,
  ANALYZER_ERROR
} AnalyzerError;

typedef uint8_t AnalyzerPredictionMode;
typedef int8_t AnalyzerReferenceFrame;
typedef int8_t AnalyzerInterpFilter;

typedef uint8_t AnalyzerTransformType;
typedef uint8_t AnalyzerTransformSize;

typedef struct AnalyzerImagePlane {
  unsigned char *buffer;
  uint32_t stride;
  size_t size;
} AnalyzerImagePlane;

typedef struct AnalyzerImage {
  AnalyzerImagePlane planes[4];
} AnalyzerImage;

/**
 * Motion Vector (MV)
 */
typedef struct AnalyzerMV {
  int16_t row;
  int16_t col;
} AnalyzerMV;

/**
 * Mode Info (MI)
 */
typedef struct AnalyzerMI {
  AnalyzerMV mv[2];
  AnalyzerReferenceFrame reference_frame[2];
  AnalyzerInterpFilter filter;
  AnalyzerPredictionMode mode;
  int8_t dering_gain;
  int8_t skip;
  uint8_t block_size;

  AnalyzerTransformType transform_type;
  AnalyzerTransformSize transform_size;
  size_t bits;
} AnalyzerMI;

typedef struct AnalyzerMIBuffer {
  AnalyzerMI *buffer;
  uint32_t length;
} AnalyzerMIBuffer;

/**
 * Holds everything that is needed by the stream analyzer.
 */
typedef struct AnalyzerData {
  AnalyzerImage image;
  AnalyzerImage predicted_image;
  AnalyzerMIBuffer mi_grid;
  uint32_t mi_rows;
  uint32_t mi_cols;
  uint32_t tile_rows_log2;
  uint32_t tile_cols_log2;
} AnalyzerData;

AnalyzerError analyzer_record_predicted_block(struct AV1Decoder *pbi,
                                              uint32_t plane,
                                              uint32_t subsampling_x,
                                              uint32_t subsampling_y,
                                              uint32_t mi_col,
                                              uint32_t mi_row,
                                              uint32_t col,
                                              uint32_t row,
                                              unsigned char* dst,
                                              uint32_t dst_stride,
                                              uint32_t transform_size);

AnalyzerError analyzer_record_frame(struct AV1Decoder *pbi);
AnalyzerError analyzer_record_mi_bits(struct AV1Decoder *pbi,
                                      uint32_t plane,
                                      uint32_t mi_col,
                                      uint32_t mi_row,
                                      uint32_t bits);
#endif  // ANALYZER_H_