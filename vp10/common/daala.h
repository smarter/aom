#if !defined(_daala_H)
# define _daala_H

# include "vpx/vpx_integer.h"
//# include "vp10/encoder/encint.h"

# define OD_DIVU(_x, _d) \
  (((_d) < OD_DIVU_DMAX)?(OD_DIVU_SMALL((_x),(_d))):((_x)/(_d)))

// from codec.h
/**The maximum number of color planes allowed in a single frame.*/
# define OD_NPLANES_MAX (4)

// fromn filter.h
typedef int32_t od_coeff;

// from daaladec.h
/**The decoder context.*/
typedef struct daala_dec_ctx daala_dec_ctx;

// from block_size.h
/*Possible block sizes, note that OD_BLOCK_NXN = log2(N) - 2.*/
#define OD_BLOCK_4X4 (0)
#define OD_BLOCK_8X8 (1)
#define OD_BLOCK_16X16 (2)
#define OD_BLOCK_32X32 (3)
#define OD_BLOCK_64X64 (4)
#define OD_BLOCK_SIZES (OD_BLOCK_64X64 + 1)

#endif
