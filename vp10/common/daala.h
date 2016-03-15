 #if !defined(_daala_H)
 # define _daala_H

#include "vpx/vpx_integer.h"

// from internal.h
/*Smallest blocks are 4*/
# define OD_LOG_BSIZE0 (2)
/*There are 5 block sizes total (4x4, 8x8, 16x16, 32x32 and 64x64).*/
# define OD_NBSIZES    (5)
/*The log of the maximum length of the side of a block.*/
# define OD_LOG_BSIZE_MAX (OD_LOG_BSIZE0 + OD_NBSIZES - 1)
/*The maximum length of the side of a block.*/
# define OD_BSIZE_MAX     (1 << OD_LOG_BSIZE_MAX)

# define OD_COEFF_SHIFT (4)
# define OD_COEFF_SCALE (1 << OD_COEFF_SHIFT)

# define OD_LIMIT_BSIZE_MIN (OD_BLOCK_4X4)
# define OD_LIMIT_BSIZE_MAX (OD_BLOCK_64X64)

# define OD_DISABLE_CFL (1)

// from codec.h
/**The maximum number of color planes allowed in a single frame.*/
# define OD_NPLANES_MAX (4)

// fromn filter.h
typedef int32_t od_coeff;

// from daalaenc.h
/**The encoder context.*/
typedef struct daala_enc_ctx daala_enc_ctx;

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
