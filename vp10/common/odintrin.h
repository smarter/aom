#ifndef VP10_COMMON_ODINTRIN_H_
#define VP10_COMMON_ODINTRIN_H_

#include <math.h>
#include "vp10/common/enums.h"
#include "vpx/vpx_integer.h"
#include "vpx_dsp/vpx_dsp_common.h"
#include "vpx_ports/bitops.h"
#include "vpx_mem/vpx_mem.h"

# if !defined(M_LOG2E)
#  define M_LOG2E (1.4426950408889634073599246810019)
# endif

/*Smallest blocks are 4x4*/
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

# define OD_DIVU_DMAX (1024)

typedef int od_coeff;

typedef int16_t od_dering_in;

#define OD_DIVU_SMALL(_x, _d) ((_x) / (_d))

# define OD_DIVU(_x, _d) \
  (((_d) < OD_DIVU_DMAX)?(OD_DIVU_SMALL((_x),(_d))):((_x)/(_d)))

#define OD_MINI VPXMIN
#define OD_MAXI VPXMAX
#define OD_CLAMPI(min, val, max) clamp((val), (min), (max))

# define OD_CLZ0 (1)
# define OD_CLZ(x) (-get_msb(x))
# define OD_ILOG_NZ(x) (OD_CLZ0 - OD_CLZ(x))
/*Note that __builtin_clz is not defined when x == 0, according to the gcc
   documentation (and that of the x86 BSR instruction that implements it), so
   we have to special-case it.
  We define a special version of the macro to use when x can be zero.*/
# define OD_ILOG(x) ((x) ? OD_ILOG_NZ(x) : 0)

# define OD_LOG2(x) (M_LOG2E*log(x))

/*Enable special features for gcc and compatible compilers.*/
# if defined(__GNUC__) && defined(__GNUC_MINOR__) && defined(__GNUC_PATCHLEVEL__)
#  define OD_GNUC_PREREQ(maj, min, pat)                             \
  ((__GNUC__ << 16) + (__GNUC_MINOR__ << 8) + __GNUC_PATCHLEVEL__ >= ((maj) << 16) + ((min) << 8) + pat)
# else
#  define OD_GNUC_PREREQ(maj, min, pat) (0)
# endif

#if OD_GNUC_PREREQ(3, 4, 0)
# define OD_WARN_UNUSED_RESULT __attribute__((__warn_unused_result__))
#else
# define OD_WARN_UNUSED_RESULT
#endif

#if OD_GNUC_PREREQ(3, 4, 0)
# define OD_ARG_NONNULL(x) __attribute__((__nonnull__(x)))
#else
# define OD_ARG_NONNULL(x)
#endif

# if defined(OD_ENABLE_ASSERTIONS)
#  if OD_GNUC_PREREQ(2, 5, 0)
__attribute__((noreturn))
#  endif
void od_fatal_impl(const char *_str, const char *_file, int _line);

#  define OD_FATAL(_str) (od_fatal_impl(_str, __FILE__, __LINE__))

#  define OD_ASSERT(_cond) \
  do { \
    if (!(_cond)) { \
      OD_FATAL("assertion failed: " # _cond); \
    } \
  } \
  while (0)

#  define OD_ASSERT2(_cond, _message) \
  do { \
    if (!(_cond)) { \
      OD_FATAL("assertion failed: " # _cond "\n" _message); \
    } \
  } \
  while (0)

#  define OD_ALWAYS_TRUE(_cond) OD_ASSERT(_cond)

# else
#  define OD_ASSERT(_cond)
#  define OD_ASSERT2(_cond, _message)
#  define OD_ALWAYS_TRUE(_cond) ((void)(_cond))
# endif

/** Copy n elements of memory from src to dst. The 0* term provides
    compile-time type checking  */
#if !defined(OVERRIDE_OD_COPY)
# define OD_COPY(dst, src, n) \
  (memcpy((dst), (src), sizeof(*(dst))*(n) + 0*((dst) - (src))))
#endif

/** Copy n elements of memory from src to dst, allowing overlapping regions.
    The 0* term provides compile-time type checking */
#if !defined(OVERRIDE_OD_MOVE)
# define OD_MOVE(dst, src, n) \
 (memmove((dst), (src), sizeof(*(dst))*(n) + 0*((dst) - (src)) ))
#endif

/** Linkage will break without this if using a C++ compiler, and will issue
 * warnings without this for a C compiler*/
#if defined(__cplusplus)
# define OD_EXTERN extern
#else
# define OD_EXTERN
#endif

/*16x16 multiplication where the result fits in 16 bits, without rounding.*/
# define OD_MULT16_16_Q15(a,b) \
  (((int16_t)(a)*((int32_t)(int16_t)(b))) >> 15)
/*16x16 multiplication where the result fits in 16 bits, without rounding.*/
# define OD_MULT16_16_Q16(a,b) \
  ((((int16_t)(a))*((int32_t)(int16_t)(b))) >> 16)
/*Shift x right by shift (with rounding)*/
# define OD_SHR_ROUND(x, shift) \
  ((int32_t)(((x) + (1 << (shift) >> 1)) >> (shift)))
/*Shift x right by shift (without rounding) or left by -shift if shift
  is negative.*/
# define OD_VSHR(x, shift) \
  ((shift) > 0 ? (int32_t)((x) >> (shift)) \
  : (int32_t)((x) << -(shift)))
/*Shift x right by shift (with rounding) or left by -shift if shift
  is negative.*/
# define OD_VSHR_ROUND(x, shift) \
  ((shift) > 0 ? (int32_t)(((x) + (1 << (shift) >> 1)) >> (shift)) \
  : (int32_t)((x) << -(shift)))

/*All of these macros should expect floats as arguments.*/
/*These two should compile as a single SSE instruction.*/
# define OD_MINF(a, b) ((a) < (b) ? (a) : (b))
# define OD_MAXF(a, b) ((a) > (b) ? (a) : (b))

# define OD_DIV_R0(x, y) (((x) + OD_FLIPSIGNI((((y) + 1) >> 1) - 1, (x)))/(y))

# define OD_SIGNMASK(a) (-((a) < 0))
# define OD_FLIPSIGNI(a, b) (((a) + OD_SIGNMASK(b)) ^ OD_SIGNMASK(b))

# define OD_MULT16_16_Q15(a,b) \
  (((int16_t)(a)*((int32_t)(int16_t)(b))) >> 15)

/* Multiplies 16-bit a by 32-bit b and keeps bits [16:47]. */
# define OD_MULT16_32_Q16(a, b) ((int16_t)(a)*(int64_t)(int32_t)(b) >> 16)

/** Copy n elements of memory from src to dst. The 0* term provides
    compile-time type checking  */
#if !defined(OVERRIDE_OD_COPY)
# define OD_COPY(dst, src, n) \
  (memcpy((dst), (src), sizeof(*(dst))*(n) + 0*((dst) - (src))))
#endif

/** Set n elements of dst to zero */
#if !defined(OVERRIDE_OD_CLEAR)
# define OD_CLEAR(dst, n) (memset((dst), 0, sizeof(*(dst))*(n)))
#endif

/** Silence unused parameter/variable warnings */
# define OD_UNUSED(expr) (void)(expr)

// from codec.h
/**The maximum number of color planes allowed in a single frame.*/
# define OD_NPLANES_MAX (4)

// fromn filter.h
typedef int32_t od_coeff;

// from block_size.h
/*Possible block sizes, note that OD_BLOCK_NXN = log2(N) - 2.*/
#define OD_BLOCK_4X4 (0)
#define OD_BLOCK_8X8 (1)
#define OD_BLOCK_16X16 (2)
#define OD_BLOCK_32X32 (3)
#define OD_BLOCK_64X64 (4)
#define OD_BLOCK_SIZES (OD_BLOCK_64X64 + 1)

#endif