/* Version 30 October 2001 */
// modified by G. Manzini for use winthin the bwtext librery on Sept. 29 2008


#ifndef RANGE_HEADER
#define RANGE_HEADER

#include <limits.h>

/* Static or adaptive model? */
#define  STATIC    0U
#define  ADAPT     1U

/* Error codes */
#define RC_OK        0U
#define RC_ERROR     1U
#define RC_IO_ERROR  2U

typedef int32 Sint;
typedef uint32 Uint;
typedef uint8 uc; 

// ======  this is how Sint Uint and uc were defined
/* #if UINT_MAX >= 0xffffffffU */
/* 	typedef   signed int  Sint; */
/* 	typedef unsigned int  Uint; */
/* #else */
/* 	typedef   signed long Sint; */
/* 	typedef unsigned long Uint; */
/* #endif */
/* typedef unsigned char uc; */




#define  TOP       ((Uint)1 << 24)
#define  BOT       ((Uint)1 << 16)

// in rc_encoder and rc_decoder I have redefined 
// passed as a uint64 since it is possible that we encode more than
// 2^32 symbols. However, that field is updated at each symbol 
// but it is never used 
typedef struct {
  Uint  low, range, error;
  FILE *fp;
  uc   *ptr, *eptr;
  uint64 passed;
} rc_encoder;

typedef struct {
  Uint  low, code, range, error;
  FILE *fp;
  uc   *ptr, *eptr;
  uint64 passed;
} rc_decoder;

typedef struct {
  Uint *freq, totFreq, incr, maxFreq, adapt, halfFreq, firstFreq, lastFreq;
  Uint lastSym, lastCumFreq;
  Uint nsym, nsym2, nsym3, nsym4, nsym23;
} rc_model;

static Uint ModelInit    (rc_model *, Uint, Uint *, Uint *, Uint, Uint, uc);
static void StartEncode  (rc_encoder *, FILE *, uc *, Uint);
static void FinishEncode (rc_encoder *);
static void StartDecode  (rc_decoder *, FILE *, uc *, Uint);
static Sint EncodeSymbol (rc_encoder *, rc_model *, Uint);
static Sint DecodeSymbol (rc_decoder *, rc_model *);

#endif
