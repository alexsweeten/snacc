/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   Constants/prototypes used internally by the bwtext library 
   Giovanni Manzini  (giovanni.manzini@unipmn.it)
   27 Nov 2007
   04 Mar 2009
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <assert.h>
#include <stdlib.h>
#include "bwtext_defs.h"    // global type definitions

// alphabet size (cannot be changed!)
#define _BW_ALPHA_SIZE 256 

// buffers size (to do: experiment with these size )
#define _bwtext_TEXT_BUFFER_SIZE (1024*1024)
#define _bwtext_BWT_BUFFER_SIZE (64*1024)


// prefix size and buffer size for kmp_class
#define _bwtext_PFX_BUFFER_SIZE 1024         // length of stored prefix 
#define _bwtext_KMP_BIT_BUF_SIZE (128*1024)  // size in bits of the buffer 


/* ---------- ridefinition of types ----------- */
typedef uint64_t         uint64;
typedef int64_t          int64;
typedef int32_t          int32;
typedef uint32_t         uint32;
typedef uint16_t         uint16;
typedef int8_t           int8;
typedef uint8_t          uint8;
typedef unsigned int     boolean;


/* ********************************************************************
   data structure used to determine which substrings of the input text 
   fully match  a given string. The structure is used in a writing phase 
   in which for each matching substring we record a bit value, and in a 
   revisiting phase in which for each matching substring we retrieve 
   the previously stored bit
   ******************************************************************** */ 
typedef struct {
  uint8 *string;
  uint32 size;            // length of string
  uint16 *kmp_shift;      // kmp shift vector
  uint32 current;         // index in kmp shift vector
  uint32 *bit_buffer;     // bit-buffer split in 32 bits blocks 
  uint32 bit_buffer_size; // size of bit-buffer in bit
  uint32 cur_bit;         // current position in buffer
  uint32 stored_bits;     // # bits stored at the end of the writing phase
  int64 pending0;
  int64 pending1;
  uint64 chars_seen;
} _bwtext_kbuf;



// constant and macro for marking groups
//#define SETMASK32 (1u << 31)
//#define CLEARMASK32 (~(SETMASK))

// prototypes for error messages
void _bwtext_fatal_error(char *s);
void _bwtext_out_of_mem(char *f);
void _bwtext_ioerror(char *file, char *fun);

// prototypes from text_class
void _bwtext_text_open();
void _bwtext_text_rewind();
int32 _bwtext_text_getc();
void _bwtext_text_close();
uint32 _bwtext_text_copy_reverse(uint8 *t, uint32 size);
boolean _bwtext_text_checkpos(uint64 p);
boolean _bwtext_text_moretocome(void);
textfile *_bwtext_text_open_write(char *fname, filter *f, int32_t bufsize);
void _bwtext_text_putc(textfile *t, int32_t c); // return c or EOF on error
void _bwtext_text_close_write(textfile *t); 



// prototypes from bwt_class
int64 _bwtext_bwt_init_temp(uint32 size, uint32 eof, uint8 *);
void _bwtext_bwt_finalize_temp(); 
void _bwtext_bwt_write_auxfile(uint64 *); 
int64 _bwtext_bwtmerge_temp(uint32 *gaps, uint32 *cur_bwt, uint32 n, 
			     uint32 cur_rank0, uint32 cur_last);


// exit codes
#define _bwtext_IOERROR 1
#define _bwtext_NOMEMERROR 2

#ifndef min
#define min(a, b) ((a)<=(b) ? (a) : (b))
#endif

#ifndef max
#define max(a, b) ((a)>=(b) ? (a) : (b))
#endif




