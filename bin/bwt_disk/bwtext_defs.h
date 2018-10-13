/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   Type defintions and prototypes of externally visible
   functions of the bwtext library.   
   Giovanni Manzini (giovanni.manzini@unipmn.it)
   4 Mar 2009
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <stdint.h>
#include <stdio.h>

/* ===============================================================
   struct representing a filter for reading/writing


   Notes:
   name is used to display the available filters at the prompt
   extension is used to create filenames
   open needs only to supports two modes: 'r' (read only) and 'w' (write only)
   read/write return the number of actually read/written uint8_t's and return
   -1 on errors.
   rewind needs to be supported only for files opened in read-only modes
   ================================================================ */
typedef struct {
  char *name;                                   // name of the filter
  char *extension;                              // associated extension
  void *(*open)(char *fname, int32_t mode);     // open returns NULL on failure
  int32_t (*read)(void *, uint8_t *b, int32_t n);   // read n bytes
  int32_t (*write)(void *, uint8_t *b, int32_t n);  // write n bytes
  void (*rewind)(void *);                      // rewind a read-only file
  int64_t (*close)(void *);      // close file. return >=0 on success
                                 // if available return current file position
} filter;


/* ================================================================
   struct representing a text file
   ================================================================ */
typedef struct {
  void *filep;                // pointer to the text file
  char *filename;             // associated filename
  filter *fil;                // filter used for reading/writing
  uint8_t *textbuffer;        // buffer used for I/O
  int32_t buf_size;           // buffer size
  int32_t buf_cur;            // pointer in the buffer
  uint64_t text_pos;          // current position in the text
} textfile;


/* ================================================================
   struct representing a bwt file
   ================================================================ */
typedef struct {
  void *filep;              // pointer to the bwt file on disk  
  char *filename;           // associated filename
  filter *fil;              // filter used for reading/writing
  uint64_t size;            // uncompressed size (text lenght+1)
  uint64_t eofpos;
} bwtfile;



// ===========================================================
// prototypes of functions provided by the bwtext library
// the full documentation of these functions is available 
// in the doc subdirectory 
// ===========================================================


// construction of the bwt using at most imem bytes of internal memory
int64_t 
_bwtext_sufsort(char *tnam, filter *tf, char *bnam, filter *bf, int32_t imem);

// functions related to filters 
filter *_bwtext_get_filter(int32_t opt, char *filename);
void _bwtext_print_available_filters(FILE *f);


//  read access to bwtfiles
bwtfile *_bwtext_bwt_open_read(char *fname,filter *f);
int32_t _bwtext_bwt_read(bwtfile *,uint8_t *, int32_t);
void _bwtext_bwt_close(bwtfile *b); 




