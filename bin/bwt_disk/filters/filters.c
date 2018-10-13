#include "filters.h"
#include <string.h>

// +++ update  if you want to use more filters 
#define NUM_TEXT_FILTERS 5


// prototypes of filter registration functions 
void _bwtext_register_plain_filter(filter *);
void _bwtext_register_gz_filter(filter *);
void _bwtext_register_rcrle_filter(filter *t);
void _bwtext_register_lzma_filter(filter *t);
void _bwtext_register_dnabwt_filter(filter *t);
// +++ add here the prototypes for your filter registration functions



// These filters work but they are very slow and require popen/pclose
// source code in in piped[gz|lzma].c
void _bwtext_register_pipedgz_filter(filter *t);
void _bwtext_register_pipedlzma_filter(filter *t);



/* *****************************************************************
   create and return an array containing pointers to the 
   available text filters
   ***************************************************************** */ 
static filter *_bwtext_register_text_filters()
{
  void _bwtext_fatal_error(char *s);
  void _bwtext_out_of_mem(char *f);
  filter *tf;
  int i=0;

  tf = (filter *) malloc(NUM_TEXT_FILTERS*sizeof(*tf));
  if(tf==NULL)
    _bwtext_out_of_mem("_bwtext_register_text_filters");

  // init text filters (they need to support read/rewind operations)
  _bwtext_register_plain_filter(&tf[i++]); // this must be the first filter
  _bwtext_register_gz_filter(&tf[i++]);
  _bwtext_register_rcrle_filter(&tf[i++]);
  _bwtext_register_dnabwt_filter(&tf[i++]);
  _bwtext_register_lzma_filter(&tf[i++]);
    
  // +++ add here additional calls to registration functions
  // +++ don't forget to update NUM_TEXT_FILTERS above


  if(i!=NUM_TEXT_FILTERS) 
    _bwtext_fatal_error("_bwtext_register_text_filters");
  return tf;
}


// variable storing the arrays of available filters
static filter *Text_filters;


// output list of available text filters
void _bwtext_print_available_filters(FILE *f)
{
  int i;

  // if necessary init array of available text filters
  if(Text_filters==NULL)
    Text_filters = _bwtext_register_text_filters();
  assert(Text_filters!=NULL);

  fprintf(f,"Available filters:\n");
  for(i=0;i<NUM_TEXT_FILTERS;i++)  {
    fprintf(f,"\t%2d.  %-36s",i,Text_filters[i].name);
    fprintf(f,"  [extension: %-5s]\n",strlen(Text_filters[i].extension)
	    ?Text_filters[i].extension:"None");
  }
}


// guess text filter according to the filename suffix 
// if none matches use filter # 0 (that should be the plain filter)
static filter *guess_filter(char *fname)
{
  char *suffix;
  int i; size_t fname_length=strlen(fname);

  // if necessary init array of available text filters
  if(Text_filters==NULL)
    Text_filters = _bwtext_register_text_filters();
  assert(Text_filters!=NULL);

  // check filters looking at extensions
  for(i=1;i<NUM_TEXT_FILTERS;i++) {
    // create suffix
    asprintf(&suffix, "%s",Text_filters[i].extension);
    if(fname_length>strlen(suffix)) {
      char *fsuffix = fname + (fname_length-strlen(suffix));
      assert(strlen(fsuffix)==strlen(suffix));
      if(strcmp(fsuffix,suffix)==0) { // compare with file suffix
	free(suffix);
	return &Text_filters[i];       // found
      }
    }
    free(suffix); // free suffix
  }
  return &Text_filters[0]; // this must be the plain filter 
}

// get filter using opt (if >=0) or filename 
filter *_bwtext_get_filter(int32 opt, char *fnam)
{

  // if necessary init array of available text filters
  if(Text_filters==NULL)
    Text_filters = _bwtext_register_text_filters();
  assert(Text_filters!=NULL);

  // if opt<0 guess using file extension 
  if(opt<0) {
    if(fnam==NULL) return NULL;
    else return guess_filter(fnam);
  }
  // otherwise select text filter from Text_filter array
  if(opt>=NUM_TEXT_FILTERS) 
    return NULL;
  else 
    return &Text_filters[opt];
}

