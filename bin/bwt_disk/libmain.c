/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   the procedures here should be the main interface for using the 
   bwtext library, but some work still has to be done for unbwt
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "bwtext.h"


// parameters to be obtained from the command line 
extern int _bwtext_Verbose;                // verbosity level
extern int _bwtext_Use_divsufsort;         // use divsvsufsort
extern double _bwtext_FilterPar1;          // parameter #1
extern int _bwtext_FilterPar1_flag;        // !=0 if parameter #1 was defined
extern double _bwtext_FilterPar2;          // parameter #2
extern int _bwtext_FilterPar2_flag;        // !=0 if parameter #2 was defined
// variables initialized here
extern filter *_bwtext_Text_filter, *_bwtext_Bwt_filter; // filters 
extern char *_bwtext_Tname, *_bwtext_Bwtname;            // file names


// ===== main entry point to the bwt computation procedure ========
int64 _bwtext_sufsort(char *tnam, filter *tfil, char *bnam, filter *bfil, 
		      int32 imem)
{
  int64 _bwtext_sufsortmain(int32);


  // --------------- check imem
  if(imem<1) 
    _bwtext_fatal_error("Invalid amount of internal memory (_bwtext_sufsort)");

  // --------------- check filters
  if(bfil!=NULL)
    _bwtext_Bwt_filter=bfil;
  else
    _bwtext_fatal_error("Invalid bwt filter (_bwtext_sufsort)");
  if(tfil!=NULL) 
    _bwtext_Text_filter = tfil;
  else
    _bwtext_fatal_error("Invalid text filter (_bwtext_sufsort)");
  // report filters
  if(_bwtext_Verbose>1) {
    fprintf(stderr,"Using text filter: %s\n",_bwtext_Text_filter->name);
    fprintf(stderr,"Using bwt filter: %s\n",_bwtext_Bwt_filter->name);
  }

  // -------------- initialize text file name
  if(tnam!=NULL) 
    _bwtext_Tname = tnam;
  else
     _bwtext_fatal_error("Invalid text file name (_bwtext_sufsort)");   
    
  // -------------- initialize bwt file name
  if(bnam!=NULL) 
    _bwtext_Bwtname = bnam;
  else 
    asprintf(&_bwtext_Bwtname,"%s.bwt%s",tnam,_bwtext_Bwt_filter->extension);

  // -------------- do the sorting
  return _bwtext_sufsortmain(imem);
}
