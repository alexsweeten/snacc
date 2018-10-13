/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
   bwtmain.c
   main function of the bwte tool for computing the bwt in external memory

   Giovanni Manzini (giovanni.manzini@unipmn.it)
   27 Nov 2007
   04 Mar 2009
   27 Jan 2010
   

    Copyright (C) 2007,2008,2009,2010  Giovanni Manzini

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include "bwtext.h"
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>

// parameters to be obtained from the command line 
extern int _bwtext_Verbose;                // verbosity level
extern int _bwtext_Use_divsufsort;         // use divsvsufsort
extern double _bwtext_FilterPar1;          // parameter #1
extern int _bwtext_FilterPar1_flag;        // !=0 if parameter #1 was defined
extern double _bwtext_FilterPar2;          // parameter #2
extern int _bwtext_FilterPar2_flag;        // !=0 if parameter #2 was defined


int main(int argc, char *argv[])
{
  int64 tsize=0;
  time_t end,start;
  extern char *optarg;
  extern int optind, opterr, optopt;
  int32 c, text_filter, bwt_filter, num_opt, imem;
  filter *tfil, *bfil;
  char *fnam, *bwt_filename;

  /* ------------ set default values ------------- */
  _bwtext_Verbose=0;
  _bwtext_Use_divsufsort=1;
  bwt_filename = NULL;
  imem = 256;
  text_filter= -1; // default: autodetect 
  bwt_filter=0;    // default: no compression 
  _bwtext_FilterPar1_flag = _bwtext_FilterPar2_flag = 0;
  /* ------------- read options from command line ----------- */
  num_opt = opterr = 0;
  while ((c=getopt(argc, argv, "b:t:o:vlm:y:z:")) != -1) {
    switch (c) 
      {
      case 'm':
        imem = atoi(optarg); break;
      case 'v':
        _bwtext_Verbose++; break;
      case 'l':
        _bwtext_Use_divsufsort=0; break;
      case 't':
        text_filter = atoi(optarg); break;
      case 'b':
        bwt_filter = atoi(optarg); break;
      case 'o':
        bwt_filename = optarg; break;
      case 'y':
        _bwtext_FilterPar1 = atof(optarg);
	_bwtext_FilterPar1_flag =1; break;
      case 'z':
        _bwtext_FilterPar2 = atof(optarg);
	_bwtext_FilterPar2_flag =1; break;
      case '?':
        fprintf(stderr,"Unknown option: %c -main-\n", optopt);
        exit(1);
      }
    num_opt++;
  }
  if(optind<argc && argv[optind]!=NULL)
    fnam=argv[optind];
  else {
    fprintf(stderr, 
	    "Usage:\n\t%s [-v][-t fil][-b fil][-m imem]",  argv[0]);
    fprintf(stderr,"[-o bwtfile][-y val][-z val] file\n");
    fprintf(stderr,"\t-m imem     internal memory to be used "
                                  "in MB [def. %d]\n",imem);
    fprintf(stderr,"\t-t fil      text file filter [def. autodetect]\n");
    fprintf(stderr,"\t-b fil      bwt file filter [def. none]\n");
    fprintf(stderr,"\t-v          produces a verbose output\n");
    fprintf(stderr,"\t-o bwtfile  write bwt to bwtfile\n");
    fprintf(stderr,"\t-y val      bwt file filter parameter #1\n");
    fprintf(stderr,"\t-z val      bwt file filter parameter #2\n\n");
    _bwtext_print_available_filters(stderr);
    fprintf(stderr,"\n");
    return 0;
  }
  if(_bwtext_Verbose) {
    fprintf(stderr,"Command line: ");
    for(c=0;c<argc;c++)
      fprintf(stderr,"%s ",argv[c]);
    fprintf(stderr,"\n");
  }

  /* ---------- get filters -----------*/
  if(!(bfil=_bwtext_get_filter(bwt_filter,bwt_filename)))
    _bwtext_fatal_error("Invalid bwt filter (main)");
  if(!(tfil = _bwtext_get_filter(text_filter,fnam)))
    _bwtext_fatal_error("Invalid text filter (main)");
  
  /* ---------  start measuring time ------------- */
  start = time(NULL);                    /* this is calendar time */
  tsize = _bwtext_sufsort(fnam, tfil, bwt_filename, bfil,imem);
  end = time(NULL);
  if(_bwtext_Verbose>0) {
    fprintf(stderr,"Elapsed wallclock time: %ld seconds.\n", end-start);
    fprintf(stderr,"musec x symbol: %lf\n",(end-start)*(1000000.0/tsize));
  }

#if 0
  if(_bwtext_Verbose>1) {
    // print information on the IO Volume on Unix systems 
    printf("---- content of /proc/PID/io file ----\n");
    system("cat /proc/`pidof bwte`/io");
    printf("--------------------------------------\n");
  }
#endif
  return 0;
}





