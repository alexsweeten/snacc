# --------------------------------------------------------------------------
# makefile for bwtext and related tools
# Tools that only read/write text/bwt data should be linked to 
# bwtext.a. The objects in SORT_OBJS are required only for bwt construction
# --------------------------------------------------------------------------   

HEADERS = *.h makefile make_options 
OBJS = bwt_class.o text_class.o kmp_class.o global.o  
SORT_OBJS = bwtmerge.o libmain.o isufsort/isufsort.a
PROGS =  bwte unbwti text_conv text_count text_rev

include make_options


all : ${PROGS} 


# --- progs

bwte: bwtmain.o ${SORT_OBJS} bwtext.a makefile
	$(CC) $(CFLAGS) -o bwte bwtmain.o ${SORT_OBJS} bwtext.a  -lz 

unbwti: tools/unbwt.o bwtext.a makefile
	$(CC) $(CFLAGS) -o unbwti tools/unbwt.o  bwtext.a -lz

# pattern rule for the text_* tools
text_%: tools/text_%.o bwtext.a
	$(CC) $(CFLAGS) -o $@  $< bwtext.a -lz


# --- archives 

# library containing everything except sorting
bwtext.a: ${OBJS} filters/filters.a
	cp -f filters/filters.a bwtext.a
	ar rs bwtext.a ${OBJS} 

# library containing filters definitions
filters/filters.a: filters/*.c filters/filters.o ${HEADERS} \
                   filters/lzma/*.c  filters/lzma/*.h
	make -C ./filters/

# library containing internal suffix sort algorithms
isufsort/isufsort.a: isufsort/*.c ${HEADERS}  \
                     isufsort/libdivsufsort-lite/*.c \
                     isufsort/libdivsufsort-lite/*.h	
	make -C ./isufsort/


# pattern rule for all objects files
%.o: %.c ${HEADERS}
	$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f *.[oa] ${PROGS}
	rm -f ./tools/*.o
	make -C ./filters clean
	make -C ./isufsort clean

tarfile:
	tar zcvf bwtdisk.0.9.0.tgz ?akefile *.c *.h make_options \
        tools/*.c \
        filters/*.[ch] filters/makefile filters/lzma/* \
        --exclude=filters/lzma/*.[oa] \
	isufsort/*.c isufsort/makefile isufsort/libdivsufsort-lite/* \
        --exclude=isufsort/libdivsufsort-lite/*.[oa]  \
        README COPYING CHANGES doc/bwtdisk.pdf
            



