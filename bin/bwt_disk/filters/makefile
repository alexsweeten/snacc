
include ../make_options

OBJECTS= filters.o plain_filter.o gz_filter.o rc_filter.o \
         lzma_filter.o dnabwt_filter.o 

HEADERS = ../*.h ../make_options filters.h makefile

# IO filters
filters.a: ${OBJECTS} lzma/liblzmadec.a makefile 
	cp -f lzma/liblzmadec.a filters.a
	ar rs filters.a ${OBJECTS}


# lzmadec library
lzma/liblzmadec.a: lzma/*.c  lzma/*.h 
	make -C ./lzma/ 

# special tule for rc_filter that includes range.c
rc_filter.o: rc_filter.c range.c range.h ${HEADERS}
	 $(CC) -c $(CFLAGS) rc_filter.c -o rc_filter.o

# pattern rule for all objects files
%.o: %.c *.h ${HEADERS}
	$(CC) -c $(CFLAGS) $< -o $@


# delete everything that is creted here
clean: 
	rm -f *.o filters.a 
	make -C ./lzma/ clean

            
