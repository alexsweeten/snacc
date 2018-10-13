
include ../make_options

OBJS = isufsort.o qsufsort.o libdivsufsort-lite/divsufsort.o
HEADERS = ../*.h ../make_options makefile

# isufsort
isufsort.a: ${OBJS} makefile
	ar rcs isufsort.a ${OBJS}

# divsufsort library
libdivsufsort-lite/divsufsort.o: libdivsufsort-lite/*.c libdivsufsort-lite/*.h 
	make -C ./libdivsufsort-lite

# pattern rule for all objects files
%.o: %.c ${HEADERS}
	$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f *.o isufsort.a
	make -C ./libdivsufsort-lite clean
            





