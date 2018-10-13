
include ../../make_options

OBJECTS= main.o io.o buffer.o

# library for lzma decompression 
liblzmadec.a: ${OBJECTS} private.h lzmadec.h # ../../make_options
	ar rcs liblzmadec.a ${OBJECTS}

# pattern rule for all objects files
%.o: %.c *.h
	$(CC) -c $(CFLAGS) -Wno-unused-parameter $< -o $@ -DHAVE_ERRNO_H -DHAVE_LIMITS_H

clean: 
	rm -f *.o liblzmadec.a


            
