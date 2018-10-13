# Makefile for suftest

# sources
OBJS					= divsufsort.o suftest.o
TARGET				= suftest
MAKEFILE			= Makefile

# options
CC						= gcc
#CXX						= g++
#OUTPUT_OPTION	= -o $@
CFLAGS				= -O3 -fomit-frame-pointer
#CXXFLAGS			= -O3 -fomit-frame-pointer
CPPFLAGS			= -Wall -DNDEBUG
LDFLAGS				= 
LDLIBS				= 
#TARGET_ARCH		=


# targets
.PHONY: all
all: $(TARGET)
$(TARGET): $(OBJS)

distclean: clean
clean:
	$(RM) $(TARGET) $(OBJS)

# dependencies
$(OBJS): divsufsort.h $(MAKEFILE)
