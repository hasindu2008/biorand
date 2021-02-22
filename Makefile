CC		 = gcc
CXX       = g++
CFLAGS   = -g -rdynamic -Wall -O2 -std=c++11
CPPFLAGS =

#include config.mk

#LIB_LDFLAGS = -lghts

# INC = $(HDF5_INC)
#$(HTS_INC)
LDFLAGS += $(LIBS) -lpthread  -lz

#SRC = $(wildcard *.c)
SRC = main.c biorand.c filterpaf.c filterfq.c comparesam.c olp.c idat.c
OBJ = $(SRC:.c=.o)
BINARY = biorand
DEPS = misc.h biorand.h kseq.h

.PHONY: clean distclean format test

all : $(BINARY) bin/clean_fscache bin/enable_perf

$(BINARY) : $(OBJ)
	$(CXX) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@


%.o : %.c $(DEPS)
	$(CXX) $(CFLAGS) $(CPPFLAGS) $< -c


bin/clean_fscache : clean_fscache.c
	$(CC) -Wall $< $(LDFLAGS) -o $@

bin/enable_perf : enable_perf.c
	$(CC) -Wall $< $(LDFLAGS) -o $@

clean:
	rm -rf temp *.o *.out $(BINARY) bin/clean_fscache bin/enable_perf

# Delete all gitignored files (but not directories)
distclean: clean
	git clean -f -X; rm -rf ./autom4te.cache

# Autoformat code with clang format
format:
	./scripts/autoformat.sh

test: $(BINARY)
	./scripts/test.sh

valgrind : $(BINARY)
	./scripts/test.sh valgrind
