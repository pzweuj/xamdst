CC=		gcc
CFLAGS=		-g -Wall -O2 -std=gnu99 
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE

# htslib support - use pkg-config if available, otherwise use default paths
HTSLIB_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null || echo "")
HTSLIB_LIBS := $(shell pkg-config --libs htslib 2>/dev/null || echo "-lhts")

# Objects for output file writing (bgzf is now from htslib)
LOBJS=		kstring.o bedutil.o commons.o
PROG=		bamdst
INCLUDES=	-I. $(HTSLIB_CFLAGS)
LIBPATH=        -L. 

.SUFFIXES:.c .o
.PHONY: all

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:clean $(PROG) 

.PHONY:all clean

lib:libbam.a

libbam.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

bamdst:lib
		$(CC) $(CFLAGS) -o $@ $(LDFLAGS) bamdst.c $(LIBPATH) $(INCLUDES) -lm $(HTSLIB_LIBS) -lz -lpthread

kstring.o:kstring.c kstring.h
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) kstring.c -o $@

commons.o:commons.c commons.h
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) commons.c -o $@	

bedutil.o:bedutil.c bedutil.h
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) bedutil.c -o $@	

clean:
		rm -fr gmon.out *.o a.out *.exe *.dSYM  $(PROG) *~ *.a target.dep *.plot *.report *.tsv.gz uncover.bed

