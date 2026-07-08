CC=		gcc
CFLAGS=		-g -Wall -O2 -std=gnu99
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE

# Source files
SOURCES=	xamdst.c kstring.c bedutil.c commons.c
PROG=		xamdst

# ---- htslib configuration ---------------------------------------------------
# Resolution order:
#   1. HTSLIB_DIR=<path>     explicit override (a built htslib source tree)
#   2. system htslib         detected via pkg-config
#   3. sibling ./htslib      local source tree -- libhts.a is auto-built if missing
#   4. otherwise             error
#
# If static linking fails with "cannot find -lbz2/-llzma", your htslib was
# built without them -- override:
#       make HTSLIB_DEPS="-lz -lpthread -lm"
HTSLIB_DIR ?=

ifeq ($(HTSLIB_DIR),)
  HAVE_PKG_HTSLIB := $(shell pkg-config --exists htslib 2>/dev/null && echo yes)
  ifneq ($(HAVE_PKG_HTSLIB),)
    # 2. system htslib (pkg-config)
    HTSLIB_CFLAGS := $(shell pkg-config --cflags htslib)
    HTSLIB_LIBS   := $(shell pkg-config --libs htslib)
  else
    ifneq ($(wildcard htslib/htslib/hts.h),)
      # 3. sibling ./htslib source tree (libhts.a auto-built below)
      HTSLIB_CFLAGS := -Ihtslib
      HTSLIB_LIBS   := htslib/libhts.a
      HTSLIB_DEPS   ?= -lz -lbz2 -llzma -lpthread -lm
      HTSLIB_BUILD  := htslib/libhts.a
    else
      # 4. not found
      $(error htslib not found. \
Install htslib system-wide (then pkg-config detects it), \
or clone it as a sibling: `git clone https://github.com/samtools/htslib htslib`, \
or pass HTSLIB_DIR=/path/to/htslib)
    endif
  endif
else
  # 1. explicit override
  HTSLIB_CFLAGS := -I$(HTSLIB_DIR)
  HTSLIB_LIBS   := $(HTSLIB_DIR)/libhts.a
  HTSLIB_DEPS   ?= -lz -lbz2 -llzma -lpthread -lm
  HTSLIB_BUILD  := $(HTSLIB_DIR)/libhts.a
endif

INCLUDES=	-I. $(HTSLIB_CFLAGS)
LIBS=		$(HTSLIB_LIBS) $(HTSLIB_DEPS) -lz -lpthread -lm

.SUFFIXES:.c .o
.PHONY: all clean

# auto-build the local htslib (handles a git clone: autoreconf + configure + make)
ifdef HTSLIB_BUILD
$(HTSLIB_BUILD):
		cd $(dir $(HTSLIB_BUILD)) && \
		(test -f configure || autoreconf -i) && \
		(test -f config.status || ./configure) && \
		$(MAKE)
endif

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:clean $(PROG)

xamdst: $(SOURCES:.c=.o) $(HTSLIB_BUILD)
		$(CC) $(CFLAGS) -o $@ $(SOURCES:.c=.o) $(INCLUDES) $(LIBS)

clean:
		rm -fr gmon.out *.o a.out *.exe *.dSYM  $(PROG) *~ *.a target.dep *.plot *.report *.tsv.gz uncover.bed
