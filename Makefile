CC ?= cc
CPPFLAGS ?= -D_FILE_OFFSET_BITS=64 -D_POSIX_C_SOURCE=200809L -D_XOPEN_SOURCE=700
CFLAGS ?= -O2 -g -std=c99 -Wall -Wextra -Wpedantic -Wformat=2
LDFLAGS ?=
THREAD_FLAGS ?= -pthread

PROG := xamdst
SOURCES := xamdst.c config.c input.c intervals.c engine.c report.c util.c
OBJECTS := $(SOURCES:.c=.o)
DEPS := $(OBJECTS:.o=.d)

# Prefer an installed HTSlib. HTSLIB_DIR may point at a source/build tree when
# pkg-config is unavailable.
HTSLIB_DIR ?=
ifneq ($(strip $(MAKECMDGOALS)),clean)
  ifeq ($(strip $(HTSLIB_DIR)),)
    HAVE_HTSLIB := $(shell pkg-config --exists htslib 2>/dev/null && echo yes)
    ifeq ($(HAVE_HTSLIB),yes)
      HTSLIB_CFLAGS := $(shell pkg-config --cflags htslib)
      HTSLIB_LIBS := $(shell pkg-config --libs htslib)
    else
      # Some distributions ship HTSlib headers and libraries without a
      # pkg-config file (or without pkg-config installed).  Fall back to the
      # conventional system include/library search paths before failing.
      HTSLIB_HEADER := $(firstword $(wildcard /usr/include/htslib/sam.h \
                                             /usr/local/include/htslib/sam.h \
                                             /opt/homebrew/include/htslib/sam.h \
                                             /opt/local/include/htslib/sam.h))
      ifeq ($(strip $(HTSLIB_HEADER)),)
        $(error HTSlib was not found. Install HTSlib >= 1.13 or pass HTSLIB_DIR=/path/to/htslib)
      else
        HTSLIB_CFLAGS := -I$(patsubst %/htslib/sam.h,%,$(HTSLIB_HEADER))
        HTSLIB_LIBS := -lhts
      endif
    endif
  else
    HTSLIB_CFLAGS := -I$(HTSLIB_DIR)
    HTSLIB_LIBS := $(HTSLIB_DIR)/libhts.a -lz -lbz2 -llzma -lpthread -lm
  endif
endif

CPPFLAGS += -I. $(HTSLIB_CFLAGS)
LDLIBS += $(HTSLIB_LIBS) -lz -lpthread -lm

.PHONY: all clean test sanitize tsan benchmark

all: $(PROG)

$(PROG): $(OBJECTS)
	$(CC) $(CFLAGS) $(THREAD_FLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $(THREAD_FLAGS) -MMD -MP -c $< -o $@

test: $(PROG)
	sh ./tests/run_tests.sh ./$(PROG)

sanitize:
	$(MAKE) clean
	$(MAKE) CFLAGS="$(CFLAGS) -O1 -fno-omit-frame-pointer -fsanitize=address,undefined" LDFLAGS="$(LDFLAGS) -fsanitize=address,undefined" $(PROG)
	sh ./tests/run_tests.sh ./$(PROG)
	$(MAKE) clean

tsan:
	$(MAKE) clean
	$(MAKE) CFLAGS="$(CFLAGS) -O1 -fno-omit-frame-pointer -fsanitize=thread" LDFLAGS="$(LDFLAGS) -fsanitize=thread" $(PROG)
	sh ./tests/run_tests.sh ./$(PROG)
	$(MAKE) clean

benchmark: $(PROG)
	sh ./tests/benchmark.sh ./$(PROG)

clean:
	rm -f $(PROG) $(OBJECTS) $(DEPS)

-include $(DEPS)
