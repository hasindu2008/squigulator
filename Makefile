CC       = gcc
CXX      = g++
LANGFLAG = -x c++
CPPFLAGS += -I slow5lib/include/
CFLAGS   += -g -Wall -O2  -std=c99
LDFLAGS  += $(LIBS) -lpthread -lz -rdynamic -lm
BUILD_DIR = build

ifeq ($(zstd),1)
LDFLAGS		+= -lzstd
endif

BINARY = sigsim
OBJ = $(BUILD_DIR)/main.o \
      $(BUILD_DIR)/model.o \
	  $(BUILD_DIR)/genref.o \
	  $(BUILD_DIR)/misc.o \
	  $(BUILD_DIR)/sim.o \

PREFIX = /usr/local
VERSION = `git describe --tags`

ifdef asan
	CFLAGS += -fsanitize=address -fno-omit-frame-pointer
	LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
endif

.PHONY: clean distclean test

$(BINARY): $(OBJ) slow5lib/lib/libslow5.a
	$(CC) $(CFLAGS) $(OBJ) slow5lib/lib/libslow5.a $(LDFLAGS) -o $@

$(BUILD_DIR)/main.o: src/main.c src/misc.h src/error.h  src/sigsim.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/model.o: src/model.c src/model.h  src/misc.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/genref.o: src/genref.c src/ref.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/misc.o: src/misc.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/sim.o: src/sim.c src/ref.h src/misc.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

slow5lib/lib/libslow5.a:
	$(MAKE) -C slow5lib zstd=$(zstd) no_simd=$(no_simd) zstd_local=$(zstd_local)  lib/libslow5.a

clean:
	rm -rf $(BINARY) $(BUILD_DIR)/*.o
	make -C slow5lib clean

# Delete all gitignored files (but not directories)
distclean: clean
	git clean -f -X
	rm -rf $(BUILD_DIR)/* autom4te.cache

test: $(BINARY)
	scripts/test.sh

valgrind: $(BINARY)
	valgrind --leak-check=full ./sigsim test/nCoV-2019.reference.fasta -o a.blow5 -q a.fasta -n 10
