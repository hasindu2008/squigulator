CC       = gcc
CXX      = g++
CPPFLAGS += -I slow5lib/include/
CFLAGS   += -g -Wall -O2  -std=c99
LDFLAGS  += $(LIBS) -lpthread -lz -rdynamic -lm
BUILD_DIR = build

ifeq ($(zstd),1)
LDFLAGS		+= -lzstd
endif

BINARY = squigulator
OBJ = $(BUILD_DIR)/main.o \
      $(BUILD_DIR)/model.o \
	  $(BUILD_DIR)/methmodel.o \
	  $(BUILD_DIR)/misc.o \
	  $(BUILD_DIR)/sim.o \
	  $(BUILD_DIR)/thread.o \
	  $(BUILD_DIR)/format.o \
	  $(BUILD_DIR)/gensig.o \
	  $(BUILD_DIR)/genread.o \
	  $(BUILD_DIR)/ref.o \

PREFIX = /usr/local
VERSION = `git describe --tags`

ifdef asan
	CFLAGS += -fsanitize=address -fno-omit-frame-pointer
	LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
endif

.PHONY: clean distclean test

$(BINARY): $(OBJ) slow5lib/lib/libslow5.a
	$(CC) $(CFLAGS) $(OBJ) slow5lib/lib/libslow5.a $(LDFLAGS) -o $@

HEADERS = src/error.h src/format.h  src/misc.h  src/model.h  \
		  src/rand.h  src/ref.h  src/seq.h  src/sq.h  src/str.h  src/version.h \
		  src/kseq.h src/khash.h src/ksort.h

$(BUILD_DIR)/main.o: src/main.c $(HEADERS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/model.o: src/model.c src/model.h  src/misc.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/methmodel.o: src/methmodel.c src/model.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/thread.o: src/thread.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/misc.o: src/misc.c src/error.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/sim.o: src/sim.c $(HEADERS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/format.o: src/format.c $(HEADERS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/gensig.o: src/gensig.c $(HEADERS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/genread.o: src/genread.c $(HEADERS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $< -c -o $@

$(BUILD_DIR)/ref.o: src/ref.c $(HEADERS)
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

release: distclean
# make the release
	mkdir -p squigulator-$(VERSION)
	mkdir -p squigulator-$(VERSION)/scripts squigulator-$(VERSION)/docs squigulator-$(VERSION)/slow5lib
	cp -r README.md LICENSE Makefile build src squigulator-$(VERSION)
	cp -r docs/man.md docs/output.md squigulator-$(VERSION)/docs/
	cp -r slow5lib/lib slow5lib/include slow5lib/src  slow5lib/Makefile slow5lib/LICENSE slow5lib/thirdparty/ squigulator-$(VERSION)/slow5lib
	tar -zcf squigulator-$(VERSION)-release.tar.gz squigulator-$(VERSION)
	rm -rf squigulator-$(VERSION)
# make the binaries
	make -j8
	mkdir -p squigulator-$(VERSION)
	mkdir squigulator-$(VERSION)/docs squigulator-$(VERSION)/scripts
	mv squigulator squigulator-$(VERSION)/
	cp -r README.md LICENSE squigulator-$(VERSION)/
	cp -r docs/man.md docs/output.md squigulator-$(VERSION)/docs/
	tar -zcf squigulator-$(VERSION)-x86_64-linux-binaries.tar.gz squigulator-$(VERSION)
	rm -rf squigulator-$(VERSION)
	tar xf squigulator-$(VERSION)-x86_64-linux-binaries.tar.gz
	mv squigulator-$(VERSION)/squigulator squigulator
	scripts/test.sh

install: $(BINARY)
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp -f $(BINARY) $(DESTDIR)$(PREFIX)/bin

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/$(BINARY)

test: $(BINARY)
	scripts/test.sh

valgrind: $(BINARY)
	valgrind --leak-check=full ./squigulator test/nCoV-2019.reference.fasta -o a.blow5 -q a.fasta -n 10
