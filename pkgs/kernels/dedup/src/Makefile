# Makefile for dedup kernel

PREFIX=${PARSECDIR}/pkgs/kernels/dedup/inst/${PARSECPLAT}

TARGET=dedup

CFLAGS += -Wall -fno-strict-aliasing -D_XOPEN_SOURCE=600

ostype=$(findstring solaris, ${PARSECPLAT})

ifeq "$(ostype)" "solaris"
    CFLAGS += -std=gnu99
endif

LIBS += -lm

DEDUP_OBJ = hashtable.o util.o rabin.o decoder.o mbuffer.o sha.o

# Uncomment the following to enable gzip compression
CFLAGS += -DENABLE_GZIP_COMPRESSION
LIBS += -lz

# Uncomment the following to enable bzip2 compression
#CFLAGS += -DENABLE_BZIP2_COMPRESSION
#LIBS += -lbz2

COMPILER=$(CC)
COMPILER_FLAGS=$(CFLAGS)

ifdef version
  ifeq "$(version)" "pthreads"
    CFLAGS += -DENABLE_PTHREADS -DPARALLEL -pthread
    DEDUP_OBJ += dedup.o encoder.o queue.o binheap.o tree.o
  else
  ifeq "$(version)" "tbb_funcobj"
    CFLAGS += -DPARALLEL -DTBB_FUNCOBJ -DTBB -DBUDDY_REORDER
    CXXFLAGS += -DPARALLEL -DTBB_FUNCOBJ -DTBB -DBUDDY_REORDER
    COMPILER=$(CXX)
    COMPILER_FLAGS=$(CXXFLAGS)
    DEDUP_OBJ += dedup_swan.o encoder-funcobj.o
  else
  ifeq "$(version)" "swan-obj"
    CFLAGS += -DENABLE_SWAN -pthread
    CXXFLAGS += -DENABLE_SWAN -pthread
    DEDUP_OBJ := dedup_swan.o encoder_lists.o $(DEDUP_OBJ)
    encoder_lists.o: CC="$(CXX) -I/home/hvandierendonck/research/parsec-3.0/bin/../pkgs/libs/ssl/inst/amd64-linux.swan-obj/include -I/home/hvandierendonck/research/parsec-3.0/bin/../pkgs/libs/zlib/inst/amd64-linux.swan-obj/include -I/home/hvandierendonck/research/swan-queue -I/home/hvandierendonck/research/swan-queue/swan -I/home/hvandierendonck/research/swan-queue/obj -I/home/hvandierendonck/research/swan-queue/obj/scheduler $(CXXFLAGS) $(CFLAGS)"
    COMPILER=$(CXX)
    COMPILER_FLAGS=$(CXXFLAGS)
  else
    DEDUP_OBJ += encoder.o dedup.o
  endif
  endif
  endif
endif

CFLAGS+=-g

#CFLAGS+=-pg -p
#COMPILER_FLAGS+=-pg -p

all: $(TARGET)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

.cpp.o:
	$(CXX) -c $(CFLAGS) $(CXXFLAGS) $< -o $@

$(TARGET): $(DEDUP_OBJ)
	$(COMPILER) $(COMPILER_FLAGS) -o $(TARGET) $(DEDUP_OBJ) $(LDFLAGS) $(LIBS)

clean:
	rm -f *~ *.o $(TARGET)

install:
	mkdir -p $(PREFIX)/bin
	cp -f $(TARGET) $(PREFIX)/bin/$(TARGET)

