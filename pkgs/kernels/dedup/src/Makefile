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

ifdef version
  ifeq "$(version)" "pthreads"
    CFLAGS += -DENABLE_PTHREADS -pthread
    DEDUP_OBJ += queue.o binheap.o tree.o
  endif
  ifeq "$(version)" "swan-obj"
    #CFLAGS += ${dedup_flags}
    CFLAGS += -DENABLE_SWAN -pthread
    DEDUP_OBJ := dedup_swan.o encoder_lists.o $(DEDUP_OBJ)
    encoder_lists.o: CC="$(CXX) -I/home/hvandierendonck/research/parsec-3.0/bin/../pkgs/libs/ssl/inst/amd64-linux.swan-obj/include -I/home/hvandierendonck/research/parsec-3.0/bin/../pkgs/libs/zlib/inst/amd64-linux.swan-obj/include -I/home/hvandierendonck/research/swan-queue -I/home/hvandierendonck/research/swan-queue/swan -I/home/hvandierendonck/research/swan-queue/obj -I/home/hvandierendonck/research/swan-queue/obj/scheduler $(CXXFLAGS) $(CFLAGS)"
  else
    DEDUP_OBJ += encoder.o dedup.o
  endif
endif

CFLAGS+=-g

all: $(TARGET)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

.cpp.o:
	$(CXX) -c $(CFLAGS) $(CXXFLAGS) $< -o $@

ifeq "$(version)" "swan-obj"
$(TARGET): $(DEDUP_OBJ)
	$(CXX) $(CFLAGS) -o $(TARGET) $(DEDUP_OBJ) $(LDFLAGS) $(LIBS)
else
$(TARGET): $(DEDUP_OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(DEDUP_OBJ) $(LIBS)
endif

clean:
	rm -f *~ *.o $(TARGET)

install:
	mkdir -p $(PREFIX)/bin
	cp -f $(TARGET) $(PREFIX)/bin/$(TARGET)

