CC ?= gcc
CFLAGS ?= -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wvla -Wpointer-arith -O3 -fomit-frame-pointer

ROOT := ..

SOURCES = $(ROOT)/sign.c $(ROOT)/packing.c $(ROOT)/polyvec.c $(ROOT)/poly.c \
  $(ROOT)/ntt.c $(ROOT)/reduce.c $(ROOT)/rounding.c
HEADERS = $(ROOT)/config.h $(ROOT)/params.h $(ROOT)/api.h $(ROOT)/sign.h \
  $(ROOT)/packing.h $(ROOT)/polyvec.h $(ROOT)/poly.h $(ROOT)/ntt.h \
  $(ROOT)/reduce.h $(ROOT)/rounding.h $(ROOT)/symmetric.h $(ROOT)/randombytes.h
KECCAK_SOURCES = $(SOURCES) $(ROOT)/fips202.c $(ROOT)/symmetric-shake.c
KECCAK_HEADERS = $(HEADERS) $(ROOT)/fips202.h

.PHONY: all clean

all: \
  test_dilithium_client2 \
  test_dilithium_client3 \
  test_dilithium_client5 \
  test_dilithium_server2 \
  test_dilithium_server3 \
  test_dilithium_server5

test_dilithium_client2: test_dilithium_client.c $(ROOT)/randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=2 \
	  -o $@ $< $(ROOT)/randombytes.c $(KECCAK_SOURCES)

test_dilithium_client3: test_dilithium_client.c $(ROOT)/randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=3 \
	  -o $@ $< $(ROOT)/randombytes.c $(KECCAK_SOURCES)

test_dilithium_client5: test_dilithium_client.c $(ROOT)/randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=5 \
	  -o $@ $< $(ROOT)/randombytes.c $(KECCAK_SOURCES)

test_dilithium_server2: test_dilithium_server.c $(ROOT)/randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=2 \
	  -o $@ $< $(ROOT)/randombytes.c $(KECCAK_SOURCES)

test_dilithium_server3: test_dilithium_server.c $(ROOT)/randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=3 \
	  -o $@ $< $(ROOT)/randombytes.c $(KECCAK_SOURCES)

test_dilithium_server5: test_dilithium_server.c $(ROOT)/randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=5 \
	  -o $@ $< $(ROOT)/randombytes.c $(KECCAK_SOURCES)

clean:
	rm -f test_dilithium_client2
	rm -f test_dilithium_client3
	rm -f test_dilithium_client5
	rm -f test_dilithium_server2
	rm -f test_dilithium_server3
	rm -f test_dilithium_server5
