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

.PHONY: all run-server run-client stress keygen clean

MODE ?= 2
TARGET_IP ?= 192.168.4.85
CONCURRENT_SESSIONS ?= 10
BATCHES ?= 0
BATCH_DELAY_SEC ?= 0

CLIENT_BIN := test_dilithium_client$(MODE)
SERVER_BIN := test_dilithium_server$(MODE)
STRESS_BIN := test_dilithium_stress$(MODE)
KEYGEN_BIN := test_dilithium_keygen$(MODE)

all: \
  test_dilithium_client2 \
  test_dilithium_client3 \
  test_dilithium_client5 \
  test_dilithium_server2 \
  test_dilithium_server3 \
	test_dilithium_server5 \
	test_dilithium_stress2 \
	test_dilithium_stress3 \
	test_dilithium_stress5 \
	test_dilithium_keygen2 \
	test_dilithium_keygen3 \
	test_dilithium_keygen5

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

test_dilithium_stress2: test_dilithium_stress.c $(ROOT)/randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=2 \
	  -o $@ $< $(ROOT)/randombytes.c $(KECCAK_SOURCES)

test_dilithium_stress3: test_dilithium_stress.c $(ROOT)/randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=3 \
	  -o $@ $< $(ROOT)/randombytes.c $(KECCAK_SOURCES)

test_dilithium_stress5: test_dilithium_stress.c $(ROOT)/randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=5 \
	  -o $@ $< $(ROOT)/randombytes.c $(KECCAK_SOURCES)

test_dilithium_keygen2: test_dilithium_keygen.c $(ROOT)/randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=2 \
	  -o $@ $< $(ROOT)/randombytes.c $(KECCAK_SOURCES)

test_dilithium_keygen3: test_dilithium_keygen.c $(ROOT)/randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=3 \
	  -o $@ $< $(ROOT)/randombytes.c $(KECCAK_SOURCES)

test_dilithium_keygen5: test_dilithium_keygen.c $(ROOT)/randombytes.c $(KECCAK_SOURCES) \
  $(KECCAK_HEADERS)
	$(CC) $(CFLAGS) -DDILITHIUM_MODE=5 \
	  -o $@ $< $(ROOT)/randombytes.c $(KECCAK_SOURCES)

run-server: $(SERVER_BIN)
	@echo "[RUN] Server MODE=$(MODE) on port 5000"
	@./$(SERVER_BIN)

run-client: $(CLIENT_BIN)
	@echo "[RUN] Client MODE=$(MODE) -> $(TARGET_IP)"
	@./$(CLIENT_BIN) $(TARGET_IP)

stress: $(STRESS_BIN)
	@echo "[STRESS] MODE=$(MODE) TARGET_IP=$(TARGET_IP) CONCURRENT=$(CONCURRENT_SESSIONS)"
	@TARGET_IP=$(TARGET_IP) CONCURRENT_SESSIONS=$(CONCURRENT_SESSIONS) \
	  BATCHES=$(BATCHES) BATCH_DELAY_SEC=$(BATCH_DELAY_SEC) \
	  ./$(STRESS_BIN)

keygen: $(KEYGEN_BIN)
	@echo "[KEYGEN] MODE=$(MODE)"
	@./$(KEYGEN_BIN)

clean:
	rm -f test_dilithium_client2
	rm -f test_dilithium_client3
	rm -f test_dilithium_client5
	rm -f test_dilithium_server2
	rm -f test_dilithium_server3
	rm -f test_dilithium_server5
	rm -f test_dilithium_stress2
	rm -f test_dilithium_stress3
	rm -f test_dilithium_stress5
	rm -f test_dilithium_keygen2
	rm -f test_dilithium_keygen3
	rm -f test_dilithium_keygen5
