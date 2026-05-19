# Dilithium (research fork)

This repository is a fork of the upstream Dilithium implementation (pq-crystals/dilithium),
customized for post-quantum cryptography (PQC) research and benchmarking. Dilithium is
standardized as [FIPS 204](https://csrc.nist.gov/pubs/fips/204/final) (ML-DSA).

It contains:

- `ref/`: portable reference C implementation
- `avx2/`: optimized x86_64 implementation using AVX2
- `Dilithium_KAT/`: pre-generated NIST KAT request/response files
- `ref/test/`: TCP client/server demo + stress tool (POSIX/WSL)

For a list of changes in this fork, see [CHANGELOG.md](CHANGELOG.md).

## Table of contents

- [Reproducibility quick start](#reproducibility-quick-start)
- [Build](#build)
- [Correctness tests](#correctness-tests)
- [Benchmarking (cycle counts)](#benchmarking-cycle-counts)
- [Deterministic test vectors](#deterministic-test-vectors)
- [NIST KAT generator (optional)](#nist-kat-generator-optional)
- [TCP client/server demo (optional)](#tcp-clientserver-demo-optional)
- [Coverage (optional)](#coverage-optional)
- [License](#license)

## Reproducibility quick start

### Platform notes

- **Linux is recommended** for reproducible benchmarking.
- **macOS** builds the `ref/` implementation fine in most setups.
- **Windows**: use **WSL2** for the simplest build/run workflow (the stress tool under
	`ref/test/` uses POSIX APIs such as `fork()`).
- `avx2/` requires an x86_64 CPU with **AVX2**.

### Dependencies

Ubuntu/Debian:

```sh
sudo apt-get update
sudo apt-get install -y build-essential make pkg-config libssl-dev
```

Optional tools:

```sh
sudo apt-get install -y valgrind lcov
```

macOS (OpenSSL headers/libs may require flags):

```sh
brew install openssl
export CFLAGS="-I$(brew --prefix openssl)/include"
export NISTFLAGS="-I$(brew --prefix openssl)/include"
export LDFLAGS="-L$(brew --prefix openssl)/lib"
```

## Build

All commands below assume you are at the repository root.

### Reference implementation (`ref/`)

Build correctness tests:

```sh
make -C ref clean
make -C ref all
```

This produces:

- `ref/test/test_dilithium2`
- `ref/test/test_dilithium3`
- `ref/test/test_dilithium5`

### AVX2 implementation (`avx2/`)

Requires an x86_64 CPU with AVX2.

```sh
make -C avx2 clean
make -C avx2 all
```

This produces (among others):

- `avx2/test/test_dilithium2`, `avx2/test/test_vectors2`, `avx2/test/test_speed2`
- `avx2/test/test_dilithium3`, `avx2/test/test_vectors3`, `avx2/test/test_speed3`
- `avx2/test/test_dilithium5`, `avx2/test/test_vectors5`, `avx2/test/test_speed5`

## Correctness tests

Reference:

```sh
./ref/test/test_dilithium2
./ref/test/test_dilithium3
./ref/test/test_dilithium5
```

AVX2:

```sh
./avx2/test/test_dilithium2
./avx2/test/test_dilithium3
./avx2/test/test_dilithium5
```

## Benchmarking (cycle counts)

The `test_speed*` programs print median and average cycle counts (1000 iterations) using
`RDTSC`.

Reference:

```sh
make -C ref speed
./ref/test/test_speed2
./ref/test/test_speed3
./ref/test/test_speed5
```

AVX2:

```sh
make -C avx2 speed
./avx2/test/test_speed2
./avx2/test/test_speed3
./avx2/test/test_speed5
```

Reproducibility tips:

- Pin the exact commit hash: `git rev-parse HEAD`
- Record compiler versions: `gcc --version` / `clang --version`
- Record CPU model and frequency scaling settings (e.g., `lscpu`)

## Deterministic test vectors

The `test_vectors*` programs generate deterministic test vectors (10000 sets). Randomness
is derived from SHAKE128 on empty input (deterministic), and the expected SHA256 sums are
stored in `SHA256SUMS`.

Reference (note: vector binaries are **not** built by default in `ref/`):

```sh
make -C ref test/test_vectors2 test/test_vectors3 test/test_vectors5
./ref/test/test_vectors2 > tvecs2
./ref/test/test_vectors3 > tvecs3
./ref/test/test_vectors5 > tvecs5
sha256sum -c SHA256SUMS
```

macOS:

```sh
shasum -a256 -c SHA256SUMS
```

## NIST KAT generator (optional)

The NIST KAT generator (under `ref/nistkat/`) requires OpenSSL.

```sh
make -C ref nistkat
./ref/nistkat/PQCgenKAT_sign2
./ref/nistkat/PQCgenKAT_sign3
./ref/nistkat/PQCgenKAT_sign5
```

This writes `PQCsignKAT_*.req` / `PQCsignKAT_*.rsp` in the current directory.

The repository also includes pre-generated KAT files under `Dilithium_KAT/`.

## TCP client/server demo (optional)

This fork includes a simple TCP challenge/response demo under `ref/test/`:

- `test_dilithium_keygen{2,3,5}`: generate a keypair and write `client_sk.bin`,
	`client_pk.bin`, `server_pk.bin`
- `test_dilithium_server{2,3,5}`: listen on TCP port `5000`, send a challenge, verify the
	signature
- `test_dilithium_client{2,3,5}`: connect to server, sign the challenge, send the signature
- `test_dilithium_stress{2,3,5}`: concurrent client load generator (uses `fork()`)

Mode mapping:

| Mode suffix | Parameter set |
|-----------:|---------------|
| 2 | Dilithium2 |
| 3 | Dilithium3 |
| 5 | Dilithium5 |

Build:

```sh
make -C ref/test clean
make -C ref/test all
```

Run (localhost, recommended on Linux/WSL):

```sh
make -C ref/test keygen MODE=2
make -C ref/test run-server MODE=2
```

In a second terminal:

```sh
make -C ref/test run-client MODE=2 TARGET_IP=127.0.0.1
```

Stress tool:

```sh
make -C ref/test stress MODE=2 TARGET_IP=127.0.0.1 CONCURRENT_SESSIONS=10
```

Network protocol (framed):

1. server → client: `uint32_be length` + challenge bytes
2. client → server: `uint32_be length` + signature bytes

Logs and files are written in `ref/test/` (e.g., `client.log`, `server.log`, and `*.bin`).
The client log path can be overridden via `CLIENT_LOG_PATH`.

## Coverage (optional)

Generate an lcov report for the `ref/` implementation:

```sh
./runlcov.sh
```

## License

See [LICENSE](LICENSE). This repository is available under CC0 (public domain dedication),
Apache 2.0, or GPL 2.0. The NIST KAT generator sources under `ref/nistkat/` are provided by
NIST.
