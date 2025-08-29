#ifndef POLY_H
#define POLY_H
/*
 * This file defines what a `poly` is and lists all the operations we can do
 * on them, like adding, subtracting, and multiplying them.
 * =================================================================================
 */
#include <stdint.h>
#include "params.h"

/*
 * This defines the `poly` structure. It's a simple container that holds an
 * array (a list) of `N` 32-bit integers. This is our "polynomial".
 */
typedef struct {
  int32_t coeffs[N];
} poly;

/*
 * A NOTE ON `DILITHIUM_NAMESPACE`:
 * You'll see this macro everywhere. It's a clever trick to add a prefix
 * to all our function names. This helps prevent name clashes if we were to
 * combine this code with another crypto library that might also have a
 * function named `poly_add`.
 */

/*
 * =================================================================================
 * BASIC ARITHMETIC AND UTILITY FUNCTIONS
 * =================================================================================
 */

// Reduces the coefficients of a polynomial to a standard range around zero.
// This is a cleanup step after doing calculations.
#define poly_reduce DILITHIUM_NAMESPACE(poly_reduce)
void poly_reduce(poly *a);

// If a coefficient is negative, this adds Q to it to make it positive.
// This ensures all coefficients are in the range [0, Q-1].
#define poly_caddq DILITHIUM_NAMESPACE(poly_caddq)
void poly_caddq(poly *a);

// Adds two polynomials together (c = a + b).
#define poly_add DILITHIUM_NAMESPACE(poly_add)
void poly_add(poly *c, const poly *a, const poly *b);

// Subtracts one polynomial from another (c = a - b).
#define poly_sub DILITHIUM_NAMESPACE(poly_sub)
void poly_sub(poly *c, const poly *a, const poly *b);

// "Shifts left", which is a fast way to multiply all coefficients by 2^D.
#define poly_shiftl DILITHIUM_NAMESPACE(poly_shiftl)
void poly_shiftl(poly *a);


/*
 * =================================================================================
 * NUMBER THEORETIC TRANSFORM (NTT) FUNCTIONS
 * NTT is a super-fast algorithm for multiplying polynomials.
 * =================================================================================
 */

// Applies the forward NTT to a polynomial, moving it into the "NTT domain".
#define poly_ntt DILITHIUM_NAMESPACE(poly_ntt)
void poly_ntt(poly *a);

// Applies the inverse NTT, bringing a polynomial back from the NTT domain.
#define poly_invntt_tomont DILITHIUM_NAMESPACE(poly_invntt_tomont)
void poly_invntt_tomont(poly *a);

// Multiplies two polynomials together while they are in the NTT domain.
// This is much faster than regular polynomial multiplication.
#define poly_pointwise_montgomery DILITHIUM_NAMESPACE(poly_pointwise_montgomery)
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b);


/*
 * =================================================================================
 * DECOMPOSITION AND HINT FUNCTIONS
 * These are used for compressing data and helping the verifier.
 * =================================================================================
 */

// Splits a polynomial `a` into a "high part" `a1` and a "low part" `a0`.
// This is used to compress the public key.
#define poly_power2round DILITHIUM_NAMESPACE(poly_power2round)
void poly_power2round(poly *a1, poly *a0, const poly *a);

// Another way to split a polynomial into high and low parts.
// This is used during signature verification.
#define poly_decompose DILITHIUM_NAMESPACE(poly_decompose)
void poly_decompose(poly *a1, poly *a0, const poly *a);

// Creates a "hint" polynomial `h` that tells the verifier about "carries"
// that happened during a calculation.
#define poly_make_hint DILITHIUM_NAMESPACE(poly_make_hint)
unsigned int poly_make_hint(poly *h, const poly *a0, const poly *a1);

// Uses the hint `h` to correct a polynomial `a` during verification.
#define poly_use_hint DILITHIUM_NAMESPACE(poly_use_hint)
void poly_use_hint(poly *b, const poly *a, const poly *h);


/*
 * =================================================================================
 * SAMPLING AND SECURITY FUNCTIONS
 * These are for generating random polynomials and for security checks.
 * =================================================================================
 */

// Checks if the largest coefficient in a polynomial is within a security bound `B`.
#define poly_chknorm DILITHIUM_NAMESPACE(poly_chknorm)
int poly_chknorm(const poly *a, int32_t B);

// Generates a polynomial with uniformly random coefficients from a seed.
// This is used to create the public matrix 'A'.
#define poly_uniform DILITHIUM_NAMESPACE(poly_uniform)
void poly_uniform(poly *a,
                  const uint8_t seed[SEEDBYTES],
                  uint16_t nonce);

// Generates a polynomial with small, random coefficients (range is `ETA`).
// This is used to create the secret key polynomials.
#define poly_uniform_eta DILITHIUM_NAMESPACE(poly_uniform_eta)
void poly_uniform_eta(poly *a,
                      const uint8_t seed[CRHBYTES],
                      uint16_t nonce);

// Generates a polynomial with random coefficients in a larger range (`GAMMA1`).
// This is used for the masking polynomial `y` during signing.
#define poly_uniform_gamma1 DILITHIUM_NAMESPACE(poly_uniform_gamma1)
void poly_uniform_gamma1(poly *a,
                         const uint8_t seed[CRHBYTES],
                         uint16_t nonce);

// Generates the "challenge" polynomial `c` from a seed. The challenge is a
// key part of the signature protocol.
#define poly_challenge DILITHIUM_NAMESPACE(poly_challenge)
void poly_challenge(poly *c, const uint8_t seed[CTILDEBYTES]);


/*
 * =================================================================================
 * PACKING AND UNPACKING FUNCTIONS
 * These functions convert polynomials to and from byte arrays for storage or sending.
 * =================================================================================
 */

// Pack a small-coefficient polynomial (`eta`) into a byte array.
#define polyeta_pack DILITHIUM_NAMESPACE(polyeta_pack)
void polyeta_pack(uint8_t *r, const poly *a);
// Unpack a byte array back into a small-coefficient polynomial.
#define polyeta_unpack DILITHIUM_NAMESPACE(polyeta_unpack)
void polyeta_unpack(poly *r, const uint8_t *a);

// Pack the high-bits polynomial `t1` (from the public key).
#define polyt1_pack DILITHIUM_NAMESPACE(polyt1_pack)
void polyt1_pack(uint8_t *r, const poly *a);
// Unpack the high-bits polynomial `t1`.
#define polyt1_unpack DILITHIUM_NAMESPACE(polyt1_unpack)
void polyt1_unpack(poly *r, const uint8_t *a);

// Pack the low-bits polynomial `t0` (from the secret key).
#define polyt0_pack DILITHIUM_NAMESPACE(polyt0_pack)
void polyt0_pack(uint8_t *r, const poly *a);
// Unpack the low-bits polynomial `t0`.
#define polyt0_unpack DILITHIUM_NAMESPACE(polyt0_unpack)
void polyt0_unpack(poly *r, const uint8_t *a);

// Pack the signature polynomial `z`.
#define polyz_pack DILITHIUM_NAMESPACE(polyz_pack)
void polyz_pack(uint8_t *r, const poly *a);
// Unpack the signature polynomial `z`.
#define polyz_unpack DILITHIUM_NAMESPACE(polyz_unpack)
void polyz_unpack(poly *r, const uint8_t *a);

// Pack the high-bits polynomial `w1` (part of the signature).
#define polyw1_pack DILITHIUM_NAMESPACE(polyw1_pack)
void polyw1_pack(uint8_t *r, const poly *a);

#endif
