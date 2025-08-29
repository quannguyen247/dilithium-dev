/*
This header defines all the cryptographic parameters for Dilithium,
such as the modulus Q, dimensions K and L, and packing sizes.
These parameters are used throughout polyvec.c and poly.h to ensure
correct sizes and arithmetic.
*/
#ifndef PARAMS_H
#define PARAMS_H

#include "config.h"

/*
 * =================================================================================
 * GENERAL PARAMETERS
 * These are the basic building blocks for our cryptographic construction.
 * =================================================================================
 */

// A "seed" is a starting value for generating random numbers.
// SEEDBYTES is the length of this seed in bytes. Think of it like the
// seed you use to generate a specific world in a game like Minecraft.
#define SEEDBYTES 32

// CRH stands for "Collision-Resistant Hash".
// CRHBYTES is the length of it.
#define CRHBYTES 64

// 'tr' stands for "trapdoor", which is the secret information.
#define TRBYTES 64

// RNDBYTES is the number of random bytes used during the signing process
// to ensure that each signature is unique, even if you sign the same
// message twice.
#define RNDBYTES 32

// N is the degree of our polynomials. Polynomials are math expressions with
// variables (like x^2 + 3x - 5). In our case, they are much bigger!
// N=256 means our polynomials have 256 coefficients.
#define N 256

// Q is a very large prime number that we use as a "modulus".
// All our math is done "modulo Q", which means we only care about the
// remainder after dividing by Q.
#define Q 8380417

// D is a parameter used for rounding numbers during our calculations.
// ==> compressing the signature to make it smaller.
#define D 13

// A special number that helps with a super-fast math trick called NTT
// (Number Theoretic Transform), which is used to multiply our polynomials
// very quickly.
#define ROOT_OF_UNITY 1753


/*
 * =================================================================================
 * DILITHIUM SECURITY MODES
 *
 * Dilithium comes in different security levels. We use a preprocessor directive
 * (`#if`) to select the parameters for the chosen level.
 *
 * - K and L: These define the dimensions of our "vectors" and "matrices" of
 *   polynomials. (Bigger numbers = harder puzzle to solve)
 * - ETA: This controls how "small" the secret numbers in our secret key are.
 *   (They must be small to keep the key safe)
 * - TAU: A parameter related to the number of non-zero parts in the "challenge"
 *   polynomial, which is used to ensure the signature is valid.
 * - BETA: This is a limit on the size of the final signature. (more details in our doc file: week 3)
 * - GAMMA1 & GAMMA2: These are used to split numbers into "high bits" and "low bits",
 *   a key trick for making the signature and public key smaller.
 * - OMEGA: The number of "hint" bits in the signature. These hints help the
 *   verifier reconstruct part of the calculation without needing the secret key.
 * - CTILDEBYTES: The size of the seed for the "challenge" polynomial.
 * =================================================================================
 */
#if DILITHIUM_MODE == 2
#define K 4
#define L 4
#define ETA 2
#define TAU 39
#define BETA 78
#define GAMMA1 (1 << 17)
#define GAMMA2 ((Q-1)/88)
#define OMEGA 80
#define CTILDEBYTES 32

#elif DILITHIUM_MODE == 3
#define K 6
#define L 5
#define ETA 4
#define TAU 49
#define BETA 196
#define GAMMA1 (1 << 19)
#define GAMMA2 ((Q-1)/32)
#define OMEGA 55
#define CTILDEBYTES 48

#elif DILITHIUM_MODE == 5
#define K 8
#define L 7
#define ETA 2
#define TAU 60
#define BETA 120
#define GAMMA1 (1 << 19)
#define GAMMA2 ((Q-1)/32)
#define OMEGA 75
#define CTILDEBYTES 64

#endif

/*
 * =================================================================================
 * PACKED SIZES (IN BYTES)
 *
 * Our polynomials are lists of numbers, but to save them in a file or send
 * them over the internet, we need to convert them into a sequence of bytes.
 * This is called "packing". These definitions control how many bytes each
 * type of packed polynomial will take up.
 *
 * The names (T1, T0, Z, W1, ETA) refer to different kinds of polynomials
 * used in the Dilithium algorithm.
 * =================================================================================
 */
#define POLYT1_PACKEDBYTES  320
#define POLYT0_PACKEDBYTES  416
#define POLYVECH_PACKEDBYTES (OMEGA + K)

#if GAMMA1 == (1 << 17)
#define POLYZ_PACKEDBYTES   576
#elif GAMMA1 == (1 << 19)
#define POLYZ_PACKEDBYTES   640
#endif

#if GAMMA2 == (Q-1)/88
#define POLYW1_PACKEDBYTES  192
#elif GAMMA2 == (Q-1)/32
#define POLYW1_PACKEDBYTES  128
#endif

#if ETA == 2
#define POLYETA_PACKEDBYTES  96
#elif ETA == 4
#define POLYETA_PACKEDBYTES 128
#endif

/*
 * =================================================================================
 * FINAL CRYPTO SIZES (IN BYTES)
 *
 * These definitions calculate the total size of the public key, the secret key,
 * and the signature itself by adding up the sizes of all their parts.
 * =================================================================================
 */

// The size of the public key. This is the key you can share with anyone.
// They use it to check that your signature is authentic.
#define CRYPTO_PUBLICKEYBYTES (SEEDBYTES + K*POLYT1_PACKEDBYTES)

// The size of the secret key. You must KEEP THIS SECRET!
// It's used to create the signatures.
#define CRYPTO_SECRETKEYBYTES (2*SEEDBYTES \
                               + TRBYTES \
                               + L*POLYETA_PACKEDBYTES \
                               + K*POLYETA_PACKEDBYTES \
                               + K*POLYT0_PACKEDBYTES)

// The size of the final digital signature. This is the "proof" that you
// send along with your message.
#define CRYPTO_BYTES (CTILDEBYTES + L*POLYZ_PACKEDBYTES + POLYVECH_PACKEDBYTES)

#endif
