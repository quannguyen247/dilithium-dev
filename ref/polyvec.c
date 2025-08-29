/*
This file implements functions for handling vectors of polynomials (we use vector to store the coefficients of the polynomial),
It provides operations for two different lengths, `L` and `K`. L (polyvecl) and K (polyveck),
such as +, -, mod, NTT transforms(allows fast computations on polynomials modulo (p)), and matrix expansion.

Key usages:
  Implements arithmetic and transformations for polynomial vectors.
  Used in signature generation and verification routines [step 2 and 3 in the our flow].
  Relies on lower-level polynomial operations from poly.h.
*/

#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"

/*************************************************
* Name:        polyvec_matrix_expand
*
* Description: Implementation of ExpandA. Generates matrix A with uniformly
*              random coefficients a_{i,j} by performing rejection
*              sampling on the output stream of SHAKE128(rho|j|i)
*
* Arguments:   - polyvecl mat[K]: output matrix
*              - const uint8_t rho[]: byte array containing seed rho
**************************************************/
/*
 * This function builds the public matrix 'A', which is a core part of the public key.
 * It's a grid (matrix) of KxL polynomials.
 * It uses the seed `rho` to generate each polynomial `A[i][j]` in a predictable
 * but random-looking way. This way, anyone with the seed `rho` can generate the
 * exact same matrix 'A'.
 * It calls `poly_uniform` from `poly.c` for each of the K*L polynomials.
 */
void polyvec_matrix_expand(polyvecl mat[K], const uint8_t rho[SEEDBYTES]) {
  unsigned int i, j;
  //The two loops ensure every polynomial in the matrix is initialized.
  for(i = 0; i < K; ++i)
    for(j = 0; j < L; ++j)
      poly_uniform(&mat[i].vec[j], rho, (i << 8) + j);
}

/*
 * This function performs a matrix-vector multiplication: t = A * v.
 * 'A' is the public matrix, and 'v' is a vector of polynomials.
 * The result 't' is another vector of polynomials.
 * This is a fundamental operation in the Dilithium algorithm, used during signing.
 * It calls `polyvecl_pointwise_acc_montgomery` to compute each element of the
 * resulting vector 't'. This is essentially a dot product.
 */
void polyvec_matrix_pointwise_montgomery(polyveck *t, const polyvecl mat[K], const polyvecl *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    polyvecl_pointwise_acc_montgomery(&t->vec[i], &mat[i], v);
}

/**************************************************************/
/************ Vectors of polynomials of length L **************/
/**************************************************************/
/*
 * The functions in this section operate on `polyvecl`, which is a vector
 * containing `L` polynomials. These are used for the secret key parts `s1` and `s2`.
 */

/*
 * Generates a vector of `L` polynomials where the coefficients are small and
 * centered around zero. The range of the coefficients is determined by `ETA`.
 * This is used to generate the secret key vectors `s1` and `s2`, which must
 * contain small, secret numbers.
 * `nonce` ensures that we get a different vector each time we call this function.
 */
void polyvecl_uniform_eta(polyvecl *v, const uint8_t seed[CRHBYTES], uint16_t nonce) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_uniform_eta(&v->vec[i], seed, nonce++);
}

/*
 * Generates a vector of `L` polynomials with uniformly random coefficients.
 * This is used to create the "masking" vector `y` during the signing process.
 * The coefficients are in a much larger range defined by `GAMMA1`.
 * This helps to hide the secret key during signing.
 */
void polyvecl_uniform_gamma1(polyvecl *v, const uint8_t seed[CRHBYTES], uint16_t nonce) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_uniform_gamma1(&v->vec[i], seed, L*nonce + i);
}

/*
 * After doing calculations, the coefficients of our polynomials can get very large.
 * This function reduces them to a standard range around zero. It's like a cleanup step.
 * It calls `poly_reduce` for each polynomial in the vector.
 */
void polyvecl_reduce(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_reduce(&v->vec[i]);
}

/*************************************************
* Name:        polyvecl_add
*
* Description: Add vectors of polynomials of length L.
*              No modular reduction is performed.
*
* Arguments:   - polyvecl *w: pointer to output vector
*              - const polyvecl *u: pointer to first summand
*              - const polyvecl *v: pointer to second summand
**************************************************/
void polyvecl_add(polyvecl *w, const polyvecl *u, const polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyvecl_ntt
*
* Description: Forward NTT of all polynomials in vector of length L. Output
*              coefficients can be up to 16*Q larger than input coefficients.
*
* Arguments:   - polyvecl *v: pointer to input/output vector
**************************************************/
/*
 * NTT is a super-fast algorithm for multiplying polynomials. We transform the
 * polynomials into the "NTT domain" before multiplying them.
 */
void polyvecl_ntt(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_ntt(&v->vec[i]);
}

/*
 * This is the reverse of the NTT. It takes polynomials from the NTT domain
 * back to the normal coefficient representation.
 * The `_tomont` part means it also helps with Montgomery reduction, a trick for
 * efficient modular arithmetic.
 */
void polyvecl_invntt_tomont(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_invntt_tomont(&v->vec[i]);
}

/*
 * Multiplies each polynomial in vector `v` by a single polynomial `a`.
 * The result is stored in vector `r`. So, r[i] = a * v[i].
 * This is used when multiplying the matrix `A` by the challenge polynomial `c`.
 */
void polyvecl_pointwise_poly_montgomery(polyvecl *r, const poly *a, const polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_pointwise_montgomery(&r->vec[i], a, &v->vec[i]);
}

/*************************************************
* Name:        polyvecl_pointwise_acc_montgomery
*
* Description: This function computes a dot product of two polynomial vectors: w = u[0]*v[0] + u[1]*v[1] + ...
* It multiplies the polynomials element-wise (`pointwise`) and then adds (`accumulates`)
* all the results into a single output polynomial `w`.
* This is used to calculate `t = A * s1`, where `A` is a matrix and `s1` is a vector.
*
* Arguments:   - poly *w: output polynomial
*              - const polyvecl *u: pointer to first input vector
*              - const polyvecl *v: pointer to second input vector
**************************************************/
void polyvecl_pointwise_acc_montgomery(poly *w,
                                       const polyvecl *u,
                                       const polyvecl *v)
{
  unsigned int i;
  poly t;

  poly_pointwise_montgomery(w, &u->vec[0], &v->vec[0]);
  for(i = 1; i < L; ++i) {
    poly_pointwise_montgomery(&t, &u->vec[i], &v->vec[i]);
    poly_add(w, w, &t);
  }
}

/*************************************************
* Name:        polyvecl_chknorm
*
* Description: This is a security check. It checks if the largest coefficient (in absolute value)
* in any of the polynomials in the vector `v` is smaller than a given `bound`.
* In Dilithium, the signature is only valid if the final `z` vector has a small norm.
* This function is used to verify that condition.
*
* Arguments:   - const polyvecl *v: pointer to vector
*              - int32_t B: norm bound
*
* Returns 0 if norm of all polynomials is strictly smaller than B <= (Q-1)/8
* and 1 otherwise.
**************************************************/
int polyvecl_chknorm(const polyvecl *v, int32_t bound)  {
  unsigned int i;

  for(i = 0; i < L; ++i)
    if(poly_chknorm(&v->vec[i], bound))
      return 1;

  return 0;
}

/**************************************************************/
/************ Vectors of polynomials of length K **************/
/**************************************************************/
/*
 * The functions in this section are very similar to the `polyvecl` functions,
 * but they operate on `polyveck`, which is a vector containing `K` polynomials.
 * These are used for the public key part `t1` and secret key part `s2`.
 */

/*
 * Generates a vector of `K` polynomials with small coefficients (controlled by `ETA`).
 * This is used for generating the secret key vector `s2`.
 */
void polyveck_uniform_eta(polyveck *v, const uint8_t seed[CRHBYTES], uint16_t nonce) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_uniform_eta(&v->vec[i], seed, nonce++);
}

/*************************************************
* Name:        polyveck_reduce
*
* Description: Reduce coefficients of polynomials in vector of length K
*              to representatives in [-6283008,6283008].
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
/*
 * Same as `polyvecl_reduce`, but for a vector of length K.
 * It brings all polynomial coefficients into a standard, smaller range.
 */
void polyveck_reduce(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_reduce(&v->vec[i]);
}

/*************************************************
* Name:        polyveck_caddq
*
* Description: For all coefficients of polynomials in vector of length K
*              add Q if coefficient is negative.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
/*
 * For every coefficient in the vector, if it's negative, it adds Q to make it
 * positive. This ensures all coefficients are in the range [0, Q-1].
 * This is useful before packing the data to be stored or sent.
 */
void polyveck_caddq(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_caddq(&v->vec[i]);
}

/*************************************************
* Name:        polyveck_add
*
* Description: Add vectors of polynomials of length K.
*              No modular reduction is performed.
*
* Arguments:   - polyveck *w: pointer to output vector
*              - const polyveck *u: pointer to first summand
*              - const polyveck *v: pointer to second summand
**************************************************/
/*
 * Adds two vectors of length K, just like `polyvecl_add`.
 * w[i] = u[i] + v[i] for all i.
 */
void polyveck_add(polyveck *w, const polyveck *u, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyveck_sub
*
* Description: Subtract vectors of polynomials of length K.
*              No modular reduction is performed.
*
* Arguments:   - polyveck *w: pointer to output vector
*              - const polyveck *u: pointer to first input vector
*              - const polyveck *v: pointer to second input vector to be
*                                   subtracted from first input vector
**************************************************/
/*
 * Subtracts vector `v` from vector `u`.
 * w[i] = u[i] - v[i] for all i.
 */
void polyveck_sub(polyveck *w, const polyveck *u, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_sub(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyveck_shiftl
*
* Description: Multiply vector of polynomials of Length K by 2^D without modular
*              reduction. Assumes input coefficients to be less than 2^{31-D}.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
/*
 * "shiftl" means "shift left". This is a fast way to multiply by a power of 2.
 * This function multiplies every coefficient in the vector by 2^D.
 * This is part of the calculation `t = A*s1 + s2`.
 */
void polyveck_shiftl(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_shiftl(&v->vec[i]);
}

/*************************************************
* Name:        polyveck_ntt
*
* Description: Forward NTT of all polynomials in vector of length K. Output
*              coefficients can be up to 16*Q larger than input coefficients.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
/*
 * Applies the NTT to every polynomial in a vector of length K.
 */
void polyveck_ntt(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_ntt(&v->vec[i]);
}

/*************************************************
* Name:        polyveck_invntt_tomont
*
* Description: Inverse NTT and multiplication by 2^{32} of polynomials
*              in vector of length K. Input coefficients need to be less
*              than 2*Q.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
/*
 * Applies the inverse NTT to every polynomial in a vector of length K.
 */
void polyveck_invntt_tomont(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_invntt_tomont(&v->vec[i]);
}

/*
 * Multiplies each polynomial in vector `v` by a single polynomial `a`.
 * Same as `polyvecl_pointwise_poly_montgomery` but for a vector of length K.
 */
void polyveck_pointwise_poly_montgomery(polyveck *r, const poly *a, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_pointwise_montgomery(&r->vec[i], a, &v->vec[i]);
}


/*************************************************
* Name:        polyveck_chknorm
*
* Description: Check infinity norm of polynomials in vector of length K.
*              Assumes input polyveck to be reduced by polyveck_reduce().
*
* Arguments:   - const polyveck *v: pointer to vector
*              - int32_t B: norm bound
*
* Returns 0 if norm of all polynomials are strictly smaller than B <= (Q-1)/8
* and 1 otherwise.
**************************************************/
/*
 * Same as `polyvecl_chknorm`, but for a vector of length K.
 * It checks if all coefficients are within a certain security bound.
 */
int polyveck_chknorm(const polyveck *v, int32_t bound) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    if(poly_chknorm(&v->vec[i], bound))
      return 1;

  return 0;
}

/*************************************************
* Name:        polyveck_power2round
*
* Description: For all coefficients a of polynomials in vector of length K,
*              compute a0, a1 such that a mod^+ Q = a1*2^D + a0
*              with -2^{D-1} < a0 <= 2^{D-1}. Assumes coefficients to be
*              standard representatives.
*
* Arguments:   - polyveck *v1: pointer to output vector of polynomials with
*                              coefficients a1
*              - polyveck *v0: pointer to output vector of polynomials with
*                              coefficients a0
*              - const polyveck *v: pointer to input vector
**************************************************/
/*
 * This function is part of compressing the public key.
 * It splits each coefficient `a` in the input vector `v` into two parts:
 * a "high part" `a1` and a "low part" `a0`.
 * `a = a1 * 2^D + a0`.
 * Only the high part `a1` is published in the public key, which saves space.
 * The low part `a0` is discarded.
 */
void polyveck_power2round(polyveck *v1, polyveck *v0, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_power2round(&v1->vec[i], &v0->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyveck_decompose
*
* Description: For all coefficients a of polynomials in vector of length K,
*              compute high and low bits a0, a1 such a mod^+ Q = a1*ALPHA + a0
*              with -ALPHA/2 < a0 <= ALPHA/2 except a1 = (Q-1)/ALPHA where we
*              set a1 = 0 and -ALPHA/2 <= a0 = a mod Q - Q < 0.
*              Assumes coefficients to be standard representatives.
*
* Arguments:   - polyveck *v1: pointer to output vector of polynomials with
*                              coefficients a1
*              - polyveck *v0: pointer to output vector of polynomials with
*                              coefficients a0
*              - const polyveck *v: pointer to input vector
**************************************************/
/*
 * This is another decomposition function, used during signature verification.
 * It splits a polynomial vector `v` into a high-bits vector `v1` and a
 * low-bits vector `v0`. This is used to reconstruct the high bits of `w`
 * during verification to see if they match the public key.
 */
void polyveck_decompose(polyveck *v1, polyveck *v0, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_decompose(&v1->vec[i], &v0->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyveck_make_hint
*
* Description: Compute hint vector.
*
* Arguments:   - polyveck *h: pointer to output vector
*              - const polyveck *v0: pointer to low part of input vector
*              - const polyveck *v1: pointer to high part of input vector
*
* Returns number of 1 bits.
**************************************************/
/*
 * This function creates a "hint" during signing.
 * The hint tells the verifier whether a "carry" happened during the
 * `w = y + c*s2` calculation.
 * The verifier uses this hint to correctly reconstruct the high bits of `w`
 * without knowing the secret `s2`.
 * It returns the number of "1s" in the hint, which must be less than OMEGA.
 */
unsigned int polyveck_make_hint(polyveck *h,
                                const polyveck *v0,
                                const polyveck *v1)
{
  unsigned int i, s = 0;

  for(i = 0; i < K; ++i)
    s += poly_make_hint(&h->vec[i], &v0->vec[i], &v1->vec[i]);

  return s;
}

/*************************************************
* Name:        polyveck_use_hint
*
* Description: Use hint vector to correct the high bits of input vector.
*
* Arguments:   - polyveck *w: pointer to output vector of polynomials with
*                             corrected high bits
*              - const polyveck *u: pointer to input vector
*              - const polyveck *h: pointer to input hint vector
**************************************************/
/*
 * During verification, this function uses the `hint` from the signature
 * to correct the calculated high bits of `w`.
 * This allows the verifier to check if `w` was constructed correctly,
 * proving the signer knew the secret key.
 */
void polyveck_use_hint(polyveck *w, const polyveck *u, const polyveck *h) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_use_hint(&w->vec[i], &u->vec[i], &h->vec[i]);
}

/*
 * This function "packs" a polynomial vector `w1` into a byte array `r`.
 * Packing converts the list of numbers (coefficients) into a compact
 * sequence of bytes that can be easily stored or sent over the internet.
 * This is used to form part of the signature.
 */
void polyveck_pack_w1(uint8_t r[K*POLYW1_PACKEDBYTES], const polyveck *w1) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    polyw1_pack(&r[i*POLYW1_PACKEDBYTES], &w1->vec[i]);
}
