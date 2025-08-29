/*
 * =================================================================================
 * This file is a super important part of our crypto project.
 * It implements the Keccak algorithm, which is like a super-secure blender for data.
 * Keccak is the foundation for the modern SHA-3 hash functions and SHAKE extendable-output
 * functions, as defined in the official U.S. government standard FIPS 202.
 *
 * HASH FUNCTION in general:
 * A hash function takes any input (like a message, a file, or a number) and "hashes"
 * it into a short, fixed-size string of characters. This output is called a "hash"
 * or "digest". It's like a unique fingerprint for the data.
 *    1. The output is always the same length.
 *    2. The same input always produces the same output.
 *    3. A tiny change in the input creates a completely different output.
 *    4. You can't go backward from the output to figure out the input (it's a one-way street).
 *
 * WHAT ARE SHAKE and SHA-3?
 * - SHA-3: This is a standard hash function. `sha3_256` gives you a 256-bit (32-byte)
 *   fingerprint. It's used for creating a unique and secure summary of data.
 * - SHAKE: This stands for "Secure Hash Algorithm and Keccak". It's an "eXtendable-Output
 *   Function" (XOF). This means that unlike SHA-3, you can ask it for an output of
 *   *any length*. If you need a lot of random-looking data generated from a single,
 *   small seed, SHAKE is the perfect tool. In Dilithium, we use SHAKE all the time
 *   to expand seeds into large matrices and vectors of polynomials.
 *
 * HOW DOES IT WORK (THE SPONGE METAPHOR)?
 * Keccak works like a "sponge".
 * 1. ABSORBING: You feed your input data into the sponge. The sponge "absorbs" the
 *    data, mixing it into its internal state.
 * 2. SQUEEZING: After you've put all your data in, you "squeeze" the sponge to get
 *    the output hash. For SHAKE, you can keep squeezing to get as much output as you need.
 *
 * This file contains all the low-level functions to make this sponge work.
 *
 * Based on the public domain implementation in crypto_hash/keccakc512/simple/
 * from http://bench.cr.yp.to/supercop.html by Ronny Van Keer and the public domain
 * "TweetFips202" implementation from https://twitter.com/tweetfips202 by
 * Gilles Van Assche, Daniel J. Bernstein, and Peter Schwabe.
 * =================================================================================
 */

#include <stddef.h>
#include <stdint.h>
#include "fips202.h"

#define NROUNDS 24
#define ROL(a, offset) ((a << offset) ^ (a >> (64-offset)))
// The below pair of functions designed to convert between a 64-bit integer and its byte array
// representation. This pair of functions ensures that the byte order is always correctly
// arranged in little-endian format, making data exchange between different systems consistent.
/*************************************************
* Name:        load64
*
* Description: Load 8 bytes into a 64-bit integer (uint64_t).
*              This function handles "endianness" by making sure the bytes
*              are arranged in the correct order (little-endian).
*
* Arguments:   - const uint8_t *x: pointer to input byte array
*
* Returns the loaded 64-bit unsigned integer
**************************************************/
static uint64_t load64(const uint8_t x[8]) {
  unsigned int i;
  uint64_t r = 0;

  for(i=0;i<8;i++)
    r |= (uint64_t)x[i] << 8*i;

  return r;
}

/*************************************************
* Name:        store64
*
* Description: Store a 64-bit integer into an array of 8 bytes (little-endian).
*              This is the reverse of `load64`.
*
* Arguments:   - uint8_t *x: pointer to the output byte array
*              - uint64_t u: input 64-bit unsigned integer
**************************************************/
static void store64(uint8_t x[8], uint64_t u) {
  unsigned int i;

  for(i=0;i<8;i++)
    x[i] = u >> 8*i;
}

/* These are the "round constants" for Keccak. They are special, predefined
 * numbers that are mixed into the state during each round of permutation.
 * They are essential for the security of the algorithm because they break
 * up symmetries that might otherwise create weaknesses. */
const uint64_t KeccakF_RoundConstants[NROUNDS] = {
  (uint64_t)0x0000000000000001ULL,
  (uint64_t)0x0000000000008082ULL,
  (uint64_t)0x800000000000808aULL,
  (uint64_t)0x8000000080008000ULL,
  (uint64_t)0x000000000000808bULL,
  (uint64_t)0x0000000080000001ULL,
  (uint64_t)0x8000000080008081ULL,
  (uint64_t)0x8000000000008009ULL,
  (uint64_t)0x000000000000008aULL,
  (uint64_t)0x0000000000000088ULL,
  (uint64_t)0x0000000080008009ULL,
  (uint64_t)0x000000008000000aULL,
  (uint64_t)0x000000008000808bULL,
  (uint64_t)0x800000000000008bULL,
  (uint64_t)0x8000000000008089ULL,
  (uint64_t)0x8000000000008003ULL,
  (uint64_t)0x8000000000008002ULL,
  (uint64_t)0x8000000000000080ULL,
  (uint64_t)0x000000000000800aULL,
  (uint64_t)0x800000008000000aULL,
  (uint64_t)0x8000000080008081ULL,
  (uint64_t)0x8000000000008080ULL,
  (uint64_t)0x0000000080000001ULL,
  (uint64_t)0x8000000080008008ULL
};

/*************************************************
* Name:        KeccakF1600_StatePermute
*
* Description: This is the heart of the Keccak algorithm. It's the "blender"
*              function that scrambles the internal state. It takes the 1600-bit
*              state (25x 64-bit words) and mixes it up thoroughly over 24 rounds.
*              Each round involves a series of bitwise operations (XOR, AND, NOT,
*              and rotations) that are designed to be fast on computers and provide
*              strong cryptographic security.
*
*              The 5 steps of each round are:
*              - θ (theta): A linear mixing step that provides diffusion.
*              - ρ (rho):   Bitwise rotation to further spread changes.
*              - π (pi):    Permutes the positions of the words.
*              - χ (chi):   The main non-linear step, providing confusion.
*              - ι (iota):  Adds a round-specific constant to break symmetry.
*               (It's a bitwise XOR operation that breaks symmetry)
*               (--> prevents the algorithm from having fixed points or repetitive
*                    cycles that could be exploited by attackers --> ensures every round
*                    is distinct, making the process much harder to reverse)
*
*              This implementation is highly optimized: it unrolls the loop to
*              process two rounds at a time and uses local variables to encourage
*              the compiler to use CPU registers.
*
* Arguments:   - uint64_t *state: pointer to input/output Keccak state
**************************************************/
static void KeccakF1600_StatePermute(uint64_t state[25])
{
        int round;

        // These local variables represent the 5x5 state matrix of Keccak.
        // The naming convention corresponds to the matrix coordinates, e.g.,
        // Aba = A[0,0], Abe = A[1,0], Aga = A[0,1], etc.
        // Using local variables helps the compiler optimize by keeping them in registers.
        uint64_t Aba, Abe, Abi, Abo, Abu;
        uint64_t Aga, Age, Agi, Ago, Agu;
        uint64_t Aka, Ake, Aki, Ako, Aku;
        uint64_t Ama, Ame, Ami, Amo, Amu;
        uint64_t Asa, Ase, Asi, Aso, Asu;
        // Intermediate variables for calculations
        uint64_t BCa, BCe, BCi, BCo, BCu;
        uint64_t Da, De, Di, Do, Du;
        // Temporary state variables for the second round in the unrolled loop
        uint64_t Eba, Ebe, Ebi, Ebo, Ebu;
        uint64_t Ega, Ege, Egi, Ego, Egu;
        uint64_t Eka, Eke, Eki, Eko, Eku;
        uint64_t Ema, Eme, Emi, Emo, Emu;
        uint64_t Esa, Ese, Esi, Eso, Esu;

        // Copy the input state array into local variables (the 'A' state matrix).
        Aba = state[ 0];
        Abe = state[ 1];
        Abi = state[ 2];
        Abo = state[ 3];
        Abu = state[ 4];
        Aga = state[ 5];
        Age = state[ 6];
        Agi = state[ 7];
        Ago = state[ 8];
        Agu = state[ 9];
        Aka = state[10];
        Ake = state[11];
        Aki = state[12];
        Ako = state[13];
        Aku = state[14];
        Ama = state[15];
        Ame = state[16];
        Ami = state[17];
        Amo = state[18];
        Amu = state[19];
        Asa = state[20];
        Ase = state[21];
        Asi = state[22];
        Aso = state[23];
        Asu = state[24];

        // The loop is unrolled to process two rounds at a time (round and round+1).
        // Unrolled loop = rewriting a loop as repeated statements to run faster
        for(round = 0; round < NROUNDS; round += 2) {
            // --- ROUND (round) ---

            // θ (Theta) step - Part 1: Compute column parities
            BCa = Aba^Aga^Aka^Ama^Asa;
            BCe = Abe^Age^Ake^Ame^Ase;
            BCi = Abi^Agi^Aki^Ami^Asi;
            BCo = Abo^Ago^Ako^Amo^Aso;
            BCu = Abu^Agu^Aku^Amu^Asu;

            // θ (Theta) step - Part 2: Compute mixing values (D) and apply to state A
            Da = BCu^ROL(BCe, 1);
            De = BCa^ROL(BCi, 1);
            Di = BCe^ROL(BCo, 1);
            Do = BCi^ROL(BCu, 1);
            Du = BCo^ROL(BCa, 1);

            Aba ^= Da;
            Age ^= De;
            Aki ^= Di;
            Amo ^= Do;
            Asu ^= Du;
            Abo ^= Do;
            Agu ^= Du;
            Aka ^= Da;
            Ame ^= De;
            Asi ^= Di;
            Abe ^= De;
            Agi ^= Di;
            Ako ^= Do;
            Amu ^= Du;
            Asa ^= Da;
            Abu ^= Du;
            Aga ^= Da;
            Ake ^= De;
            Ami ^= Di;
            Aso ^= Do;
            Abi ^= Di;
            Ago ^= Do;
            Aku ^= Du;
            Ama ^= Da;
            Ase ^= De;

            // ρ (Rho) and π (Pi) steps are combined and interleaved.
            // The results are then used in the χ (Chi) step.
            // The ROL (rotate left) is the ρ step.
            // The rearrangement of variables (e.g., from Age to BCe) is the π step.
            BCa = Aba;
            BCe = ROL(Age, 44);
            BCi = ROL(Aki, 43);
            BCo = ROL(Amo, 21);
            BCu = ROL(Asu, 14);
            // χ (Chi) step: The main non-linear operation. Result stored in E.
            // ι (Iota) step: XOR with round constant.
            Eba =   BCa ^((~BCe)&  BCi );
            Eba ^= (uint64_t)KeccakF_RoundConstants[round];
            Ebe =   BCe ^((~BCi)&  BCo );
            Ebi =   BCi ^((~BCo)&  BCu );
            Ebo =   BCo ^((~BCu)&  BCa );
            Ebu =   BCu ^((~BCa)&  BCe );

            BCa = ROL(Abo, 28);
            BCe = ROL(Agu, 20);
            BCi = ROL(Aka,  3);
            BCo = ROL(Ame, 45);
            BCu = ROL(Asi, 61);
            Ega =   BCa ^((~BCe)&  BCi );
            Ege =   BCe ^((~BCi)&  BCo );
            Egi =   BCi ^((~BCo)&  BCu );
            Ego =   BCo ^((~BCu)&  BCa );
            Egu =   BCu ^((~BCa)&  BCe );

            BCa = ROL(Abe,  1);
            BCe = ROL(Agi,  6);
            BCi = ROL(Ako, 25);
            BCo = ROL(Amu,  8);
            BCu = ROL(Asa, 18);
            Eka =   BCa ^((~BCe)&  BCi );
            Eke =   BCe ^((~BCi)&  BCo );
            Eki =   BCi ^((~BCo)&  BCu );
            Eko =   BCo ^((~BCu)&  BCa );
            Eku =   BCu ^((~BCa)&  BCe );

            BCa = ROL(Abu, 27);
            BCe = ROL(Aga, 36);
            BCi = ROL(Ake, 10);
            BCo = ROL(Ami, 15);
            BCu = ROL(Aso, 56);
            Ema =   BCa ^((~BCe)&  BCi );
            Eme =   BCe ^((~BCi)&  BCo );
            Emi =   BCi ^((~BCo)&  BCu );
            Emo =   BCo ^((~BCu)&  BCa );
            Emu =   BCu ^((~BCa)&  BCe );

            BCa = ROL(Abi, 62);
            BCe = ROL(Ago, 55);
            BCi = ROL(Aku, 39);
            BCo = ROL(Ama, 41);
            BCu = ROL(Ase,  2);
            Esa =   BCa ^((~BCe)&  BCi );
            Ese =   BCe ^((~BCi)&  BCo );
            Esi =   BCi ^((~BCo)&  BCu );
            Eso =   BCo ^((~BCu)&  BCa );
            Esu =   BCu ^((~BCa)&  BCe );

            // --- ROUND (round+1) ---
            // The same steps are repeated, but this time the input is the 'E' matrix
            // from the previous round, and the output is the 'A' matrix.

            // θ (Theta) step
            BCa = Eba^Ega^Eka^Ema^Esa;
            BCe = Ebe^Ege^Eke^Eme^Ese;
            BCi = Ebi^Egi^Eki^Emi^Esi;
            BCo = Ebo^Ego^Eko^Emo^Eso;
            BCu = Ebu^Egu^Eku^Emu^Esu;

            Da = BCu^ROL(BCe, 1);
            De = BCa^ROL(BCi, 1);
            Di = BCe^ROL(BCo, 1);
            Do = BCi^ROL(BCu, 1);
            Du = BCo^ROL(BCa, 1);

            Eba ^= Da;
            Ege ^= De;
            Eki ^= Di;
            Emo ^= Do;
            Esu ^= Du;
            Ebo ^= Do;
            Egu ^= Du;
            Eka ^= Da;
            Eme ^= De;
            Esi ^= Di;
            Ebe ^= De;
            Egi ^= Di;
            Eko ^= Do;
            Emu ^= Du;
            Esa ^= Da;
            Ebu ^= Du;
            Ega ^= Da;
            Eke ^= De;
            Emi ^= Di;
            Eso ^= Do;
            Ebi ^= Di;
            Ego ^= Do;
            Eku ^= Du;
            Ema ^= Da;
            Ese ^= De;

            // ρ (Rho), π (Pi), χ (Chi), and ι (Iota) steps
            BCa = Eba;
            BCe = ROL(Ege, 44);
            BCi = ROL(Eki, 43);
            BCo = ROL(Emo, 21);
            BCu = ROL(Esu, 14);
            Aba =   BCa ^((~BCe)&  BCi );
            Aba ^= (uint64_t)KeccakF_RoundConstants[round+1];
            Abe =   BCe ^((~BCi)&  BCo );
            Abi =   BCi ^((~BCo)&  BCu );
            Abo =   BCo ^((~BCu)&  BCa );
            Abu =   BCu ^((~BCa)&  BCe );

            BCa = ROL(Ebo, 28);
            BCe = ROL(Egu, 20);
            BCi = ROL(Eka, 3);
            BCo = ROL(Eme, 45);
            BCu = ROL(Esi, 61);
            Aga =   BCa ^((~BCe)&  BCi );
            Age =   BCe ^((~BCi)&  BCo );
            Agi =   BCi ^((~BCo)&  BCu );
            Ago =   BCo ^((~BCu)&  BCa );
            Agu =   BCu ^((~BCa)&  BCe );

            BCa = ROL(Ebe, 1);
            BCe = ROL(Egi, 6);
            BCi = ROL(Eko, 25);
            BCo = ROL(Emu, 8);
            BCu = ROL(Esa, 18);
            Aka =   BCa ^((~BCe)&  BCi );
            Ake =   BCe ^((~BCi)&  BCo );
            Aki =   BCi ^((~BCo)&  BCu );
            Ako =   BCo ^((~BCu)&  BCa );
            Aku =   BCu ^((~BCa)&  BCe );

            BCa = ROL(Ebu, 27);
            BCe = ROL(Ega, 36);
            BCi = ROL(Eke, 10);
            BCo = ROL(Emi, 15);
            BCu = ROL(Eso, 56);
            Ama =   BCa ^((~BCe)&  BCi );
            Ame =   BCe ^((~BCi)&  BCo );
            Ami =   BCi ^((~BCo)&  BCu );
            Amo =   BCo ^((~BCu)&  BCa );
            Amu =   BCu ^((~BCa)&  BCe );

            BCa = ROL(Ebi, 62);
            BCe = ROL(Ego, 55);
            BCi = ROL(Eku, 39);
            BCo = ROL(Ema, 41);
            BCu = ROL(Ese, 2);
            Asa =   BCa ^((~BCe)&  BCi );
            Ase =   BCe ^((~BCi)&  BCo );
            Asi =   BCi ^((~BCo)&  BCu );
            Aso =   BCo ^((~BCu)&  BCa );
            Asu =   BCu ^((~BCa)&  BCe );
        }

        // Copy the final local state variables back to the state array.
        state[ 0] = Aba;
        state[ 1] = Abe;
        state[ 2] = Abi;
        state[ 3] = Abo;
        state[ 4] = Abu;
        state[ 5] = Aga;
        state[ 6] = Age;
        state[ 7] = Agi;
        state[ 8] = Ago;
        state[ 9] = Agu;
        state[10] = Aka;
        state[11] = Ake;
        state[12] = Aki;
        state[13] = Ako;
        state[14] = Aku;
        state[15] = Ama;
        state[16] = Ame;
        state[17] = Ami;
        state[18] = Amo;
        state[19] = Amu;
        state[20] = Asa;
        state[21] = Ase;
        state[22] = Asi;
        state[23] = Aso;
        state[24] = Asu;
}

/*************************************************
* Name:        keccak_init
*
* Description: Initializes the Keccak state by setting it all to zeros.
*              This is the first step before you start absorbing data.
*
* Arguments:   - uint64_t *s: pointer to Keccak state
**************************************************/
static void keccak_init(uint64_t s[25])
{
  unsigned int i;
  for(i=0;i<25;i++)
    s[i] = 0;
}

/*************************************************
* Name:        keccak_absorb
*
* Description: The "absorbing" part of the sponge. It takes input data and
*              XORs it into the Keccak state. If the state gets full, it
*              calls the `KeccakF1600_StatePermute` function to scramble
*              everything before absorbing more data.
*
* Arguments:   - uint64_t *s: pointer to Keccak state
*              - unsigned int r: the "rate" in bytes (how much data to absorb per block)
*              - const uint8_t *in: pointer to input data
*              - size_t inlen: length of input in bytes
*
* Returns new position in the current block.
**************************************************/
static unsigned int keccak_absorb(uint64_t s[25],
                                  unsigned int pos,
                                  unsigned int r,
                                  const uint8_t *in,
                                  size_t inlen)
{
  unsigned int i;

  while(pos+inlen >= r) {
    for(i=pos;i<r;i++)
      s[i/8] ^= (uint64_t)*in++ << 8*(i%8);
    inlen -= r-pos;
    KeccakF1600_StatePermute(s);
    pos = 0;
  }

  for(i=pos;i<pos+inlen;i++)
    s[i/8] ^= (uint64_t)*in++ << 8*(i%8);

  return i;
}

/*************************************************
* Name:        keccak_finalize
*
* Description: Finishes the absorbing phase. It adds special "padding" to the
*              last block of data. This padding is crucial for security and
*              ensures that messages of different lengths produce different hashes.
*              The byte `p` is for "domain separation", which means we use
*              different values of `p` for SHAKE vs. SHA-3 to make sure their
*              outputs are independent.
*
* Arguments:   - uint64_t *s: pointer to Keccak state
*              - unsigned int pos: current position in the block
*              - unsigned int r: the rate in bytes
*              - uint8_t p: the domain separation byte
**************************************************/
static void keccak_finalize(uint64_t s[25], unsigned int pos, unsigned int r, uint8_t p)
{
  s[pos/8] ^= (uint64_t)p << 8*(pos%8);
  s[r/8-1] ^= 1ULL << 63;
}

/*************************************************
* Name:        keccak_squeeze
*
* Description: The "squeezing" part of the sponge. It extracts output data
*              from the state. If you need more output than is currently in
*              the state, it calls `KeccakF1600_StatePermute` to generate a
*              new block of output. This can be called multiple times to get
*              an output of any length (which is what SHAKE does).
*
* Arguments:   - uint8_t *out: pointer to output buffer
*              - size_t outlen: number of bytes to squeeze
*              - uint64_t *s: pointer to Keccak state
*              - unsigned int pos: position in the current block
*              - unsigned int r: the rate in bytes
*
* Returns new position in the current block.
**************************************************/
static unsigned int keccak_squeeze(uint8_t *out,
                                   size_t outlen,
                                   uint64_t s[25],
                                   unsigned int pos,
                                   unsigned int r)
{
  unsigned int i;

  while(outlen) {
    if(pos == r) {
      KeccakF1600_StatePermute(s);
      pos = 0;
    }
    for(i=pos;i < r && i < pos+outlen; i++)
      *out++ = s[i/8] >> 8*(i%8);
    outlen -= i-pos;
    pos = i;
  }

  return pos;
}


/*************************************************
* Name:        keccak_absorb_once
*
* Description: A helper function that does the whole absorb process in one go
*              (initializes, absorbs all data, and finalizes). This is faster
*              for when you have all the input data at once.
*
* Arguments:   - uint64_t *s: pointer to Keccak state
*              - unsigned int r: the rate in bytes
*              - const uint8_t *in: pointer to input data
*              - size_t inlen: length of input in bytes
*              - uint8_t p: domain separation byte
**************************************************/
static void keccak_absorb_once(uint64_t s[25],
                               unsigned int r,
                               const uint8_t *in,
                               size_t inlen,
                               uint8_t p)
{
  unsigned int i;

  for(i=0;i<25;i++)
    s[i] = 0;

  while(inlen >= r) {
    for(i=0;i<r/8;i++)
      s[i] ^= load64(in+8*i);
    in += r;
    inlen -= r;
    KeccakF1600_StatePermute(s);
  }

  for(i=0;i<inlen;i++)
    s[i/8] ^= (uint64_t)in[i] << 8*(i%8);

  s[i/8] ^= (uint64_t)p << 8*(i%8);
  s[(r-1)/8] ^= 1ULL << 63;
}

/*************************************************
* Name:        keccak_squeezeblocks
*
* Description: A faster version of `keccak_squeeze` that is optimized for when
*              you need to squeeze many full blocks of output at once.
*
* Arguments:   - uint8_t *out: pointer to output buffer
*              - size_t nblocks: number of full blocks to squeeze
*              - uint64_t *s: pointer to Keccak state
*              - unsigned int r: the rate in bytes
**************************************************/
static void keccak_squeezeblocks(uint8_t *out,
                                 size_t nblocks,
                                 uint64_t s[25],
                                 unsigned int r)
{
  unsigned int i;

  while(nblocks) {
    KeccakF1600_StatePermute(s);
    for(i=0;i<r/8;i++)
      store64(out+8*i, s[i]);
    out += r;
    nblocks -= 1;
  }
}

/*
 * =================================================================================
 * SHAKE128 Functions
 *
 * These functions provide an API for using SHAKE128. SHAKE128 has a security
 * level of 128 bits and can produce an output of any length.
 * The `_RATE` is 168 bytes (1344 bits).
 * =================================================================================
 */

/*************************************************
* Name:        shake128_init
*
* Description: Initializes a SHAKE128 state.
*
* Arguments:   - keccak_state *state: pointer to Keccak state
**************************************************/
void shake128_init(keccak_state *state)
{
  keccak_init(state->s);
  state->pos = 0;
}

/*************************************************
* Name:        shake128_absorb
*
* Description: Absorbs input data into the SHAKE128 state. You can call this
*              multiple times if your data comes in chunks.
*
* Arguments:   - keccak_state *state: pointer to Keccak state
*              - const uint8_t *in: pointer to input data
*              - size_t inlen: length of input in bytes
**************************************************/
void shake128_absorb(keccak_state *state, const uint8_t *in, size_t inlen)
{
  state->pos = keccak_absorb(state->s, state->pos, SHAKE128_RATE, in, inlen);
}

/*************************************************
* Name:        shake128_finalize
*
* Description: Finalizes the absorbing phase for SHAKE128.
*
* Arguments:   - keccak_state *state: pointer to Keccak state
**************************************************/
void shake128_finalize(keccak_state *state)
{
  keccak_finalize(state->s, state->pos, SHAKE128_RATE, 0x1F);
  state->pos = SHAKE128_RATE;
}

/*************************************************
* Name:        shake128_squeeze
*
* Description: Squeezes output from the SHAKE128 state. Can be called
*              multiple times to get more output.
*
* Arguments:   - uint8_t *out: pointer to output buffer
*              - size_t outlen : number of bytes to squeeze
*              - keccak_state *state: pointer to Keccak state
**************************************************/
void shake128_squeeze(uint8_t *out, size_t outlen, keccak_state *state)
{
  state->pos = keccak_squeeze(out, outlen, state->s, state->pos, SHAKE128_RATE);
}

/*************************************************
* Name:        shake128_absorb_once
*
* Description: A helper that combines init, absorb, and finalize for SHAKE128.
*
* Arguments:   - keccak_state *state: pointer to Keccak state
*              - const uint8_t *in: pointer to input data
*              - size_t inlen: length of input in bytes
**************************************************/
void shake128_absorb_once(keccak_state *state, const uint8_t *in, size_t inlen)
{
  keccak_absorb_once(state->s, SHAKE128_RATE, in, inlen, 0x1F);
  state->pos = SHAKE128_RATE;
}

/*************************************************
* Name:        shake128_squeezeblocks
*
* Description: Squeezes full blocks of output from a SHAKE128 state.
*
* Arguments:   - uint8_t *out: pointer to output buffer
*              - size_t nblocks: number of blocks to squeeze
*              - keccak_state *state: pointer to Keccak state
**************************************************/
void shake128_squeezeblocks(uint8_t *out, size_t nblocks, keccak_state *state)
{
  keccak_squeezeblocks(out, nblocks, state->s, SHAKE128_RATE);
}

/*
 * =================================================================================
 * SHAKE256 Functions
 *
 * These functions provide an API for using SHAKE256. SHAKE256 has a security
 * level of 256 bits and can produce an output of any length.
 * The `_RATE` is 136 bytes (1088 bits).
 * =================================================================================
 */

/*************************************************
* Name:        shake256_init
*
* Description: Initializes a SHAKE256 state.
*
* Arguments:   - keccak_state *state: pointer to Keccak state
**************************************************/
void shake256_init(keccak_state *state)
{
  keccak_init(state->s);
  state->pos = 0;
}

/*************************************************
* Name:        shake256_absorb
*
* Description: Absorbs input data into the SHAKE256 state.
*
* Arguments:   - keccak_state *state: pointer to Keccak state
*              - const uint8_t *in: pointer to input data
*              - size_t inlen: length of input in bytes
**************************************************/
void shake256_absorb(keccak_state *state, const uint8_t *in, size_t inlen)
{
  state->pos = keccak_absorb(state->s, state->pos, SHAKE256_RATE, in, inlen);
}

/*************************************************
* Name:        shake256_finalize
*
* Description: Finalizes the absorbing phase for SHAKE256.
*
* Arguments:   - keccak_state *state: pointer to Keccak state
**************************************************/
void shake256_finalize(keccak_state *state)
{
  keccak_finalize(state->s, state->pos, SHAKE256_RATE, 0x1F);
  state->pos = SHAKE256_RATE;
}

/*************************************************
* Name:        shake256_squeeze
*
* Description: Squeezes output from the SHAKE256 state.
*
* Arguments:   - uint8_t *out: pointer to output buffer
*              - size_t outlen : number of bytes to squeeze
*              - keccak_state *state: pointer to Keccak state
**************************************************/
void shake256_squeeze(uint8_t *out, size_t outlen, keccak_state *state)
{
  state->pos = keccak_squeeze(out, outlen, state->s, state->pos, SHAKE256_RATE);
}

/*************************************************
* Name:        shake256_absorb_once
*
* Description: A helper that combines init, absorb, and finalize for SHAKE256.
*
* Arguments:   - keccak_state *state: pointer to Keccak state
*              - const uint8_t *in: pointer to input data
*              - size_t inlen: length of input in bytes
**************************************************/
void shake256_absorb_once(keccak_state *state, const uint8_t *in, size_t inlen)
{
  keccak_absorb_once(state->s, SHAKE256_RATE, in, inlen, 0x1F);
  state->pos = SHAKE256_RATE;
}

/*************************************************
* Name:        shake256_squeezeblocks
*
* Description: Squeezes full blocks of output from a SHAKE256 state.
*
* Arguments:   - uint8_t *out: pointer to output buffer
*              - size_t nblocks: number of blocks to squeeze
*              - keccak_state *state: pointer to Keccak state
**************************************************/
void shake256_squeezeblocks(uint8_t *out, size_t nblocks, keccak_state *state)
{
  keccak_squeezeblocks(out, nblocks, state->s, SHAKE256_RATE);
}

/*
 * =================================================================================
 * Simple "One-Shot" API Functions
 *
 * These functions are for when you have all your input data at once and just
 * want to get the hash in a single function call.
 * =================================================================================
 */

/*************************************************
* Name:        shake128
*
* Description: A simple, one-shot function to compute a SHAKE128 hash.
*
* Arguments:   - uint8_t *out: pointer to output buffer
*              - size_t outlen: requested output length in bytes
*              - const uint8_t *in: pointer to input data
*              - size_t inlen: length of input in bytes
**************************************************/
void shake128(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen)
{
  size_t nblocks;
  keccak_state state;

  shake128_absorb_once(&state, in, inlen);
  nblocks = outlen/SHAKE128_RATE;
  shake128_squeezeblocks(out, nblocks, &state);
  outlen -= nblocks*SHAKE128_RATE;
  out += nblocks*SHAKE128_RATE;
  shake128_squeeze(out, outlen, &state);
}

/*************************************************
* Name:        shake256
*
* Description: A simple, one-shot function to compute a SHAKE256 hash.
*
* Arguments:   - uint8_t *out: pointer to output buffer
*              - size_t outlen: requested output length in bytes
*              - const uint8_t *in: pointer to input data
*              - size_t inlen: length of input in bytes
**************************************************/
void shake256(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen)
{
  size_t nblocks;
  keccak_state state;

  shake256_absorb_once(&state, in, inlen);
  nblocks = outlen/SHAKE256_RATE;
  shake256_squeezeblocks(out, nblocks, &state);
  outlen -= nblocks*SHAKE256_RATE;
  out += nblocks*SHAKE256_RATE;
  shake256_squeeze(out, outlen, &state);
}

/*************************************************
* Name:        sha3_256
*
* Description: A simple, one-shot function to compute a SHA3-256 hash.
*              This produces a fixed-size 32-byte output.
*
* Arguments:   - uint8_t *h: pointer to output (must be 32 bytes)
*              - const uint8_t *in: pointer to input data
*              - size_t inlen: length of input in bytes
**************************************************/
void sha3_256(uint8_t h[32], const uint8_t *in, size_t inlen)
{
  unsigned int i;
  uint64_t s[25];

  keccak_absorb_once(s, SHA3_256_RATE, in, inlen, 0x06);
  KeccakF1600_StatePermute(s);
  for(i=0;i<4;i++)
    store64(h+8*i,s[i]);
}

/*************************************************
* Name:        sha3_512
*
* Description: A simple, one-shot function to compute a SHA3-512 hash.
*              This produces a fixed-size 64-byte output.
*
* Arguments:   - uint8_t *h: pointer to output (must be 64 bytes)
*              - const uint8_t *in: pointer to input data
*              - size_t inlen: length of input in bytes
**************************************************/
void sha3_512(uint8_t h[64], const uint8_t *in, size_t inlen)
{
  unsigned int i;
  uint64_t s[25];

  keccak_absorb_once(s, SHA3_512_RATE, in, inlen, 0x06);
  KeccakF1600_StatePermute(s);
  for(i=0;i<8;i++)
    store64(h+8*i,s[i]);
}
