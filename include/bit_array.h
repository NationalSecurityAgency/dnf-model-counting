// Paul Burkhardt, Ph.D.
// March 07, 2019
// March 20, 2019: added atomic write
// March 22, 2019: simplified return value for atomic write
//
// Bit array for compacting binary values.
//
// The smallest data type is one gmt::byte, e.g. char, which wastes seven bits when
// only one bit is need to denote a binary decision or value.
//
// This bit array can compactly store a single bit for each integer into an
// array of integers. The intuition is as follows.
//
// Let A be a gmt::byte array of n integers. The first 8-bit element of A can hold a
// bit for the first eight integers, and likewise each subsquent element holds a
// bit for the next multiple of eight integers. But since a gmt::byte holds only
// eight positions, then integers greater than eight will not fit. Thus all
// integers need to be transformed or mapped into a number space of eight.
//
// Recall that dividing rational numbers have remainders that can only range by
// the radix of the divisor. This is because division subtracts a multiple of
// the divisor that is closest to the number being divided and therefore the
// difference isn't larger than the divisor. Let b be our divisor, then dividing
// a number i at the last stage of the division by a multiple p of b will then
// have a remainder (i-p*b). This remainder is the modulo of the divisor.
//
// Mapping an integer i into a space of eight is simply taking the mod of i with
// eight. The index of the gmt::byte array A where i should be placed is the integer
// multiple of eight into i. If counting by one then a number less than eight
// will fit in the first element of A. If the number is more than two times
// eight but less than three times eight then it will go into the third element
// of A.
//
// Now let A and each 8-bit element in A be zero-indexed and unsigned. Also let
// each element of A be initialized to zero. In simple arithmetic the bit array
// operations are as follows:
//
// Given a number i, the index of A where i should be stored is then p=i/8.
//
// The bit position of i in A[i/8] is given by (i mod 8).
//
// The bit for number i in A[i/8] can be flipped by left shifting 1 by (i mod 8)
// and taking the bitwise OR with A[i/8]. But (i mod 8) can be replaced by
// (i-p*8).
//
// Thus the sequence of operations is then,
//
// A[p=i/8] = A[p=i/8] OR 1<<(i-p*8)
//
// Observe that 8 is a power of two. Hence expensive integer division and
// multiplication can be replaced by the more efficient method of shifting by
// three to the right for division and three to the left for multiplication. The
// leads to,
//
// p=i>>3
// A[p] = A[p] OR 1<<(i-(p<<3))
//
// Testing if a bit was flipped for a number i is simply taking the bitwise AND
// with the same operations for setting the bit, as follows.
//
// Is (A[p] & 1<<(i-(p<<3)) true?
//
// If the test is true then (1<<(i-(p<<3)) is returned, otherwise it is zero.
//------------------------------------------------------------------------------
#ifndef __BIT_ARRAY_H
#define __BIT_ARRAY_H
// using namespace std;
namespace gmt{typedef unsigned char byte;}
// using namespace gmt;
// typedef unsigned char gmt::byte;
//------------------------------------------------------------------------------
template<class B, class T>
class bit_array_write
{
public:
   bit_array_write() {}
   ~bit_array_write() {}
   void operator()(gmt::byte* a, T i)
   {
      int p = i>>3;
      a[p] |= 1<<(i-(p<<3));
   }
   void clear(gmt::byte* a, T i)
   {
      int p = i >> 3;
      a[p] &= ~(1<<(i-(p<<3)));
   }
    void clear_byte(gmt::byte* a, T i)
    {
        int p = i >> 3;
        a[p] = 0;
    }
};
template<class B, class T>
class bit_array_read
{
public:
   bit_array_read() {}
   ~bit_array_read() {}
   bool operator()(gmt::byte* a, T i)
   {
      int p = i>>3;
      return (a[p] & 1<<(i-(p<<3)));
   }
};
template<class B, class T>
class bit_array_write_atomic
{
public:
   bit_array_write_atomic() {}
   ~bit_array_write_atomic() {}
   bool operator()(gmt::byte* a, T i)
   {
      int p = i>>3;
      gmt::byte x = 1<<(i-(p<<3));
      gmt::byte y =__atomic_fetch_or(&a[p], x, __ATOMIC_RELAXED);
      return (x & y);
   }
};
#endif
