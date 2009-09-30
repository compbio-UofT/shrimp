#ifndef __BITMAP_H
#define __BITMAP_H

#include <stdio.h>
#include <assert.h>

#ifdef NDEBUG
#define INLINE_ND inline
#else
#define INLINE_ND
#endif


#define BITMAP_TYPE uint64_t


#define BITMAP_ALL1_FIELD(_bitsPerField) \
  ((((BITMAP_TYPE)1) << (_bitsPerField)) - 1)

#define BITMAP32_ALL1_FIELD(_bitsPerField) \
  ((((uint32_t)1) << (_bitsPerField)) - 1)



INLINE_ND void
bitmapInsertClean(BITMAP_TYPE *a, uint bitsPerField, uint index, uint value) {
  assert(index < 64/bitsPerField);
  assert((value & ~BITMAP_ALL1_FIELD(bitsPerField)) == 0);

  (*a) |= (((BITMAP_TYPE)value) << index*bitsPerField);
}


INLINE_ND void
bitmapClear(BITMAP_TYPE *a, uint bitsPerField, uint index) {
  assert(index < 64/bitsPerField);

  (*a) &= ~(BITMAP_ALL1_FIELD(bitsPerField) << index*bitsPerField);
}


INLINE_ND void
bitmapPrepend(BITMAP_TYPE *a, uint bitsPerField, uint value) {
  assert((value & ~BITMAP_ALL1_FIELD(bitsPerField)) == 0);

  (*a) <<= bitsPerField;
  (*a) |= value;
}



INLINE_ND uint
bitmapExtract(BITMAP_TYPE *a, uint bitsPerField, uint index) {
  assert(index < 64/bitsPerField);

  return ((*a) >> bitsPerField*index) & BITMAP_ALL1_FIELD(bitsPerField);
}


INLINE_ND uint
bitmap32Extract(uint32_t *a, uint bitsPerField, uint index) {
  assert(index < 32/bitsPerField);

  return ((*a) >> bitsPerField*index) & BITMAP32_ALL1_FIELD(bitsPerField);
}


char *
bitmap32VString(uint32_t *a, uint nFields, uint bitsPerField,
		bool reverse, bool useBits, bool useCommas) {
  static char buffer[1001];

  int i, j, k;
  uint32_t val;
  uint fieldsPerByte = 32/bitsPerField;

  int first, last, increment;

  if (!reverse) {
    first = nFields - 1;
    last = 0;
    increment = -1;
  } else {
    first = 0;
    last = nFields - 1;
    increment = 1;
  }

  for (i = first, j = 0;
       i != last + increment && j < 990;
       i += increment) {
    val = ((a[i/fieldsPerByte] >> bitsPerField*(i%fieldsPerByte))
	   & ((((uint32_t)1) << bitsPerField) - 1));
    if (useBits) {
      for (k = bitsPerField - 1; k >= 0; k--) {
	buffer[j++] = (bitmap32Extract(&val, 1, k) == 1? '1' : '0');
      }
    } else {
      j += sprintf(&buffer[j], "%u", val); // unsafe
    }
    if (useCommas && i != last)
      buffer[j++] = ',';
  }
  buffer[j] = 0;

  return buffer;
}
      

#endif
