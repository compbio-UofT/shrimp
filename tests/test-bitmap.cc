#include <iostream>

#include "../common/bitmap.h"

using namespace std;

int main() {
  char c;
  int bitsPerField, index, value;
  BITMAP_TYPE a;

  while (cin >> c, !cin.eof()) {
    switch (c) {
    case 'i':
      cin >> bitsPerField >> index >> value;
      bitmapInsertClean(&a, bitsPerField, index, value);
      break;
    case 'c':
      cin >> bitsPerField >> index;
      bitmapClear(&a, bitsPerField, index);
      break;
    case 'q':
      cin >> bitsPerField >> index;
      cout << bitmapExtract(&a, bitsPerField, index) << endl;
      break;
    }
  }

  return 0;
}
