#include <iostream>
#include <cassert>

using namespace std;

#include "../common/anchors.h"

struct anchor a[100];
struct anchor dest;
char M[100][100];
uint xlen, ylen;

void set_m(int x, int y, char c) {
  if (0 <= x && x <= xlen-1 && 0 <= y && y <= ylen-1)
    M[y][x] = c;
}

int main() {
  uint n;

  memset(M, ' ', 100*100);

  cin >> xlen >> ylen;
  assert(xlen < 100 && ylen < 100);
  for (uint y = 0; y < ylen; y++) {
    M[y][xlen] = 0;
  }

  cin >> n;
  assert(n > 0 && n < 100);

  for (uint i = 0; i < n; i++) {
    int x, y, l, w, more, alt;

    cout << "anchor " << i << ": ";
    cin >> x >> y >> l >> w;
    a[i].x = x;
    a[i].y = y;
    a[i].length = l;
    a[i].width = w;
    a[i].more_than_once = 0;

    for (uint k = 0; k <= l-1; k++) {
      set_m(x+k, y+k, 'o');
      set_m(x+(w-1)+k, y-(w-1)+k, 'o');
    }
    for (uint k = 1; k < w-1; k++) {
      set_m(x+k, y-k, 'o');
      set_m(x+(l-1)+k, y+(l-1)-k, 'o');
    }
  }

  join_anchors(a, n, &dest);
  widen_anchor(&dest, 2);

  for (uint k = 0; k <= dest.length-1; k++) {
    set_m(dest.x+k, dest.y+k, '\\');
    set_m(dest.x+(dest.width-1)+k, dest.y-(dest.width-1)+k, '\\');
  }

  for (uint y = 0; y < ylen; y++) {
    cout << M[y] << endl;
  }
  cout << "---" << endl;

  cerr << "anchor x:" << dest.x << " y:" << dest.y
       << " l:" << (int)dest.length << " w:" << (int)dest.width << endl;

  for (uint y = 0; y < ylen; y++) {
    int x_min, x_max;

    get_x_range(&dest, xlen, ylen, y, &x_min, &x_max);
    cerr << "y:" << y << " x_min:" << x_min << " x_max:" << x_max << endl;
    for (uint x = 0; x < x_min; x++) {
      set_m(x, y, 'X');
    }
    for (uint x = x_max + 1; x < xlen; x++) {
      set_m(x, y, 'X');
    }
  }

  for (uint y = 0; y < ylen; y++) {
    cout << M[y] << endl;
  }

  return 0;
}
