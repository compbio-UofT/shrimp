#include <iostream>
#include <cassert>

using namespace std;

#include "../common/anchors.h"

struct uw_anchor a[100];
struct uw_anchor dest;
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
    int x, y, l, weight;

    cout << "anchor " << i << ": ";
    cin >> x >> y >> l >> weight;
    a[i].x = x;
    a[i].y = y;
    a[i].length = l;
    a[i].n_kmers = weight;

    for (uint k = 0; k <= l-1; k++) {
      set_m(x+k, y+k, 'o');
    }
  }

  cout << "---" << endl;
  for (uint y = 0; y < ylen; y++) {
    cout << M[y] << endl;
  }
  cout << "---" << endl;

  for (uint i = 0; i < n; i++) {
    for (uint j = i+1; j < n; j++) {
      if (uw_anchors_colinear(&a[i], &a[j]))
	cout << "anchors " << i << " and " << j << " are colinear" << endl;
    }
  }

  for (uint i = 0; i < n; i++) {
    for (uint j = i+1; j < n; j++) {
      if (uw_anchors_intersect(&a[i], &a[j]))
	cout << "anchors " << i << " and " << j << " intersect" << endl;
    }
  }

  // bubble sort/merge
  bool none_found;
  uint m = n;
  do {
    none_found = true;
    for (uint i = 0; none_found && i < m; i++) {
      for (uint j = i + 1; none_found && j < m; j++) {
	if (uw_anchors_intersect(&a[i], &a[j])) {
	  none_found = false;
	  uw_anchors_join(&a[i], &a[j]);
	  a[j] = a[m-1];
	  m--;
	}
      }
    }
  } while (!none_found);

  memset(M, ' ', 100*100);
  for (uint y = 0; y < ylen; y++) {
    M[y][xlen] = 0;
  }
  cout << "---" << endl;
  for (uint i = 0; i < m; i++) {
    cout << "anchor i:" << i << " x:" << (uint)a[i].x << " y:" << (uint)a[i].y
	 << " l:" << (uint)a[i].length << " weight:" << (uint)a[i].n_kmers << endl;
    for (uint k = 0; k <= a[i].length - 1; k++) {
      set_m(a[i].x + k, a[i].y + k, 'o');
    }
  }

  cout << "---" << endl;
  for (uint y = 0; y < ylen; y++) {
    cout << M[y] << endl;
  }
  cout << "---" << endl;

  return 0;
}
