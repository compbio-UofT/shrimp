#include <iostream>
#include <cstring>

using namespace std;

extern "C" {

#include "../common/util.h"

}

int baseFromChar(char c) {
  switch (c) {
  case 'a':
  case 'A':
    return BASE_A;
  case 'c':
  case 'C':
    return BASE_C;
  case 'g':
  case 'G':
    return BASE_G;
  case 't':
  case 'T':
    return BASE_T;
  }

  return -1;
}


char charFromBase(int b) {
  switch (b) {
  case BASE_A:
    return 'a';
  case BASE_C:
    return 'c';
  case BASE_G:
    return 'g';
  case BASE_T:
    return 't';
  }

  return -1;
}


int main() {
  char buffer[10000];

  cin >> buffer;

  int n = strlen(buffer);
  cout << 'T' << lstocs(BASE_T, baseFromChar(buffer[0]), false);
  for (int i = 1; i < n; i++) {
    cout << lstocs(baseFromChar(buffer[i-1]), baseFromChar(buffer[i]), false);
  }

  cout << endl;

  return 0;
}
