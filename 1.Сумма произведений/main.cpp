#include <iostream>

int rec(int iterator, int* used, int k, int* lens, int** elems) {
  if (iterator == k) return 1;
  else {
    int s = 0;
    for (int j = 0; j < lens[iterator]; j++) {
      int flag = 0;
      for (int i = 0; i < iterator; i++) {
        if (j == used[i]) {
          flag = 1;
        }
      }
      if (flag == 0) {
        used[iterator] = j;
        s += elems[iterator][j] * rec(iterator + 1, used, k, lens, elems);
      }
    }
    return s;
  }
}

int main(int argc, char* argv[]) {
  int k = argc - 1;
  int** elems = new int* [k];
  int* lens = new int[k];
  for (int i = 0; i < k; i++) {
    int n = atoi(argv[i + 1]);
    lens[i] = n;
    elems[i] = new int[n];
    for (int j = 0; j < n; j++) {
      std::cin >> elems[i][j];
    }
  }
  int* used = new int[k];
  int ans = rec(0, used, k, lens, elems);
  std::cout << ans;
  delete[] used;
  delete[] lens;
  for (int i = 0; i < k; i++) {
    delete[] elems[i];
  }
  delete[] elems;
  return 0;
}
