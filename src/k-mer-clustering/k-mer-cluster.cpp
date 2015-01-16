#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <set>
#include <limits>
#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

#define ALPHABETSIZE 8
#define KMER 6

double pairwise_sim[ALPHABETSIZE][ALPHABETSIZE] = { { 1.88889, -0.8, -1,
    -0.666667, -0.666667, -1.08333, -2.22222, -1 }, { -0.8, 1.52, -0.1, -3.2,
    -1.8, -2.35, -2.66667, -1.2 }, { -1, -0.1, 4, -3, -1, -2.75, -1.66667, -2 },
    { -0.666667, -3.2, -3, 9, -3, -1, -2, -3 }, { -0.666667, -1.8, -1, -3, 6,
        -3.5, -2.66667, -2 }, { -1.08333, -2.35, -2.75, -1, -3.5, 2.3125,
        -1.16667, -2.5 }, { -2.22222, -2.66667, -1.66667, -2, -2.66667,
        -1.16667, 4, -3.66667 }, { -1, -1.2, -2, -3, -2, -2.5, -3.66667, 7 } };

void GetDigits(const size_t& num, vector<size_t>& digits) {
  size_t n = num;
  size_t i = 0;
  while (n) {
    digits[i++] = n % ALPHABETSIZE;
    n /= ALPHABETSIZE;
  }
  while (i < KMER) {
    digits[i++] = 0;
  }
}

double GetDistance(const size_t& a, const size_t& b) {
  vector<size_t> av(KMER);
  vector<size_t> bv(KMER);
  GetDigits(a, av);
  GetDigits(b, bv);

  double dis = 0.0;
  for (size_t i = 0; i < KMER; i++) {
    dis += pairwise_sim[av[i]][bv[i]];
  }
  return dis;
}

int main() {
  size_t r = static_cast<size_t>(pow(ALPHABETSIZE, KMER));
  for(size_t i = 0;i < r;i++) {
    vector<size_t> av(KMER);
    GetDigits(i, av);
    for(size_t j = 0;j < KMER;j++) {
      cout << av[j];
    }
    cout << endl;
  }
  vector<vector<double> > sim_matrix(r, vector<double>(r, 0));
  for(size_t i = 0;i < r;i++) {
    for(size_t j = 0;j <= i;j++) {
      double sim = GetDistance(i, j);
      sim_matrix[i][j] = sim;
      sim_matrix[j][i] = sim;
    }
  }

}
