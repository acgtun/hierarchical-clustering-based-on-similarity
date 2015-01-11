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

#include "reduced_alphabet.h"

using namespace std;

const size_t KMER = 10;

struct multisetcmp {
  bool operator()(const pair<double, int> & a,
                  const pair<double, int> & b) const {
    return a.first < b.first;
  }
};

struct Cluster {
  bool indicator;
  vector<size_t> cluster_ids;
  vector<double> pairwise_distance;
};

void GetDigits(const size_t& num, vector<size_t>& digits) {
  size_t n = num;
  size_t i = 0;
  while (n) {
    digits[i++] = n % 5;
    n /= 5;
  }
  while (i < KMER) {
    digits[i++] = 0;
  }
}
double GetDistance(const size_t& a, const size_t& b,
                   const vector<vector<double> >& pairwise_dis) {
  vector<size_t> av(KMER);
  vector<size_t> bv(KMER);
  GetDigits(a, av);
  GetDigits(b, bv);

  double dis = 0.0;
  for (size_t i = 0; i < KMER; i++) {
    dis += pairwise_dis[av[i]][bv[i]];
  }
  return dis;
}

void HierarchicalClustering(vector<Cluster>& clusters,
                            const vector<vector<double> >& pairwise_dis) {
  cout << "Hierarchical Clustering..." << endl;
  ofstream fout("test.txt");
  size_t r = clusters.size();
  for (size_t t = 1; t <= r - 1; t++) {
    cout << "t = " << t << " " << r - 1 << endl;
    double min_dis = std::numeric_limits<double>::max();
    size_t k1 = 0, k2 = 0;
    for (size_t i = 0; i < r; i++) {
      for (size_t j = 0; j < i; j++) {
        if (clusters[i].indicator && clusters[j].indicator
            && clusters[i].pairwise_distance[j] < min_dis) {
          min_dis = clusters[i].pairwise_distance[j];
          k1 = i;
          k2 = j;
        }
      }
    }

    /*merge k2 to k1*/
    for (size_t j = 0; j < clusters[k2].cluster_ids.size(); j++) {
      clusters[k1].cluster_ids.push_back(clusters[k2].cluster_ids[j]);
    }
    clusters[k2].indicator = false;
    clusters[k2].cluster_ids.clear();
    clusters[k2].pairwise_distance.clear();

    for (size_t j = 0; j < r; j++) {
      if (clusters[j].indicator && j != k1) {
        double max_dis = std::numeric_limits<double>::min();
        for (size_t p = 0; p < clusters[k1].cluster_ids.size(); p++) {
          for (size_t q = 0; q < clusters[j].cluster_ids.size(); q++) {
            double d = GetDistance(clusters[k1].cluster_ids[p],
                                   clusters[j].cluster_ids[q], pairwise_dis);
            if (d > max_dis) {
              d = max_dis;
            }
          }
        }
        clusters[j].pairwise_distance[k1] = max_dis;
        clusters[k1].pairwise_distance[j] = max_dis;
      }
    }
    fout << "t = " << t << " " << r - 1 << endl;
    for (size_t i = 0; i < r; i++) {
      if (clusters[i].cluster_ids.size() == 0)
        continue;
      fout << "c " << i << ":";
      for (size_t j = 0; j < clusters[i].cluster_ids.size(); j++) {
        fout << " " << clusters[i].cluster_ids[j];
      }
      fout << endl;
    }
  }
  fout.close();
}

int main() {
  vector<vector<double> > pairwise_dis;
  GetReducedAphabetDistance(pairwise_dis);

  size_t r = static_cast<size_t>(pow(5, 6));

  vector<Cluster> clusters(r);
  for (size_t i = 0; i < r; i++) {
    clusters[i].indicator = true;
    clusters[i].cluster_ids.push_back(i);
    clusters[i].pairwise_distance.resize(r);
  }

  for (size_t i = 0; i < r; i++) {
    for (size_t j = 0; j < i; j++) {
      double d = GetDistance(i, j, pairwise_dis);
      clusters[i].pairwise_distance[j] = d;
      clusters[j].pairwise_distance[i] = d;
    }
  }
  HierarchicalClustering(clusters, pairwise_dis);
}
