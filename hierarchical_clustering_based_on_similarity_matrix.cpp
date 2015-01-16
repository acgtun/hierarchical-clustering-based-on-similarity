#include "hierarchical_clustering_based_on_similarity_matrix.hpp"

#include <string>
using std::string;

namespace clustering {

void HierarchicalClusteringBasedonSimilarityMatrix::GetDigits(
    const size_t& num, vector<size_t>& digits) {
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

double HierarchicalClusteringBasedonSimilarityMatrix::GetDistance(
    const size_t& a, const size_t& b) {
  vector < size_t > av(KMER);
  vector < size_t > bv(KMER);
  GetDigits(a, av);
  GetDigits(b, bv);

  double dis = 0.0;
  for (size_t i = 0; i < KMER; i++) {
    dis += pairwise_sim[av[i]][bv[i]];
  }
  return dis;
}

double HierarchicalClusteringBasedonSimilarityMatrix::GetSimilarityScore(
    const Cluster& a, const Cluster& b) {
  double avg = 0.0;
  for (size_t p = 0; p < a.members.size(); p++) {
    for (size_t q = 0; q < b.members.size(); q++) {
      avg += GetDistance(a.members[p], b.members[q]);
    }
  }
  avg /= (a.members.size() * b.members.size());

  return avg;
}

void HierarchicalClusteringBasedonSimilarityMatrix::HierarchicalClustering() {
  size_t r = num_of_points;
  size_t cluster_id = r;
  for (size_t t = 1; t <= r - 1; t++) {
    double max_sim = std::numeric_limits<int>::min();
    size_t k1 = 0, k2 = 0;
    for (size_t i = 0; i < r; i++) {
      for (size_t j = 0; j < r; j++) {
        if (i != j && clusters[i].indicator && clusters[j].indicator) {
          double sim = GetSimilarityScore(clusters[i], clusters[j]);
          if (sim > max_sim) {
            max_sim = sim;
            k1 = i;
            k2 = j;
          }
        }
      }
    }

    cout << "Iteration " << t << ": " << clusters[k1].id << " "
        << clusters[k2].id << " "
        << GetSimilarityScore(clusters[k1], clusters[k2]) << endl;
    /*merge k2 to k1*/
    for (size_t j = 0; j < clusters[k2].members.size(); j++) {
      clusters[k1].members.push_back(clusters[k2].members[j]);
    }
    clusters[k2].indicator = false;
    clusters[k2].members.clear();
    clusters[k1].id = cluster_id;
    cluster_id++;

    int c = 1;
    for (size_t i = 0; i < r; i++) {
      if (!clusters[i].indicator)
        continue;
      cout << "c" << c++ << ":";
      for (size_t j = 0; j < clusters[i].members.size(); j++) {
        //cout << " " << AA20[clusters[i].members[j]];
        cout << " " << clusters[i].members[j];
      }
      cout << endl;
    }
  }
}

} /* namespace clustering */
