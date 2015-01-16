#include "hierarchical_clustering_based_on_similarity_matrix.hpp"

#include <string>
using std::string;

namespace clustering {

void HierarchicalClusteringBasedonSimilarityMatrix::BuildID2Kmer() {
  for (size_t i = 0; i < num_of_points; ++i) {
    string kmer;
    size_t n = i;
    size_t j = 0;
    while (n) {
      kmer += 48 + n % ALPHABETSIZE;
      j++;
      n /= ALPHABETSIZE;
    }
    while (j < KMER) {
      kmer += 48;
      j++;
    }
    id_kmer.insert(make_pair(i, kmer));
    cout << i << " " << kmer << endl;
  }
}

double HierarchicalClusteringBasedonSimilarityMatrix::KmerSimilarityScore(
    const size_t& a, const size_t& b) {
  string sa = id_kmer[a];
  string sb = id_kmer[b];

  double sim = 0.0;
  for (size_t i = 0; i < KMER; i++) {
    sim += pairwise_sim[sa[i] - 48][sb[i] - 48];
  }
  return sim;
}

double HierarchicalClusteringBasedonSimilarityMatrix::ClusterSimilarityScore(
    const Cluster& a, const Cluster& b) {
  double avg = 0.0;
  for (size_t p = 0; p < a.members.size(); p++) {
    for (size_t q = 0; q < b.members.size(); q++) {
      avg += KmerSimilarityScore(a.members[p], b.members[q]);
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
      if (!clusters[i].indicator)
        continue;
      for (size_t j = i + 1; j < r; j++) {
        if (!clusters[j].indicator)
          continue;
        double sim = ClusterSimilarityScore(clusters[i], clusters[j]);
        if (sim > max_sim) {
          max_sim = sim;
          k1 = i;
          k2 = j;
        }
      }
    }

    cout << "Iteration " << t << ": " << clusters[k1].id << " "
        << clusters[k2].id << " "
        << ClusterSimilarityScore(clusters[k1], clusters[k2]) << endl;
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
