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

  cout << "similarity matrix..." << endl;
  similarity_matrix.resize(num_of_points);
  for (size_t i = 0; i < num_of_points; ++i) {
    similarity_matrix[i].resize(num_of_points);
    for (size_t j = 0; j < num_of_points; ++j) {
      similarity_matrix[i][j] = KmerSimilarityScore(i, j);
    }
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
      avg += similarity_matrix[a.members[p]][b.members[q]];
    }
  }
  avg /= (a.members.size() * b.members.size());

  return avg;
}

void HierarchicalClusteringBasedonSimilarityMatrix::HierarchicalClustering() {
  size_t r = num_of_points;
  size_t cluster_id = r;
  for (size_t t = 1; t <= r - 1; t++) {
    cout << "t = " << t << endl;
    double max_sim = std::numeric_limits<int>::min();
    size_t k1 = 0, k2 = 0;
    for (size_t i = 0; i < r; i++) {
      for (size_t j = 0; j < r; j++) {
        if (i != j && clusters[i].indicator && clusters[j].indicator
            && clusters[i].similarity[j] > max_sim) {
          max_sim = clusters[i].similarity[j];
          k1 = i;
          k2 = j;
        }
      }
    }

    cout << "Iteration " << t << ": " << clusters[k1].id << " "
        << clusters[k2].id << " " << clusters[k1].similarity[k2] << endl;
    /*merge k2 to k1*/
    for (size_t j = 0; j < clusters[k2].members.size(); j++) {
      clusters[k1].members.push_back(clusters[k2].members[j]);
    }
    clusters[k2].indicator = false;
    clusters[k2].members.clear();
    clusters[k2].similarity.clear();
    clusters[k1].id = cluster_id;
    cluster_id++;

    for (size_t i = 0; i < r; ++i) {
      if (clusters[i].indicator) {
        double sim = ClusterSimilarityScore(clusters[i], clusters[k1]);
        clusters[i].similarity[k1] = sim;
        clusters[k1].similarity[i] = sim;
      }
    }

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
