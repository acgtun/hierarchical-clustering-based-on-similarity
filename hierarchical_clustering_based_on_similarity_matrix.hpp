#ifndef HIERARCHICALCLUSTERINGBASEDONSIMILARITYMATRIX_H_
#define HIERARCHICALCLUSTERINGBASEDONSIMILARITYMATRIX_H_

#include <tr1/unordered_map>
#include <cmath>
#include <vector>
#include <limits>
#include <string>
#include <iostream>

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::tr1::unordered_map;

namespace clustering {

#define ALPHABETSIZE 8
#define KMER 5

const double pairwise_sim[ALPHABETSIZE][ALPHABETSIZE] = { { 1.88889, -0.8, -1,
    -0.666667, -0.666667, -1.08333, -2.22222, -1 }, { -0.8, 1.52, -0.1, -3.2,
    -1.8, -2.35, -2.66667, -1.2 }, { -1, -0.1, 4, -3, -1, -2.75, -1.66667, -2 },
    { -0.666667, -3.2, -3, 9, -3, -1, -2, -3 }, { -0.666667, -1.8, -1, -3, 6,
        -3.5, -2.66667, -2 }, { -1.08333, -2.35, -2.75, -1, -3.5, 2.3125,
        -1.16667, -2.5 }, { -2.22222, -2.66667, -1.66667, -2, -2.66667,
        -1.16667, 4, -3.66667 }, { -1, -1.2, -2, -3, -2, -2.5, -3.66667, 7 } };

struct Cluster {
  bool indicator;
  size_t id;
  vector<size_t> members;
  vector<double> similarity;
};

class HierarchicalClusteringBasedonSimilarityMatrix {
 public:
  HierarchicalClusteringBasedonSimilarityMatrix() {
    num_of_points = static_cast<size_t>(pow(ALPHABETSIZE, KMER));
    BuildID2Kmer();
    clusters.resize(num_of_points);
    for (size_t i = 0; i < num_of_points; ++i) {
      clusters[i].indicator = true;
      clusters[i].id = i;
      clusters[i].members.push_back(i);
    }
    for (size_t i = 0; i < num_of_points; ++i) {
      clusters[i].similarity.resize(num_of_points);
      for (size_t j = 0; j < num_of_points; ++j) {
        clusters[i].similarity[j] = ClusterSimilarityScore(clusters[i],
                                                           clusters[j]);
      }
    }

    HierarchicalClustering();
  }
  ~HierarchicalClusteringBasedonSimilarityMatrix() {
  }

 private:
  void HierarchicalClustering();
  void BuildID2Kmer();
  double KmerSimilarityScore(const size_t& a, const size_t& b);
  double ClusterSimilarityScore(const Cluster& a, const Cluster& b);

  vector<Cluster> clusters;
  size_t num_of_points;
  unordered_map<size_t, string> id_kmer;
  vector<vector<double> > similarity_matrix;
};

} /* namespace clustering */

#endif /* HIERARCHICALCLUSTERINGBASEDONSIMILARITYMATRIX_H_ */
