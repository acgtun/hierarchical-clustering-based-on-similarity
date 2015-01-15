#ifndef HIERARCHICALCLUSTERINGBASEDONSIMILARITYMATRIX_H_
#define HIERARCHICALCLUSTERINGBASEDONSIMILARITYMATRIX_H_

#include <vector>
#include <limits>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

namespace clustering {

struct Cluster {
  bool indicator;
  size_t id;
  vector<size_t> members;
  vector<double> similarity;
};

class HierarchicalClusteringBasedonSimilarityMatrix {
 public:
  HierarchicalClusteringBasedonSimilarityMatrix(
      const vector<vector<double> >& _similarity_matrix)
      : similarity_matrix(_similarity_matrix) {
    num_of_points = similarity_matrix.size();
    clusters.resize(num_of_points);
    for (size_t i = 0; i < num_of_points; ++i) {
      clusters[i].indicator = true;
      clusters[i].id = i;
      clusters[i].members.push_back(i);
      clusters[i].similarity.resize(num_of_points);
    }
    for (size_t i = 0; i < num_of_points; ++i) {
      for (size_t j = 0; j <= i; ++j) {
        clusters[i].similarity[j] = similarity_matrix[i][j];
        clusters[j].similarity[i] = similarity_matrix[j][i];
      }
    }

    HierarchicalClustering();
  }
  ~HierarchicalClusteringBasedonSimilarityMatrix() {
  }

 private:
  void HierarchicalClustering();

  vector<vector<double> > similarity_matrix;
  vector<Cluster> clusters;
  size_t num_of_points;
};

} /* namespace clustering */

#endif /* HIERARCHICALCLUSTERINGBASEDONSIMILARITYMATRIX_H_ */
