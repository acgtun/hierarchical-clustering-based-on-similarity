#include "hierarchical_clustering_based_on_similarity_matrix.h"

#include <string>
using std::string;

namespace clustering {

const string AA20 = "ARNDCEQGHILKMFPSTWYV";

void HierarchicalClusteringBasedonSimilarityMatrix::HierarchicalClustering() {
  size_t r = num_of_points;
  size_t cluster_id = r;
  for (size_t t = 1; t <= r - 1; t++) {
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

    for (size_t j = 0; j < r; j++) {
      if (clusters[j].indicator) {
        double avg = 0.0;
        for (size_t p = 0; p < clusters[k1].members.size(); p++) {
          for (size_t q = 0; q < clusters[j].members.size(); q++) {
            avg += similarity_matrix[clusters[k1].members[p]][clusters[j]
                .members[q]];
          }
        }
        avg /= (clusters[k1].members.size() * clusters[j].members.size());
        clusters[j].similarity[k1] = avg;
        clusters[k1].similarity[j] = avg;
      }
    }

    for (size_t i = 0; i < r; i++) {
      if (!clusters[i].indicator)
        continue;
      for (size_t j = 0; j < r; j++) {
        if (!clusters[j].indicator)
          continue;
        cout << clusters[i].similarity[j] << "\t";
      }
      cout << endl;
    }

  }
}

} /* namespace clustering */
