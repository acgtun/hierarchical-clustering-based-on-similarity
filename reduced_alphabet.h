#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <vector>
#include <iostream>


using namespace std;

const int BLOSUM62[][20] = {
//A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0 },  //A
{ -1, 5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3 },  //R
{ -2, 0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3 },  //N
{ -2,-2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3 },  //D
{  0,-3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1 },  //C
{ -1, 1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2 },  //Q
{ -1, 0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2 },  //E
{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3 },  //G
{ -2, 0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3 },  //H
{ -1,-3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3 },  //I
{ -1,-2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1 },  //L
{ -1, 2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2 },  //K
{ -1,-1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1 },  //M
{ -2,-3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1 },  //F
{ -1,-2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2 },  //P
{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2 },  //S
{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0 },  //T
{ -3,-3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3 },  //W
{ -2,-2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1 },  //Y
{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3,  -1, 4 } };  //V

void GetMetricDistance(vector<vector<double> >& distance) {
  double lamda = -0.17;
  int cnt = 0;
  for (int i = 0; i < 20; i++) {
    for (int j = 0; j < 20; j++) {
      for (int k = 0; k < 20; k++) {
        if (exp(lamda * BLOSUM62[i][j]) + exp(lamda * BLOSUM62[j][k])
            < exp(lamda * BLOSUM62[i][k])) {
          cnt++;
        } else {
        }
      }
    }
  }
  if (cnt > 0) {
    cerr << "error~!~" << endl;
    exit (EXIT_FAILURE);
  }

  distance.resize(20);
  for (int i = 0; i < 20; i++) {
    distance[i].resize(20);
    for (int j = 0; j < 20; j++) {
      if (i == j) {
        distance[i][j] = 0.000;
      } else {
        distance[i][j] = exp(lamda * BLOSUM62[i][j]);
      }
    }
  }
}

void GetReducedAphabetDistance(vector<vector<double> >& pairwise_dis) {
  vector<vector<double> > distance;
  GetMetricDistance(distance);
  vector<vector<int> > cluster(5);
  cluster[0].push_back(9);
  cluster[0].push_back(19);
  cluster[0].push_back(10);
  cluster[0].push_back(12);
  cluster[0].push_back(4);

  cluster[1].push_back(13);
  cluster[1].push_back(18);
  cluster[1].push_back(17);

  cluster[2].push_back(0);
  cluster[2].push_back(15);
  cluster[2].push_back(16);
  cluster[2].push_back(14);

  cluster[3].push_back(1);
  cluster[3].push_back(11);
  cluster[3].push_back(5);
  cluster[3].push_back(6);
  cluster[3].push_back(8);

  cluster[4].push_back(2);
  cluster[4].push_back(3);
  cluster[4].push_back(7);

  pairwise_dis.resize(5);
  for (int i = 0; i < 5; i++) {
    pairwise_dis[i].resize(5);
    for (int j = 0; j < 5; j++) {
      pairwise_dis[i][j] = 0.0;
      if(i == j) continue;
      for (size_t p = 0; p < cluster[i].size(); p++) {
        for (size_t q = 0; q < cluster[j].size(); q++) {
          if (distance[cluster[i][p]][cluster[j][q]] > pairwise_dis[i][j]) {
            pairwise_dis[i][j] = distance[cluster[i][p]][cluster[j][q]];
          }
        }
      }
    }
  }

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      cout << pairwise_dis[i][j] << " ";
    }
    cout << endl;
  }
}

