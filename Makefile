
CXX := g++
CXXFLAGS := -O3 -Wall -fmessage-length=50

SRCS := hierarchical_clustering_based_on_similarity_matrix.cpp 
OBJS := $(SRCS:.cc=.o)

.cc.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	
hierarchical_clustering_based_on_similarity_matrix: hierarchical_clustering_based_on_similarity_matrix_main.o $(OBJS)
	 $(CXX) $(CXXFLAGS) -o $@ $^
                

all: hierarchical_clustering_based_on_similarity_matrix

clean:
	rm -rf hierarchical_clustering_basedon_similarity_matrix  *.exe *.o  
