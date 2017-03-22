PROGS = hierarchical_clustering_based_on_similarity_matrix_main

CXX = g++
CFLAGS = -Wall -O3 -fPIC -fmessage-length=50


ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

all:	$(PROGS)

SRC = hierarchical_clustering_based_on_similarity_matrix.cpp
OBJ = $(patsubst %.cpp,%.o,$(SRC))

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $<

hierarchical_clustering_based_on_similarity_matrix_main :  $(OBJ)


%: %.cpp
		$(CXX) $(CFLAGS) -o $@ $^


clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean