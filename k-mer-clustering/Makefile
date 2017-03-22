PROGS = k-mer-cluster

CXX = g++
CFLAGS = -Wall -O3 -fPIC -fmessage-length=50

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

SRC = hierarchical_clustering_based_on_similarity_matrix.cpp
OBJ = $(patsubst %.cpp,%.o,$(SRC))

all:	$(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $< 

k-mer-cluster : $(OBJ) 

%: %.cpp
		$(CXX) $(CFLAGS) -o $@ $^ 

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
