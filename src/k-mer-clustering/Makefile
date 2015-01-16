ifndef SRC_ROOT
SRC_ROOT=./../
endif

ifndef Hierarchical_Clustering
Hierarchical_Clustering=./../hierarchical_clustering
endif

PROGS = k-mer-cluster

CXX = g++
CFLAGS = -Wall -O3 -fPIC -fmessage-length=50

TEST_DIR = $(SRC_ROOT)/test

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

SRC = 
OBJ = $(patsubst %.cpp,%.o,$(SRC))

all:	$(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $<  -I$(Hierarchical_Clustering)

k-mer-cluster : $(OBJ) \
	$(addprefix $(Hierarchical_Clustering)/, hierarchical_clustering_based_on_similarity_matrix.o)

%: %.cpp
		$(CXX) $(CFLAGS) -o $@ $^  -I$(Hierarchical_Clustering)
	
install: all
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

test_%:	%
	@$(TEST_DIR)/$@ $(TEST_DIR)

test:	$(addprefix test_, $(PROGS))

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
