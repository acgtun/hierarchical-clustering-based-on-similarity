ifndef Hierarchical_Clustering
Hierarchical_Clustering=$(K-MER-CLUSTERING)/src/hierarchical_clustering
endif

app_subdirs=k-mer-clustering

all:
	@for i in $(app_subdirs); do \
		make -C $${i} Hierarchical_Clustering=$(Hierarchical_Clustering) SRC_ROOT=$(K-MER-CLUSTERING) OPT=1; \
	done;

install:
	@for i in $(app_subdirs); do \
		make -C $${i} Hierarchical_Clustering=$(Hierarchical_Clustering) SRC_ROOT=$(K-MER-CLUSTERING) OPT=1 install; \
	done;

test:
	@for i in $(app_subdirs); do \
		make -C $${i} Hierarchical_Clustering=$(Hierarchical_Clustering) SRC_ROOT=$(K-MER-CLUSTERING) test; \
	done;

clean:
	@make -C hierarchical_clustering OPT=1 clean
	@for i in $(app_subdirs); do \
		make -C $${i} Hierarchical_Clustering=$(Hierarchical_Clustering) SRC_ROOT=$(K-MER-CLUSTERING) clean; \
	done;
