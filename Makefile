K-MER-CLUSTERING = $(shell pwd)

all:
	@make -C src K-MER-CLUSTERING=$(K-MER-CLUSTERING) OPT=1

install:
	@make -C src K-MER-CLUSTERING=$(K-MER-CLUSTERING) OPT=1 install

test:
	@make -C src K-MER-CLUSTERING=$(K-MER-CLUSTERING) test
.PHONY: test

clean:
	@make -C src K-MER-CLUSTERING=$(K-MER-CLUSTERING) clean
.PHONY: clean

distclean: clean
	@rm -rf $(K-MER-CLUSTERING)/bin
	@rm -rf $(K-MER-CLUSTERING)/include
.PHONY: distclean
