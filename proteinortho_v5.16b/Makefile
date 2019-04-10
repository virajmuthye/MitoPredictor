
INSTALLDIR=/usr/local/bin

CPP = g++
CPPFLAGS += -Wall -O2 -std=gnu++11 -Wno-unused-result

all: proteinortho5_clustering proteinortho5_tree

proteinortho5_clustering: proteinortho5_clustering.cpp
	$(CPP) $(CPPFLAGS) $(LDFLAGS) -o $@ $<

proteinortho5_tree: po_tree.c
	$(CPP) $(CPPFLAGS) $(LDFLAGS) -o $@ $<

install: proteinortho5.pl proteinortho5_clustering proteinortho5_singletons.pl proteinortho5_clean_edges2.pl ffadj_mcs.py po2tree.pl proteinortho5_tree
	install -v $^ $(INSTALLDIR)

test: proteinortho5.pl proteinortho5_clustering
	./proteinortho5.pl -project=test -synteny -singles test/*.faa
	@./chk_test.pl test.proteinortho
	@./chk_test.pl test.poff
	rm test.*
	@echo "Test okay"
