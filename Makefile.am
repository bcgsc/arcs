SUBDIRS = \
	Common \
	DataLayer \
	Graph \
	Arcs \
	Test

dist_doc_DATA = \
	README.md

EXTRA_DIST = autogen.sh \
	     Examples/arcs-make \
	     Examples/long-to-linked-pe \
	     Examples/makeTSVfile.py

CXXFLAGS=-O3 -Wall -Werror -Wno-unknown-pragmas -std=c++11 -pthread -Isrc/btllib/include -Isrc/include

all: src/long-to-linked-pe                                                                                                  
clean:	rm -f src/long-to-linked-pe

src/long-to-linked-pe:
	$(CXX) $(CXXFLAGS) $@.cpp -o $@

