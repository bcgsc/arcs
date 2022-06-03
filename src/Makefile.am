bin_PROGRAMS = long-to-linked-pe

long_to_linked_pe_CXXFLAGS = -O3 -Wall -Werror -Wno-unknown-pragmas -std=c++17 -pthread -I$(top_srcdir)/src/include 

CXXFLAGS += $(long_to_linked_pe_CXXFLAGS)

long_to_linked_pe_LDADD = -lbtllib

long_to_linked_pe_SOURCES = long-to-linked-pe.cpp 

# Check the C++ source code for errors.
lint: clang-format clang-tidy

# Check the C++ source code for errors with clang-tidy.
clang-tidy:
	clang-tidy -warnings-as-errors='*' *.cpp -- -std=c++11 -x c++  ${CXXFLAGS}

# Check the C++ source code for white-space errors with clang-format.
clang-format:
	for i in *.cpp; do clang-format -style=file $$i >$$i.fixed; done
	for i in *.cpp; do diff -su $$i $$i.fixed && rm -f $$i.fixed; done
	if ls *.fixed; then exit 1; fi
