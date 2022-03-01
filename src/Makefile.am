bin_PROGRAMS = long-to-linked-pe

long_to_linked_pe_CXXFLAGS = -O3 -Wall -Werror -Wno-unknown-pragmas -std=c++11 -pthread -I$(top_srcdir)/src/btllib/include

CXXFLAGS =  $(long_to_linked_pe_CXXFLAGS)

long_to_linked_pe_SOURCES = long-to-linked-pe.cpp \
	include/btllib/seq_reader.hpp \
	include/btllib/status.hpp \
	include/btllib/util.hpp \
	include/btllib/cstring.hpp \
	include/btllib/data_stream.hpp \
	include/btllib/order_queue.hpp \
	include/btllib/seq.hpp

# Check the C++ source code for errors.
lint: clang-format clang-tidy

# Check the C++ source code for errors with clang-tidy.
clang-tidy:
	clang-tidy -warnings-as-errors='*' *.cpp -- -std=c++11 -x c++ -Ibtllib/include -Iinclude

# Check the C++ source code for white-space errors with clang-format.
clang-format:
	for i in *.cpp; do clang-format -style=file $$i >$$i.fixed; done
	for i in *.cpp; do diff -su $$i $$i.fixed && rm -f $$i.fixed; done
	if ls *.fixed; then exit 1; fi
