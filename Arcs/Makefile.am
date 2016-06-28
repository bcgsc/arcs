bin_PROGRAMS = arcs

arcs_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

arcs_CPPFLAGS = -I$(top_srcdir)/Arcs \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)
	
arcs_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a -lz ~/lib/libboost_graph.a ~/lib/libboost_regex.a
	
arcs_LDFLAGS = $(OPENMP_CXXFLAGS)

arcs_SOURCES = Arcs.h Arcs.cpp
