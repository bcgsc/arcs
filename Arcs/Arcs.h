#ifndef ARCS_H
#define ARCS_H 1

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <map>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <utility> 
#include <vector>
#include <iterator>
#include <time.h> 
//#include <boost/graph/graph_traits.hpp>
//#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graphviz.hpp>
#include "Common/Uncompress.h"
#include "DataLayer/FastaReader.h"
#include "DataLayer/FastaReader.cpp"


namespace ARCS {

    /**
     * Parameters controlling ARCS run
     */
    struct ArcsParams {

        std::string file;
        std::string fofName;
        int seq_id;
        int min_reads;
        int min_links;
        std::string base_name;
        std::string original_file;
        int min_mult;
        int max_mult;
        int max_grpSize;
        int max_degree;
        int end_length;
        int indexLen;
        int verbose;

        ArcsParams() : file(), fofName(), seq_id(90), min_reads(2), min_links(5), 
        base_name(""), min_mult(1000), max_mult(2000), 
        max_grpSize(100), max_degree(0), indexLen(14), verbose(0) {}

    };

    /* Scaffold counts: <pair(scaffold, bool), count>, cout =  # times index maps to scaffold (c), bool = true-head, false-tail*/
    typedef std::map<std::pair<int, bool>, int> ScafMap;
    /* indexMap: key = index sequence, value = scaffold counts */
    typedef std::unordered_map<std::string, ScafMap> IndexMap; 
    /* PairMap: key = pair(first < second) of scaf sequences, value = num links*/
    typedef std::map<std::pair<int, int>, std::vector<int>> PairMap; 

    struct VertexProperties {
        int id;
    };

    struct EdgeProperties {
        int orientation;
        int weight;
        EdgeProperties(): orientation(0), weight(0) {}
    };

	// Define the type of the graph - this specifies a bundled property for vertices
	typedef boost::undirected_graph<VertexProperties, EdgeProperties> Graph;
    typedef std::map<Graph::vertex_descriptor, std::size_t> VertexDescMap;
    typedef std::unordered_map<int, Graph::vertex_descriptor> VidVdesMap;
    typedef boost::graph_traits<ARCS::Graph>::vertex_descriptor VertexDes;
}

#endif
