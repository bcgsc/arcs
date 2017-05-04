#ifndef ARCS_H
#define ARCS_H 1

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <iostream>
#include <utility>
#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <utility> 
#include <vector>
#include <iterator>
#include <time.h> 
#include <boost/graph/undirected_graph.hpp>
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
        int min_size;
        std::string base_name;
        std::string dist_graph_name;
        std::string tsv_name;
        std::string barcode_counts_name;
        unsigned gap;
        int min_mult;
        int max_mult;
        int max_degree;
        int end_length;
        float error_percent;
        int verbose;

        ArcsParams() :
            seq_id(98),
            min_reads(5),
            min_links(0),
            min_size(500),
            gap(100),
            min_mult(50),
            max_mult(10000),
            max_degree(0),
            end_length(0),
            error_percent(0.05),
            verbose(0) {
        }

    };

    /** One segment of a scaffold. */
    typedef std::pair<std::string, unsigned> Segment;

    /** Hash a Segment. */
    struct HashSegment {
        size_t operator()(const Segment& key) const {
            return std::hash<std::string>()(key.first) ^ std::hash<unsigned>()(key.second);
        }
    };

    /* ScafMap: <pair(scaffold id, bool), count>, count = Number of times index maps to scaffold */
    typedef std::map<std::pair<std::string, unsigned>, int> ScafMap;
    /* IndexMap: key = index sequence, value = ScafMap */
    typedef std::unordered_map<std::string, ScafMap> IndexMap; 
    /* PairMap: key = pair of scaf sequence id, value = num links*/
    typedef std::map<std::pair<Segment, Segment>, unsigned> PairMap;

    struct VertexProperties {
        std::string id;
    };

    struct EdgeProperties {
        int weight;
        EdgeProperties(): weight(0) {}
    };

	typedef boost::undirected_graph<VertexProperties, EdgeProperties> Graph;
    typedef std::unordered_map<std::string, Graph::vertex_descriptor> VidVdesMap;
    typedef boost::graph_traits<ARCS::Graph>::vertex_descriptor VertexDes;
}

#endif
