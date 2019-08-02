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
#include "Common/ReadsProcessor.h"
// using sparse hash maps for k-merization
#include <google/sparse_hash_map>
#include "city.h"


namespace ARCS {

    /** value to use for 'd' in ABySS dist.gv */
    enum DistMode { DIST_MEDIAN=0, DIST_UPPER };

    /**
     * Parameters controlling ARCS run
     */
    struct ArcsParams {

        std::string file;
        std::string fofName;
        int seq_id;
        int min_reads;
        /** enable/disable distance estimation on graph edges */
        bool dist_est;
        /** bin size when computing distance estimates */
        unsigned dist_bin_size;
        /** output path for intra-contig distance/barcode samples (TSV) */
        std::string dist_samples_tsv;
        /** output path for inter-contig distance estimates (TSV) */
        std::string dist_tsv;
        /** chooses median or upper bound for `d` in ABySS dist.gv */
        DistMode dist_mode;
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
        std::string program;
        std::string multfile; 
        int k_value;
        double j_index;
	    unsigned threads;
        /** arks=true when user runs arks method. */
        bool arks;

        ArcsParams() :
            seq_id(98),
            min_reads(5),
            dist_est(false),
            dist_bin_size(20),
            dist_mode(DIST_MEDIAN),
            min_links(0),
            min_size(500),
            gap(100),
            min_mult(50),
            max_mult(10000),
            max_degree(0),
            end_length(30000),
            error_percent(0.05),
            verbose(0), 
            k_value(30), 
            j_index(0.55),
            threads(1),
            arks(false) {
        }

    };
    /** A contig end: (FASTA ID, head?) */
    typedef std::pair<std::string, bool> CI;

    //typedef const char* Kmer;

    /* ScafMap: <pair(scaffold id, bool), count>, cout =  # times index maps to scaffold (c), bool = true-head, false-tail*/
    typedef std::map<CI, int> ScafMap;
    typedef typename ScafMap::const_iterator ScafMapConstIt;
    /* IndexMap: key = index sequence, value = ScafMap */
    typedef std::unordered_map<std::string, ScafMap> IndexMap;
    /* PairMap: key = pair(first < second) of scaf sequence id, value = num links*/
    typedef std::map<std::pair<std::string, std::string>, std::vector<unsigned>> PairMap;

    

    /** a pair of contig IDs */
    typedef std::pair<std::string, std::string> ContigPair;

    /**
     * a list of the input scaffolds and their lengths, in the order
     * that they appear in the input contigs FASTA file
     */
    //typedef std::vector< std::pair<std::string, int> > ScaffSizeList;   // only really used in createABySS graph for iteration map can do the work easily.

    /**
     * a list of the input scaffolds and their lengths, in the order
     * that they appear in the input contigs FASTA file
     */
    typedef std::unordered_map<std::string, int> ScaffSizeMap;

    /** maps contig FASTA ID to contig length (bp) */
    typedef std::unordered_map<std::string, int> ContigToLength;
    typedef typename ContigToLength::const_iterator ContigToLengthIt;

    //* Hashmap data structure

    /* ContigKMap: <k-mer, pair(contig id, bool), hash<k-mer>, eqstr>
     * 	k-mer = string sequence
     *  contig id = string
     *  bool = True for Head; False for Tail
     *  eqstr = equal key
     */

    // simple hash adapter for types without pointers
    template<typename T>
    struct CityHasher {
	    size_t operator()(const T& t) const {
		    return CityHash64(t, sizeof(t));
	    }
    };

    // specialization for strings
    template<>
    struct CityHasher<std::string> {
	    size_t operator()(const string t) const {
		    return CityHash64(t.c_str(), t.size());
	    }
    };

    struct eqstr {
	    bool operator()(std::string s1, std::string s2) const {
		    return (s1 == s2);
	    }
    };
    
    typedef google::sparse_hash_map<std::string, int, CityHasher<std::string>, eqstr> ContigKMap;

    struct VertexProperties {
        std::string id;
    };

    /* Orientation: 0-HH, 1-HT, 2-TH, 3-TT */
    struct EdgeProperties {
        int orientation;
        int weight;
        int minDist;
        int dist;
        int maxDist;
        float jaccard;
        EdgeProperties() :
            orientation(0), weight(0),
            minDist(std::numeric_limits<int>::min()),
            dist(std::numeric_limits<int>::max()),
            maxDist(std::numeric_limits<int>::max()),
            jaccard(-1.0f)
            {}
    };

    template <class GraphT>
    struct EdgePropertyWriter
    {
        typedef typename boost::graph_traits<GraphT>::edge_descriptor E;
        typedef typename boost::edge_property<GraphT>::type EP;

        GraphT& m_g;

        EdgePropertyWriter(GraphT& g) : m_g(g) {}

        void operator()(std::ostream& out, const E& e) const
        {
            EP ep = m_g[e];
            out << '['
                << "label=" << ep.orientation << ", "
                << "weight=" << ep.weight;

            if (ep.minDist != std::numeric_limits<int>::min()) {
                assert(ep.dist != std::numeric_limits<int>::max());
                assert(ep.maxDist != std::numeric_limits<int>::max());
                assert(ep.jaccard >= 0.0f);
                out << ", "
                    << "d=" << ep.dist << ", "
                    << "maxd=" << ep.maxDist;
            }
            out << ']';
        }
    };

    template <class GraphT>
    struct VertexPropertyWriter
    {
        typedef typename boost::graph_traits<GraphT>::vertex_descriptor V;
        typedef typename boost::vertex_property<GraphT>::type VP;

        GraphT& m_g;

        VertexPropertyWriter(GraphT& g) : m_g(g) {}
        void operator()(std::ostream& out, const V& v) const
        {
            out << " [id=" << m_g[v].id << "]";
        }
    };

	typedef boost::undirected_graph<VertexProperties, EdgeProperties> Graph;
    typedef std::unordered_map<std::string, Graph::vertex_descriptor> VidVdesMap;
    typedef boost::graph_traits<ARCS::Graph>::vertex_descriptor VertexDes;
}

#endif
