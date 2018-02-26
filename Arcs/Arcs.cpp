#include "config.h"
#include "Arcs.h"
#include "Arcs/DistanceEst.h"
#include "Common/Barcode.h"
#include "Common/ContigProperties.h"
#include "Common/Estimate.h"
#include "Common/ParseUtil.h"
#include "Common/SAM.h"
#include "Common/Segment.h"
#include "Graph/ContigGraph.h"
#include "Graph/DirectedGraph.h"
#include "Graph/DotIO.h"
#include <algorithm>
#include <cassert>
#include <utility>

#define PROGRAM "arcs"

static const char VERSION_MESSAGE[] =
PROGRAM " " PACKAGE_VERSION "\n"
"\n"
"http://www.bcgsc.ca/platform/bioinfo/software/links \n"
"We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca.\n"
"If you use LINKS, ARCS code or ideas, please cite our work. \n"
"\n"
"LINKS and ARCS Copyright (c) 2014-2016 Canada's Michael Smith Genome Science Centre.  All rights reserved. \n";

static const char USAGE_MESSAGE[] =
PROGRAM " " PACKAGE_VERSION "\n"
"Usage: arcs [OPTION]... -f FASTA_FILE BAM_FILE...\n"
"\n"
"NOTE 1: BAM_FILE must be sorted in order of name\n"
"NOTE 2: read names in BAM_FILE must be formatted as <READ_NAME>_<BARCODE_SEQ>\n"
"        unless --bx is used\n"
"\n"
" Options:\n"
"\n"
"       --bx              extract barcode sequences from BX tag\n"
"   -f, --file=FILE       input sequences to scaffold [required]\n"
"   -a, --fofName=FILE    text file listing input BAM filenames\n"
"   -s, --seq_id=N        min sequence identity for read alignments [98]\n"
"   -c, --min_reads=N     min aligned read pairs per barcode mapping [5]\n"
"   -l, --min_links=N     min shared barcodes between contigs [0]\n"
"   -z, --min_size=N      min contig length [500]\n"
"   -b, --base_name=STR   output file prefix\n"
"   -g, --graph=FILE      write the ABySS dist.gv to FILE\n"
"       --gap=N           fixed gap size for ABySS dist.gv file [100]\n"
"       --tsv=FILE        write graph in TSV format to FILE\n"
"       --barcode-counts=FILE       write number of reads per barcode to FILE\n"
"   -m, --index_multiplicity=RANGE  barcode multiplicity range [50-10000]\n"
"   -d, --max_degree=N    max node degree in scaffold graph [0]\n"
"   -e, --end_length=N    contig head/tail length for masking alignments [30000]\n"
"   -r, --error_percent=N p-value for head/tail assignment and link orientation\n"
"                         (lower is more stringent) [0.05]\n"
"   -v, --run_verbose     verbose logging\n"
"\n"
" Distance Estimation Options:\n"
"\n"
"   -B, --bin_size=N        estimate distance using N closest Jaccard scores [20]\n"
"   -D, --dist_est          enable distance estimation\n"
"       --no_dist_est       disable distance estimation [default]\n"
"       --dist_median       use median distance in ABySS dist.gv [default]\n"
"       --dist_upper        use upper bound distance in ABySS dist.gv\n"
"       --dist_tsv=FILE     write min/max distance estimates to FILE\n"
"       --samples_tsv=FILE  write intra-contig distance/barcode samples to FILE\n";

static const char shortopts[] = "f:a:B:s:c:Dl:z:b:g:m:d:e:r:v";

enum {
    OPT_HELP = 1,
    OPT_VERSION,
    OPT_BX,
    OPT_GAP,
    OPT_TSV,
    OPT_BARCODE_COUNTS,
    OPT_SAMPLES_TSV,
    OPT_DIST_TSV,
    OPT_NO_DIST_EST,
    OPT_DIST_MEDIAN,
    OPT_DIST_UPPER
};

static const struct option longopts[] = {
    {"file", required_argument, NULL, 'f'},
    {"fofName", required_argument, NULL, 'a'},
    {"bin_size", required_argument, NULL, 'B'},
    {"bx", no_argument, NULL, OPT_BX },
    {"samples_tsv", required_argument, NULL, OPT_SAMPLES_TSV},
    {"dist_tsv", required_argument, NULL, OPT_DIST_TSV},
    {"seq_id", required_argument, NULL, 's'},
    {"min_reads", required_argument, NULL, 'c'},
    {"min_reads", required_argument, NULL, 'c'},
    {"dist_est", no_argument, NULL, 'D'},
    {"no_dist_est", no_argument, NULL, OPT_NO_DIST_EST},
    {"dist_median", no_argument, NULL, OPT_DIST_MEDIAN},
    {"dist_upper", no_argument, NULL, OPT_DIST_UPPER},
    {"min_links", required_argument, NULL, 'l'},
    {"min_size", required_argument, NULL, 'z'},
    {"base_name", required_argument, NULL, 'b'},
    {"graph", required_argument, NULL, 'g'},
    {"tsv", required_argument, NULL, OPT_TSV},
    {"barcode-counts", required_argument, NULL, OPT_BARCODE_COUNTS},
    {"gap", required_argument, NULL, OPT_GAP },
    {"index_multiplicity", required_argument, NULL, 'm'},
    {"max_degree", required_argument, NULL, 'd'},
    {"end_length", required_argument, NULL, 'e'},
    {"error_percent", required_argument, NULL, 'r'},
    {"run_verbose", required_argument, NULL, 'v'},
    {"version", no_argument, NULL, OPT_VERSION},
    {"help", no_argument, NULL, OPT_HELP},
    { NULL, 0, NULL, 0 }
};

/** Command line parameters. */
static ARCS::ArcsParams params;

// Declared in Graph/Options.h and used by DotIO
namespace opt {
    /** The size of a k-mer. */
    unsigned k;

    /** The file format of the graph when writing. */
    int format;
}

/** A distance estimate graph. */
typedef ContigGraph<DirectedGraph<Length, DistanceEst>> DistGraph;

/** A dictionary of contig names. */
Dictionary g_contigNames;
unsigned g_nextContigName;

/**
 * One end of a scaffold.
 * The left (head) is true and the right (tail) is false.
 */
typedef std::pair<std::string, bool> ScaffoldEnd;

/** Hash a ScaffoldEnd. */
struct HashScaffoldEnd {
    size_t operator()(const ScaffoldEnd& key) const {
        return std::hash<std::string>()(key.first) ^ key.second;
    }
};

/* Returns true if seqence only contains ATGC and is of length indexLen */
bool checkIndex(std::string seq) {
    for (int i = 0; i < static_cast<int>(seq.length()); i++) {
        char c = toupper(seq[i]);
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C')
            return false;
    }
    //return (static_cast<int>(seq.length()) == params.indexLen);
    return true;
}

/*
 * Check if SAM flag is one of the accepted ones.
 */
bool checkFlag(int flag) {
    return (flag == 99 || flag == 163 || flag == 83 || flag == 147);
}

/*
 * Check if character is one of the accepted ones.
 */
bool checkChar(char c) {
    return (c == 'M' || c == '=' || c == 'X' || c == 'I');
}

/*
 * Calculate the sequence identity from the cigar string
 * sequence length, and tags.
 */
double calcSequenceIdentity(const std::string& line, const std::string& cigar, const std::string& seq) {

    int qalen = 0;
    std::stringstream ss;
    for (auto i = cigar.begin(); i != cigar.end(); ++i) {
        if (!isdigit(*i)) {
            if (checkChar(*i)) {
                ss << "\t";
                int value = 0;
                ss >> value;
                qalen += value;
                ss.str("");
            } else {
                ss.str("");
            }
        } else {
            ss << *i;
        }
    }

    int edit_dist = 0;
    std::size_t found = line.find("NM:i:");
    if (found!=std::string::npos) {
        edit_dist = std::strtol(&line[found + 5], 0, 10);
    }

    double si = 0;
    if (qalen != 0) {
        double mins = qalen - edit_dist;
        double div = mins/seq.length();
        si = div * 100;
    }

    return si;
}


/* Get all scaffold sizes from FASTA file */
void getScaffSizes(std::string file, ARCS::ScaffSizeList& scaffSizes) {

    int counter = 0;
    FastaReader in(file.c_str(), FastaReader::FOLD_CASE);
    for (FastaRecord rec; in >> rec;) {
        counter++;
        std::string  scafName = rec.id;
        int size = rec.seq.length();
        scaffSizes.push_back(std::make_pair(rec.id, size));
    }

    if (params.verbose)
        std::cout << "Saw " << counter << " sequences.\n";
}

/*
 * Read BAM file, if sequence identity greater than threashold
 * update indexMap. IndexMap also stores information about
 * contig number index algins with and counts.
 */
void readBAM(const std::string bamName, ARCS::IndexMap& imap,
    std::unordered_map<std::string, int>& indexMultMap,
    BarcodeToIndexMap& barcodeToIndex,
    SegmentToBarcode& segmentToBarcode,
    std::unordered_map<std::string, int> sMap)
{

    /* Open BAM file */
    std::ifstream bamName_stream;
    bamName_stream.open(bamName.c_str());
    assert_good(bamName_stream, bamName);

    std::string prevRN = "", readyToAddIndex = "", prevRef = "", readyToAddRefName = "";
    int prevSI = 0, prevFlag = 0, prevMapq = 0, prevPos = -1,
        prevTSpan = -1, readyToAddPos = -1, readyToAddTStart = -1,
        readyToAddTEnd = -1;
    int ct = 1;

    std::string line;
    size_t linecount = 0;

    // Number of unpaired reads.
    size_t countUnpaired = 0;

    /* Read each line of the BAM file */
    while (getline(bamName_stream, line)) {

        /* Check to make sure it is not the header */
        if (line.substr(0, 1).compare("@") != 0) {
            linecount++;

            std::stringstream ss(line);
            std::string readName, scafName, cigar, rnext, seq, qual, tags;
            int flag, pos, mapq, pnext, tlen;

            ss >> readName >> flag >> scafName >> pos >> mapq >> cigar
                >> rnext >> pnext >> tlen >> seq >> qual >> std::ws;

            getline(ss, tags);

            /* calculate req/query alignment coords from CIGAR string */
            SAMAlignment::CigarCoord cigarCoords(cigar);

            /* Parse the index from the readName */
            std::string index;
            if (params.bx) {
                index = getBarcodeSeq(tags);
            } else {
                std::size_t found = readName.rfind("_");
                if (found!=std::string::npos)
                    index = readName.substr(found+1);
            }

            if (!checkIndex(index))
                index.clear();

            /* Keep track of index multiplicity */
            if (!index.empty())
                indexMultMap[index]++;

            /* Calculate the sequence identity */
            int si = calcSequenceIdentity(line, cigar, seq);

            if (ct == 2 && readName != prevRN) {
                if (countUnpaired == 0)
                    std::cerr << "Warning: Skipping an unpaired read. BAM file should be sorted in order of read name.\n"
                        "  Prev read: " << prevRN << "\n"
                        "  Curr read: " << readName << std::endl;
                ++countUnpaired;
                if (countUnpaired % 1000000 == 0)
                    std::cerr << "Warning: Skipped " << countUnpaired << " unpaired reads." << std::endl;
                ct = 1;
            }

            if (ct >= 3)
                ct = 1;
            if (ct == 1) {
                if (readName.compare(prevRN) != 0) {
                    prevRN = readName;
                    prevSI = si;
                    prevFlag = flag;
                    prevMapq = mapq;
                    prevRef = scafName;
                    prevPos = pos;
                    prevTSpan = cigarCoords.tspan;

                    /*
                     * Read names are different so we can add the previous index and scafName as
                     * long as there were only two mappings (one for each read)
                     */
                    if (!readyToAddIndex.empty() && !readyToAddRefName.empty() && readyToAddRefName.compare("*") != 0 && readyToAddPos != -1) {

                        int size = sMap[readyToAddRefName];
                        if (size >= params.min_size) {

                           /*
                            * If length of sequence is less than 2 x end_length, split
                            * the sequence in half to determing head/tail
                            */
                           int cutOff = params.end_length;
                           if (cutOff == 0 || size <= cutOff * 2)
                               cutOff = size/2;

                           /*
                            * pair <X, true> indicates read pair aligns to head,
                            * pair <X, false> indicates read pair aligns to tail
                            */
                           ScaffoldEnd key(readyToAddRefName, true);
                           ScaffoldEnd keyR(readyToAddRefName, false);

                           /* Aligns to head */
                           if (readyToAddPos <= cutOff) {
                               imap[readyToAddIndex][key]++;

                               if (imap[readyToAddIndex].count(keyR) == 0)
                                   imap[readyToAddIndex][keyR] = 0;

                            /* Aligns to tail */
                           } else if (readyToAddPos > size - cutOff) {
                               imap[readyToAddIndex][keyR]++;

                               if (imap[readyToAddIndex].count(key) == 0)
                                   imap[readyToAddIndex][key] = 0;
                           }

                           BarcodeIndex barcodeIndex =
                               barcodeToIndex.getIndex(readyToAddIndex);

                           assert(readyToAddTStart > 0);
                           assert(readyToAddTEnd > 0);

                           SegmentCalc calc(params.segment_length);
                           std::pair<SegmentIndex, SegmentIndex> range;
                           bool valid;
                           boost::tie(range, valid) = calc.indexRange(
                               readyToAddTStart, readyToAddTEnd, size);

                           if (valid) {
                               for (unsigned i = range.first;
                                    i <= range.second; ++i) {
                                   Segment segment(readyToAddRefName, i);
                                   segmentToBarcode[segment][barcodeIndex]++;
                               }
                           }

                        }

                        readyToAddIndex = "";
                        readyToAddRefName = "";
                        readyToAddPos = -1;
                    }
                } else {
                    ct = 0;
                    readyToAddIndex = "";
                    readyToAddRefName = "";
                    readyToAddPos = -1;
                }
            } else if (ct == 2) {
                assert(readName == prevRN);
                if (!seq.empty() && checkFlag(flag) && checkFlag(prevFlag)
                        && mapq != 0 && prevMapq != 0 && si >= params.seq_id && prevSI >= params.seq_id) {
                    if (prevRef.compare(scafName) == 0 && scafName.compare("*") != 0 && !scafName.empty() && !index.empty()) {

                        readyToAddIndex = index;
                        readyToAddRefName = scafName;
                        /* Take average read alignment position between read pairs */
                        readyToAddPos = (prevPos + pos)/2;

                        assert(cigarCoords.tspan > 0);
                        assert(prevTSpan > 0);
                        readyToAddTStart = std::min(pos, prevPos);
                        readyToAddTEnd = std::max(pos + cigarCoords.tspan - 1,
                            unsigned(prevPos + prevTSpan - 1));
                    }
                }
            }
           ct++;

            if (params.verbose && linecount % 10000000 == 0)
                std::cout << "On line " << linecount << std::endl;

        }
        assert(bamName_stream);
    }

    /* Close BAM file */
    assert_eof(bamName_stream, bamName);
    bamName_stream.close();

    if (countUnpaired > 0)
        std::cerr << "Warning: Skipped " << countUnpaired << " unpaired reads. BAM file should be sorted in order of read name.\n";
}

/**
 * Read the file of file names.
 */
static std::vector<std::string> readFof(const std::string& fofName)
{
    std::vector<std::string> filenames;
    if (fofName.empty())
      return filenames;
    std::ifstream fin(fofName);
    assert_good(fin, fofName);
    for (std::string bamName; getline(fin, bamName);)
        filenames.push_back(bamName);
    assert_eof(fin, fofName);
    if (filenames.empty()) {
        cerr << PROGRAM ": error: " << fofName << " is empty" << std::endl;
        exit(EXIT_FAILURE);
    }
    for (const auto& filename : filenames)
      assert_readable(filename);
    return filenames;
}

/**
 * Read the BAM files.
 */
void readBAMS(const std::vector<std::string> bamNames, ARCS::IndexMap& imap,
    std::unordered_map<std::string, int>& indexMultMap,
    BarcodeToIndexMap& barcodeToIndex,
    SegmentToBarcode& segmentToBarcode,
    std::unordered_map<std::string, int> sMap) {
    assert(!bamNames.empty());
    for (const auto& bamName : bamNames) {
        if (params.verbose)
            std::cout << "Reading bam " << bamName << std::endl;
        readBAM(bamName, imap, indexMultMap,
            barcodeToIndex, segmentToBarcode, sMap);
    }
}

/** Count barcodes. */
static size_t countBarcodes(ARCS::IndexMap& imap, const std::unordered_map<std::string, int>& indexMultMap)
{
    size_t barcodeCount = 0;
    for (auto x : indexMultMap)
        if (x.second >= params.min_mult && x.second <= params.max_mult)
            ++barcodeCount;

    std::cout
        << "{ \"All_barcodes_unfiltered\":" << indexMultMap.size()
        << ", \"All_barcodes_filtered\":" << barcodeCount
        << ", \"Scaffold_end_barcodes\":" << imap.size()
        << ", \"Min_barcode_reads_threshold\":" << params.min_mult
        << ", \"Max_barcode_reads_threshold\":" << params.max_mult
        << " }\n";

    return barcodeCount;
}

/* Normal approximation to the binomial distribution */
float normalEstimation(int x, float p, int n) {
    float mean = n * p;
    float sd = std::sqrt(n * p * (1 - p));
    return 0.5 * (1 + std::erf((x - mean)/(sd * std::sqrt(2))));
}

/*
 * Based on number of read pairs that align to the
 * head or tail of scaffold, determine if is significantly
 * different from a uniform distribution (p=0.5)
 */
std::pair<bool, bool> headOrTail(int head, int tail) {
    int max = std::max(head, tail);
    int sum = head + tail;
    if (sum < params.min_reads) {
        return std::pair<bool, bool> (false, false);
    }
    float normalCdf = normalEstimation(max, 0.5, sum);
    if (1 - normalCdf < params.error_percent) {
        bool isHead = (max == head);
        return std::pair<bool, bool> (true, isHead);
    } else {
        return std::pair<bool, bool> (false, false);
    }
}

/*
 * Iterate through IndexMap and for every pair of scaffolds
 * that align to the same index, store in PairMap. PairMap
 * is a map with a key of pairs of saffold names, and value
 * of number of links between the pair. (Each link is one index).
 */
void pairContigs(ARCS::IndexMap& imap, ARCS::PairMap& pmap, std::unordered_map<std::string, int>& indexMultMap) {

    /* Iterate through each index in IndexMap */
    for(auto it = imap.begin(); it != imap.end(); ++it) {

        /* Get index multiplicity from indexMultMap */
        std::string index = it->first;
        int indexMult = indexMultMap[index];

        if (indexMult >= params.min_mult && indexMult <= params.max_mult) {

           /* Iterate through all the scafNames in ScafMap */
            for (auto o = it->second.begin(); o != it->second.end(); ++o) {
                for (auto p = it->second.begin(); p != it->second.end(); ++p) {
                    std::string scafA, scafB;
                    bool scafAflag, scafBflag;
                    std::tie (scafA, scafAflag) = o->first;
                    std::tie (scafB, scafBflag) = p->first;

                    /* Only insert into pmap if scafA < scafB to avoid duplicates */
                    if (scafA < scafB && scafAflag && scafBflag) {
                        bool validA, validB, scafAhead, scafBhead;

                        std::tie(validA, scafAhead) = headOrTail(it->second[ScaffoldEnd(scafA, true)], it->second[ScaffoldEnd(scafA, false)]);
                        std::tie(validB, scafBhead) = headOrTail(it->second[ScaffoldEnd(scafB, true)], it->second[ScaffoldEnd(scafB, false)]);

                        if (validA && validB) {
                            std::pair<std::string, std::string> pair (scafA, scafB);
                            if (pmap.count(pair) == 0)
                                pmap[pair].resize(4);
                            // Head - Head
                            if (scafAhead && scafBhead)
                                pmap[pair][0]++;
                            // Head - Tail
                            else if (scafAhead && !scafBhead)
                                pmap[pair][1]++;
                            // Tail - Head
                            else if (!scafAhead && scafBhead)
                                pmap[pair][2]++;
                            // Tail - Tail
                            else if (!scafAhead && !scafBhead)
                                pmap[pair][3]++;
                        }
                    }
                }
            }
        }
    }
}

/*
 * Return the max value and its index position
 * in the vector
 */
std::pair<unsigned, unsigned> getMaxValueAndIndex(const ARCS::PairMap::value_type::second_type& array) {
    unsigned max = 0;
    unsigned index = 0;
    for (unsigned i = 0; i < array.size(); ++i) {
        if (array[i] > max) {
            max = array[i];
            index = i;
        }
    }
    return std::make_pair(max, index);
}

/*
 * Return true if the link orientation with the max support
 * is dominant
 */
bool checkSignificance(int max, int second) {
    if (max < params.min_links) {
        return false;
    }
    float normalCdf = normalEstimation(max, 0.5, second);
    return (1 - normalCdf < params.error_percent);
}

/*
 * Construct a boost graph from PairMap. Each pair represents an
 * edge in the graph. The weight of each edge is the number of links
 * between the scafNames.
 * VidVdes is a mapping of vertex descriptors to scafNames (vertex id).
 */
void createGraph(const ARCS::PairMap& pmap, ARCS::Graph& g) {

    ARCS::VidVdesMap vmap;

    ARCS::PairMap::const_iterator it;
    for(it = pmap.begin(); it != pmap.end(); ++it) {
        std::string scaf1, scaf2;
        std::tie (scaf1, scaf2) = it->first;

        unsigned max, index;
        const auto& count = it->second;
        std::tie(max, index) = getMaxValueAndIndex(count);

        unsigned second = 0;
        for (unsigned i = 0; i < count.size(); ++i) {
            if (count[i] != max && count[i] > second)
                second = count[i];
        }

        /* Only insert edge if orientation with max links is dominant */
        if (checkSignificance(max, second)) {

            /* If scaf1 is not a node in the graph, add it */
            if (vmap.count(scaf1) == 0) {
                ARCS::Graph::vertex_descriptor v = boost::add_vertex(g);
                g[v].id = scaf1;
                vmap[scaf1] = v;
            }

            /* If scaf2 is not a node in the graph, add it */
            if (vmap.count(scaf2) == 0) {
                ARCS::Graph::vertex_descriptor v = boost::add_vertex(g);
                g[v].id = scaf2;
                vmap[scaf2] = v;
            }

            ARCS::Graph::edge_descriptor e;
            bool inserted;

            /* Add the edge representing the pair */
            std::tie (e, inserted) = boost::add_edge(vmap[scaf1], vmap[scaf2], g);
            if (inserted) {
                g[e].weight = max;
                g[e].orientation = index;
            }
        }
    }
}

/*
 * Write out the boost graph in a .dot file.
 */
void writeGraph(const std::string& graphFile_dot, ARCS::Graph& g)
{
	std::ofstream out(graphFile_dot.c_str());
	assert(out);

	ARCS::VertexPropertyWriter<ARCS::Graph> vpWriter(g);
	ARCS::EdgePropertyWriter<ARCS::Graph> epWriter(g);

	boost::write_graphviz(out, g, vpWriter, epWriter);
	assert(out);
	out.close();
}

/*
 * Remove all nodes from graph wich have a degree
 * greater than max_degree
 */
void removeDegreeNodes(ARCS::Graph& g, int max_degree) {

    boost::graph_traits<ARCS::Graph>::vertex_iterator vi, vi_end, next;
    boost::tie(vi, vi_end) = boost::vertices(g);

    std::vector<ARCS::VertexDes> dVertex;
    for (next = vi; vi != vi_end; vi = next) {
        ++next;
        if (static_cast<int>(boost::degree(*vi, g)) > max_degree) {
            dVertex.push_back(*vi);
        }
    }

    for (unsigned i = 0; i < dVertex.size(); i++) {
        boost::clear_vertex(dVertex[i], g);
        boost::remove_vertex(dVertex[i], g);
    }
    boost::renumber_indices(g);
}

/*
 * Remove nodes that have a degree greater than max_degree
 * Write graph
 */
void writePostRemovalGraph(ARCS::Graph& g, const std::string graphFile) {
    if (params.max_degree != 0) {
        std::cout << "      Deleting nodes with degree > " << params.max_degree <<"... \n";
        removeDegreeNodes(g, params.max_degree);
    } else {
        std::cout << "      Max Degree (-d) set to: " << params.max_degree << ". Will not delete any vertices from graph.\n";
    }

    std::cout << "      Writing graph file to " << graphFile << "...\n";
    writeGraph(graphFile, g);
}

/*
 * Construct an ABySS distance estimate graph from a boost graph.
 */
void createAbyssGraph(const ARCS::ScaffSizeList& scaffSizes, const ARCS::Graph& gin, DistGraph& gout) {
    // Add the vertices.
    for (const auto& it : scaffSizes) {
        vertex_property<DistGraph>::type vp;
        vp.length = it.second;
        const auto u = add_vertex(vp, gout);
        put(vertex_name, gout, u, it.first);
    }

    // Add the edges.
    for (const auto ein : boost::make_iterator_range(boost::edges(gin))) {
        const auto einp = gin[ein];
        const auto u = find_vertex(gin[source(ein, gin)].id, einp.orientation < 2, gout);
        const auto v = find_vertex(gin[target(ein, gin)].id, einp.orientation % 2, gout);

        edge_property<DistGraph>::type ep;
        ep.distance = params.gap;
        ep.stdDev = params.gap;
        ep.numPairs = einp.weight;

        /* use distance estimates, if enabled */
        if (params.dist_est) {
            if (params.dist_mode == ARCS::DIST_MEDIAN) {
                ep.distance = einp.dist;
            } else {
                assert(params.dist_mode == ARCS::DIST_UPPER);
                ep.distance = einp.maxDist;
            }
        }

        graph_traits<DistGraph>::edge_descriptor e;
        bool inserted;
        std::tie(e, inserted) = add_edge(u, v, ep, gout);
        if (!inserted) {
            std::cerr << "error: Duplicate edge: \"" <<
                get(vertex_name, gout, u) << "\" -> \"" <<
                get(vertex_name, gout, v) << '"' << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

/*
 * Write out an ABySS distance estimate graph to a .dist.gv file.
 */
void writeAbyssGraph(const std::string& path, const DistGraph& g) {
    if (path.empty())
        return;
    ofstream out(path.c_str());
    assert_good(out, path);
    write_dot(out, g, "arcs");
    assert_good(out, path);
}

/** Write a TSV file of the number of reads per barcode.
 * - Barcode: the barcode
 * - Reads: the number of reads
 */
void writeBarcodeCountsTSV(
        const std::string& tsvFile,
        const std::unordered_map<std::string, int>& indexMultMap)
{
    if (tsvFile.empty())
        return;

    // Sort the barcodes by their counts and then their sequence.
    typedef std::vector<std::pair<std::string, unsigned>> Sorted;
    Sorted sorted(indexMultMap.begin(), indexMultMap.end());
    sort(sorted.begin(), sorted.end(),
            [](const Sorted::value_type& a, const Sorted::value_type& b) {
                return a.second != b.second ? a.second > b.second : a.first < b.first;
            });

    std::ofstream f(tsvFile);
    assert_good(f, tsvFile);
    f << "Barcode\tReads\n";
    assert_good(f, tsvFile);
    for (auto x : sorted)
        f << x.first << '\t' << x.second << '\n';
    assert_good(f, tsvFile);
}

/** Write a TSV file to calculate a hypergeometric test.
 * For each pair of scaffolds...
 * - CountBoth: the number of barcodes shared between scaffold ends U and V
 * - CountU: the number of barcodes on scaffold end U
 * - CountV: the number of barcodes on scaffold end V
 * - CountAll: the total number of barcodes observed
 */
void writeTSV(
        const std::string& tsvFile,
        const ARCS::IndexMap& imap,
        const ARCS::PairMap& pmap,
        size_t barcodeCount)
{
    if (tsvFile.empty())
        return;

    // Count the number of barcodes seen per scaffold end.
    std::unordered_map<ScaffoldEnd, unsigned, HashScaffoldEnd> barcodes_per_scaffold_end;
    for (const auto& it : imap) {
        for (const auto& scaffold_count : it.second) {
            const auto& scaffold = scaffold_count.first;
            const auto& count = scaffold_count.second;
            if (count >= params.min_reads)
               ++barcodes_per_scaffold_end[scaffold];
        }
    }

    std::ofstream f(tsvFile);
    assert_good(f, tsvFile);
    f << "U\tV\tBest_orientation\tShared_barcodes\tU_barcodes\tV_barcodes\tAll_barcodes\n";
    assert_good(f, tsvFile);
    for (const auto& it : pmap) {
        const auto& u = it.first.first;
        const auto& v = it.first.second;
        const auto& counts = it.second;
        assert(!counts.empty());
        unsigned max_counts = *std::max_element(counts.begin(), counts.end());
        for (unsigned i = 0; i < counts.size(); ++i) {
            if (counts[i] == 0)
                continue;
            bool usense = i < 2;
            bool vsense = i % 2;
            f << u << (usense ? '-' : '+')
                << '\t' << v << (vsense ? '-' : '+')
                << '\t' << (counts[i] == max_counts ? "T" : "F")
                << '\t' << counts[i]
                << '\t' << barcodes_per_scaffold_end[std::make_pair(u, usense)]
                << '\t' << barcodes_per_scaffold_end[std::make_pair(v, !vsense)]
                << '\t' << barcodeCount
                << '\n';
            f << v << (vsense ? '+' : '-')
                << '\t' << u << (usense ? '+' : '-')
                << '\t' << (counts[i] == max_counts ? "T" : "F")
                << '\t' << counts[i]
                << '\t' << barcodes_per_scaffold_end[std::make_pair(v, !vsense)]
                << '\t' << barcodes_per_scaffold_end[std::make_pair(u, usense)]
                << '\t' << barcodeCount
                << '\n';
        }
    }
    assert_good(f, tsvFile);
}

/** Return NA if the specified string is empty, and the string itself otherwise. */
static const char* maybeNA(const std::string& s)
{
    return s.empty() ? "NA" : s.c_str();
}

/**
 * calculate distance estimates for edges, based on number of shared
 * barcodes between contig ends
 */
static inline void calcDistanceEstimates(
    const ARCS::IndexMap& imap,
    const std::unordered_map<std::string, int> &indexMultMap,
    const ARCS::ContigToLength& contigToLength,
    ARCS::Graph& g)
{
    std::time_t rawtime;

    time(&rawtime);
    std::cout << "\n\t=>Measuring intra-contig distances / shared barcodes... "
        << ctime(&rawtime);
    DistSampleMap distSamples;
    calcDistSamples(imap, contigToLength, indexMultMap, params, distSamples);

    time(&rawtime);
    std::cout << "\n\t=>Writing intra-contig distance samples to TSV... "
        << ctime(&rawtime);
    writeDistSamplesTSV(params.dist_samples_tsv, distSamples);

    time(&rawtime);
    std::cout << "\n\t=>Building Jaccard => distance map... "
        << ctime(&rawtime);
    JaccardToDist jaccardToDist;
    buildJaccardToDist(distSamples, jaccardToDist);

    time(&rawtime);
    std::cout << "\n\t=>Calculating barcode stats for scaffold pairs... "
        << ctime(&rawtime);
    PairToBarcodeStats pairToStats;
    buildPairToBarcodeStats(imap, indexMultMap, contigToLength, params, pairToStats);

    time(&rawtime);
    std::cout << "\n\t=>Adding edge distances... " << ctime(&rawtime);
    addEdgeDistances(pairToStats, jaccardToDist, params, g);

    time(&rawtime);
    std::cout << "\n\t=>Writing distance estimates to TSV... "
        << ctime(&rawtime);
    writeDistTSV(params.dist_tsv, pairToStats, g);
}

/** Run ARCS. */
void runArcs(const std::vector<std::string>& filenames) {

    std::cout << "Running: " << PROGRAM << " " << PACKAGE_VERSION
        << "\n pid " << ::getpid()
        // Options
        << "\n -c " << params.min_reads
        << "\n -d " << params.max_degree
        << "\n -e " << params.end_length
        << "\n -l " << params.min_links
        << "\n -m " << params.min_mult << '-' << params.max_mult
        << "\n -r " << params.error_percent
        << "\n -s " << params.seq_id
        << "\n -v " << params.verbose
        << "\n -z " << params.min_size
        << "\n --gap=" << params.gap
        // Output files
        << "\n -b " << maybeNA(params.base_name)
        << "\n -g " << maybeNA(params.dist_graph_name)
        << "\n --barcode-counts=" << maybeNA(params.barcode_counts_name)
        << "\n --tsv=" << maybeNA(params.tsv_name)
        // Input files
        << "\n -a " << maybeNA(params.fofName)
        << "\n -f " << maybeNA(params.file)
        << '\n';
    for (const auto& filename : filenames)
        std::cout << ' ' << filename << '\n';
    std::cout.flush();

    std::string graphFile = params.base_name + "_original.gv";
    std::string distGraphFile = !params.dist_graph_name.empty() ? params.dist_graph_name : params.base_name + ".dist.gv";

    ARCS::IndexMap imap;
    ARCS::PairMap pmap;
    ARCS::Graph g;

    std::time_t rawtime;

    ARCS::ScaffSizeList scaffSizeList;
    time(&rawtime);
    std::cout << "\n=>Getting scaffold sizes... " << ctime(&rawtime);
    getScaffSizes(params.file, scaffSizeList);
    ARCS::ScaffSizeMap scaffSizeMap(
        scaffSizeList.begin(), scaffSizeList.end());

    /* maps Chromium barcodes to integer indices */
    BarcodeToIndexMap barcodeToIndex;
    /* maps contig segments to barcodes */
    SegmentToBarcode segmentToBarcode;
    /* maps Chromium barcode sequence to number of read pair mappings */
    std::unordered_map<std::string, int> indexMultMap;
    time(&rawtime);
    std::cout << "\n=>Starting to read BAM files... " << ctime(&rawtime);
    std::vector<std::string> bamFiles = readFof(params.fofName);
    std::copy(filenames.begin(), filenames.end(), std::back_inserter(bamFiles));
    readBAMS(bamFiles, imap, indexMultMap, barcodeToIndex,
        segmentToBarcode, scaffSizeMap);

    size_t barcodeCount = countBarcodes(imap, indexMultMap);

    time(&rawtime);
    std::cout << "\n=>Starting to write reads per barcode TSV file... " << ctime(&rawtime) << "\n";
    writeBarcodeCountsTSV(params.barcode_counts_name, indexMultMap);

    time(&rawtime);
    std::cout << "\n=>Starting pairing of scaffolds... " << ctime(&rawtime);
    pairContigs(imap, pmap, indexMultMap);

    time(&rawtime);
    std::cout << "\n=>Starting to create graph... " << ctime(&rawtime);
    createGraph(pmap, g);

    if (params.dist_est) {
        std::cout << "\n=>Calculating distance estimates... " << ctime(&rawtime);
        calcDistanceEstimates(imap, indexMultMap, scaffSizeMap, g);
    }

    time(&rawtime);
    std::cout << "\n=>Starting to write graph file... " << ctime(&rawtime) << "\n";
    writePostRemovalGraph(g, graphFile);

    time(&rawtime);
    std::cout << "\n=>Starting to create ABySS graph... " << ctime(&rawtime);
    DistGraph gdist;
    createAbyssGraph(scaffSizeList, g, gdist);

    time(&rawtime);
    std::cout << "\n=>Starting to write ABySS graph file... " << ctime(&rawtime) << "\n";
    writeAbyssGraph(distGraphFile, gdist);

    time(&rawtime);
    std::cout << "\n=>Starting to write TSV file... " << ctime(&rawtime) << "\n";
    writeTSV(params.tsv_name, imap, pmap, barcodeCount);

    time(&rawtime);
    std::cout << "\n=>Done. " << ctime(&rawtime);
}

int main(int argc, char** argv) {

    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case '?':
                die = true; break;
            case 'f':
                arg >> params.file; break;
            case 'a':
                arg >> params.fofName; break;
            case 'B':
                arg >> params.dist_bin_size; break;
            case 's':
                arg >> params.seq_id; break;
            case 'c':
                arg >> params.min_reads; break;
            case 'D':
                params.dist_est = true; break;
            case 'l':
                arg >> params.min_links; break;
            case 'z':
                arg >> params.min_size; break;
            case 'b':
                arg >> params.base_name; break;
            case 'g':
                arg >> params.dist_graph_name; break;
            case OPT_BX:
                params.bx = true; break;
            case OPT_TSV:
                arg >> params.tsv_name; break;
            case OPT_GAP:
                arg >> params.gap; break;
            case OPT_BARCODE_COUNTS:
                arg >> params.barcode_counts_name; break;
            case OPT_SAMPLES_TSV:
                arg >> params.dist_samples_tsv; break;
            case OPT_DIST_TSV:
                arg >> params.dist_tsv; break;
            case OPT_NO_DIST_EST:
                params.dist_est = false; break;
            case OPT_DIST_MEDIAN:
                params.dist_mode = ARCS::DIST_MEDIAN; break;
            case OPT_DIST_UPPER:
                params.dist_mode = ARCS::DIST_UPPER; break;
            case 'm': {
                std::string firstStr, secondStr;
                std::getline(arg, firstStr, '-');
                std::getline(arg, secondStr);
                std::stringstream ss;
                ss << firstStr << "\t" << secondStr;
                ss >> params.min_mult >> params.max_mult;
                }
                break;
            case 'd':
                arg >> params.max_degree; break;
            case 'e':
                arg >> params.end_length; break;
            case 'r':
                arg >> params.error_percent; break;
            case 'v':
                ++params.verbose; break;
            case OPT_HELP:
                std::cout << USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
        if (optarg != NULL && (!arg.eof() || arg.fail())) {
            std::cerr << PROGRAM ": invalid option: `-"
                << (char)c << optarg << "'\n";
            exit(EXIT_FAILURE);
        }
    }

    if (params.file.empty()) {
        cerr << PROGRAM ": missing -f option\n";
        die = true;
    }

    std::vector<std::string> filenames(argv + optind, argv + argc);
    if (params.fofName.empty() && filenames.empty()) {
        cerr << PROGRAM ": missing either BAM files or -a option\n";
        die = true;
    }

    if (die) {
        std::cerr << "Try " PROGRAM " --help for more information.\n";
        exit(EXIT_FAILURE);
    }

    assert_readable(params.file);
    if (!params.fofName.empty())
      assert_readable(params.fofName);
    for (const auto& filename : filenames)
      assert_readable(filename);

    /* Setting base name if not previously set */
    if (params.base_name.empty()) {
        std::ostringstream filename;
        filename << params.file << ".scaff"
            << "_s" << params.seq_id
            << "_c" << params.min_reads
            << "_l" << params.min_links
            << "_d" << params.max_degree
            << "_e" << params.end_length
            << "_r" << params.error_percent;
        params.base_name = filename.str();
    }

    runArcs(filenames);

    return 0;
}
