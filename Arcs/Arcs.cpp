#include "Arcs.h"

#define PROGRAM "arcs"
#define VERSION "1.3.1"

static const char VERSION_MESSAGE[] = 
"VERSION: " PROGRAM " " VERSION "\n"
"\n"
"RAILS.readme distributed with this software @ www.bcgsc.ca \n"
"http://www.bcgsc.ca/platform/bioinfo/software/links \n"
"We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca.\n"
"If you use LINKS, RAILS, ARCS code or ideas, please cite our work. \n"
"\n"
"LINKS, RAILS and ARCS Copyright (c) 2014-2016 Canada's Michael Smith Genome Science Centre.  All rights reserved. \n";

static const char USAGE_MESSAGE[] =
"Usage: [" PROGRAM " " VERSION "]\n"
"   -f  Assembled Sequences to further scaffold (Multi-Fasta format, required)\n"
"   -a  File of File Names listing all input BAM alignment files (required). \n"
"       NOTE: alignments must be sorted in order of name\n"
"             index must be included in read name e.g read1_indexA\n"
"   -s  Minimum sequence identity (min. required to include the read's scaffold alignment in the graph file, default: 98)\n"
"   -c  Minimum number of mapping read pairs/Index required before creating edge in graph. (default: 5)\n"
"   -l  Minimum number of links to create edge in graph (default: 0)\n"
"   -z  Minimum contig length to consider for scaffolding (default: 500)\n"
"   -b  Base name for your output files (optional)\n"
"   -m  Range (in the format min-max) of index multiplicity (only reads with indices in this multiplicity range will be included in graph) (default: 50-1000)\n"
"   -d  Maximum degree of nodes in graph. All nodes with degree greater than this number will be removed from the graph prior to printing final graph. For no node removal, set to 0 (default: 0)\n"
"   -e  Length (bp) of ends of read to consider (default: 30000)\n"
"   -r  Maximum P-value for link orientation determination. Lower is more stringent (default: 0.05)\n"
"   -i  Length (bp) of index sequence. 14 for GemCode, 16 for Chromium (default: 16)\n"
"   -v  Runs in verbose mode (optional, default: 0)\n";


ARCS::ArcsParams params;

static const char shortopts[] = "f:a:s:c:l:z:b:m:d:e:r:i:v";

enum { OPT_HELP = 1, OPT_VERSION};

static const struct option longopts[] = {
    {"file", required_argument, NULL, 'f'},
    {"fofName", required_argument, NULL, 'a'},
    {"seq_id", required_argument, NULL, 's'}, 
    {"min_reads", required_argument, NULL, 'c'},
    {"min_links", required_argument, NULL, 'l'},
    {"min_size", required_argument, NULL, 'z'},
    {"base_name", required_argument, NULL, 'b'},
    {"index_multiplicity", required_argument, NULL, 'm'},
    {"max_degree", required_argument, NULL, 'd'},
    {"end_length", required_argument, NULL, 'e'},
    {"error_percent", required_argument, NULL, 'r'},
    {"index_length", required_argument, NULL, 'i'},
    {"run_verbose", required_argument, NULL, 'v'},
    {"version", no_argument, NULL, OPT_VERSION},
    {"help", no_argument, NULL, OPT_HELP},
    { NULL, 0, NULL, 0 }
};



/* Returns true if seqence only contains ATGC and is of length indexLen */
bool checkIndex(std::string seq) {
    for (int i = 0; i < static_cast<int>(seq.length()); i++) {
        char c = toupper(seq[i]);
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C')
            return false;
    }
    return (static_cast<int>(seq.length()) == params.indexLen);
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

/* Extract the integer scaffold id from scaffold name */
int getIntFromScafName(const std::string& scafName) {
    std::stringstream ss;
    for (auto i = scafName.begin(); i != scafName.end(); ++i) {
        if (isdigit(*i)) {
            ss << *i;
        }
    }

    int numb;
    ss >> numb;
    return numb;
}

/* Get all scaffold sizes from FASTA file */
void getScaffSizes(std::string file, std::unordered_map<int, int>& sMap) {

    int counter = 0;
    FastaReader in(file.c_str(), FastaReader::FOLD_CASE);
    for (FastaRecord rec; in >> rec;) {
        counter++;
        int scafName = getIntFromScafName(rec.id);
        int size = rec.seq.length();
        //assert(sMap.count(scafName) == 0);
        sMap[scafName] = size;
    }
    
    if (params.verbose)
        std::cout << "Saw " << counter << " sequences.\n";
}

/* 
 * Read BAM file, if sequence identity greater than threashold
 * update indexMap. IndexMap also stores information about
 * contig number index algins with and counts.
 */
void readBAM(const std::string bamName, ARCS::IndexMap& imap, std::unordered_map<std::string, int>& indexMultMap, std::unordered_map<int, int> sMap) {

    /* Open BAM file */
    std::ifstream bamName_stream;
    bamName_stream.open(bamName.c_str());
    if (!bamName_stream) {
        std::cerr << "Could not open " << bamName << ". --fatal.\n";
        exit(EXIT_FAILURE);
    }

    std::string prevRN = "", readyToAddIndex = "";
    int prevSI = 0, prevFlag = 0, prevMapq = 0, prevRef = 0, readyToAddRefName = 0, prevPos = -1, readyToAddPos = -1;
    int ct = 1; 

    std::string line;
    int linecount = 0;

    /* Read each line of the BAM file */
    while (getline(bamName_stream, line)) {
        
        /* Check to make sure it is not the header */
        if (line.substr(0, 1).compare("@") != 0) {
            linecount++;

            std::stringstream ss(line);
            std::string readName, scafNameString, cigar, rnext, seq, qual;
            int flag, scafName, pos, mapq, pnext, tlen;

            ss >> readName >> flag >> scafNameString >> pos >> mapq >> cigar
                >> rnext >> pnext >> tlen >> seq >> qual;

            /* If there are no numbers in name, will return 0 - ie "*" */
            scafName = getIntFromScafName(scafNameString);

            /* Parse the index from the readName */
            std::string index = "", readNameSpl = "";
            std::size_t found = readName.find("_");
            if (found!=std::string::npos)
                readNameSpl = readName.substr(found+1);
            if (checkIndex(readNameSpl))
                index = readNameSpl;

            /* Keep track of index multiplicity */
            if (!index.empty())
                indexMultMap[index]++;

            /* Calculate the sequence identity */
            int si = calcSequenceIdentity(line, cigar, seq);

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

                    /* 
                     * Read names are different so we can add the previous index and scafName as
                     * long as there were only two mappings (one for each read)
                     */
                    if (!readyToAddIndex.empty() && readyToAddRefName != 0 && readyToAddPos != -1) {

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
                           std::pair<int, bool> key(readyToAddRefName, true);
                           std::pair<int, bool> keyR(readyToAddRefName, false);

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

                        }
                        readyToAddIndex = "";
                        readyToAddRefName = 0;
                        readyToAddPos = -1;
                    }
                } else {
                    ct = 0;
                    readyToAddIndex = "";
                    readyToAddRefName = 0;
                    readyToAddPos = -1;
                }
            } else if (ct == 2) {
                if (prevRN.compare(readName) != 0) {
                    std::cerr << "ERROR! BAM file should be sorted in order of read name. Exiting... \n Prev Read: " << prevRN << "; Curr Read: " << readName << "\n";
                    exit(EXIT_FAILURE);
                }

                if (!seq.empty() && checkFlag(flag) && checkFlag(prevFlag)
                        && mapq != 0 && prevMapq != 0 && si >= params.seq_id && prevSI >= params.seq_id) {
                    if ((prevRef == scafName) && scafName != 0 && !index.empty()) {
                        readyToAddIndex = index;
                        readyToAddRefName = scafName;
                        /* Take average read alignment position between read pairs */
                        readyToAddPos = (prevPos + pos)/2;
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
    bamName_stream.close();
}

/* 
 * Reading each BAM file from fofName
 */
void readBAMS(const std::string& fofName, ARCS::IndexMap& imap, std::unordered_map<std::string, int>& indexMultMap, std::unordered_map<int, int> sMap) {

    std::ifstream fofName_stream(fofName.c_str());
    if (!fofName_stream) {
        std::cerr << "Could not open " << fofName << " ...\n";
        exit(EXIT_FAILURE);
    }

    std::string bamName;
    while (getline(fofName_stream, bamName)) {
        if (params.verbose)
            std::cout << "Reading bam " << bamName << std::endl;
        readBAM(bamName, imap, indexMultMap, sMap);
        assert(fofName_stream);
    }
    fofName_stream.close();
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
                    int scafA, scafB;
                    bool scafAflag, scafBflag;
                    std::tie (scafA, scafAflag) = o->first;
                    std::tie (scafB, scafBflag) = p->first;

                    /* Only insert into pmap if scafA < scafB to avoid duplicates */
                    if (scafA < scafB && scafAflag && scafBflag) {
                        bool validA, validB, scafAhead, scafBhead;

                        std::tie(validA, scafAhead) = headOrTail(it->second[std::pair<int, bool>(scafA, true)], it->second[std::pair<int, bool>(scafA, false)]);
                        std::tie(validB, scafBhead) = headOrTail(it->second[std::pair<int, bool>(scafB, true)], it->second[std::pair<int, bool>(scafB, false)]);

                        if (validA && validB) {
                            std::pair<int, int> pair (scafA, scafB);
                            if (pmap.count(pair) == 0) {
                                std::vector<int> init(4,0); 
                                pmap[pair] = init;
                            }
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
std::pair<int, int> getMaxValueAndIndex(const std::vector<int> array) {
    int max = 0;
    int index = 0;
    for (int i = 0; i < int(array.size()); i++) {
        if (array[i] > max) {
            max = array[i];
            index = i;
        }
    }

    std::pair<int, int> result(max, index);
    return result;
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
        int scaf1, scaf2;
        std::tie (scaf1, scaf2) = it->first;

        int max, index;
        std::vector<int> count = it->second;
        std::tie(max, index) = getMaxValueAndIndex(count);

        int second = 0;
        for (int i = 0; i < int(count.size()); i++) {
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
void writeGraph(const std::string& graphFile_dot, ARCS::Graph& g) {
    std::ofstream out(graphFile_dot.c_str());
    assert(out);

    boost::dynamic_properties dp;
    dp.property("id", get(&ARCS::VertexProperties::id, g));
    dp.property("weight", get(&ARCS::EdgeProperties::weight, g));
    dp.property("label", get(&ARCS::EdgeProperties::orientation, g));
    dp.property("node_id", get(boost::vertex_index, g));
    boost::write_graphviz_dp(out, g, dp);
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
        std::cout << "      Max Degree (-d) set to: " << params.max_degree << ". Will not delete any verticies from graph.\n";
    }

    std::cout << "      Writting graph file to " << graphFile << "...\n";
    writeGraph(graphFile, g);
}


void runArcs() {

    std::cout << "Running: " << PROGRAM << " " << VERSION 
        << "\n pid " << ::getpid()
        << "\n -f " << params.file 
        << "\n -a " << params.fofName
        << "\n -s " << params.seq_id 
        << "\n -c " << params.min_reads 
        << "\n -l " << params.min_links     
        << "\n -z " << params.min_size
        << "\n -b " << params.base_name
        << "\n Min index multiplicity: " << params.min_mult 
        << "\n Max index multiplicity: " << params.max_mult 
        << "\n -d " << params.max_degree 
        << "\n -e " << params.end_length
        << "\n -r " << params.error_percent
        << "\n -i " << params.indexLen 
        << "\n -v " << params.verbose << "\n";

    std::string graphFile = params.base_name + "_original.gv";

    ARCS::IndexMap imap;
    ARCS::PairMap pmap;
    ARCS::Graph g;

    std::time_t rawtime;

    std::unordered_map<int, int> scaffSizeMap;
    time(&rawtime);
    std::cout << "\n=>Getting scaffold sizes... " << ctime(&rawtime);
    getScaffSizes(params.file, scaffSizeMap);

    std::unordered_map<std::string, int> indexMultMap;
    time(&rawtime);
    std::cout << "\n=>Starting to read BAM files... " << ctime(&rawtime);
    readBAMS(params.fofName, imap, indexMultMap, scaffSizeMap);

    time(&rawtime);
    std::cout << "\n=>Starting pairing of scaffolds... " << ctime(&rawtime);
    pairContigs(imap, pmap, indexMultMap);

    time(&rawtime);
    std::cout << "\n=>Starting to create graph... " << ctime(&rawtime);
    createGraph(pmap, g);

    time(&rawtime);
    std::cout << "\n=>Starting to write graph file... " << ctime(&rawtime) << "\n";
    writePostRemovalGraph(g, graphFile);

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
            case 's':
                arg >> params.seq_id; break;
            case 'c':
                arg >> params.min_reads; break;
            case 'l':
                arg >> params.min_links; break;
            case 'z':
                arg >> params.min_size; break;
            case 'b':
                arg >> params.base_name; break;
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
            case 'i':
                arg >> params.indexLen; break;
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

    std::ifstream f(params.fofName.c_str());
    if (!f.good()) {
        std::cerr << "Cannot find -a " << params.fofName << ". Exiting... \n";
        die = true;
    }
    std::ifstream g(params.file.c_str());
    if (!g.good()) {
        std::cerr << "Cannot find -f " << params.file << ". Exiting... \n";
        die = true;
    }

    if (die) {
        std::cerr << "Try " << PROGRAM << " --help for more information.\n";
        exit(EXIT_FAILURE);
    }

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


    runArcs();

    return 0;
}
