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
"   -a  File of File Names listing all input BAM alignment files (required). NOTE: alignments must be sorted in order of name\n"
"   -s  Minimum sequence identity (min. required to include the read's scaffold alignment in the graph file, default: 90)\n"
"   -c  Minimum number of mapping read pairs/Index required before creating edge in graph. (default: 2)\n"
"   -l  Minimum number of links to connect scaffolds (default: 5)\n"
"   -b  Base name for your output files (optional)\n"
"   -m  Range (in the format min-max) of index multiplicity (only reads with indices in this multiplicity range will be included in graph) (default: 1000-2000)\n"
"   -g  Maximum number of scaffolds in a group (default: 100)\n"
"   -o  Path to previously generated original graph file (optional)\n"
"   -d  Maximum degree of nodes in graph. All nodes with degree greater than this number will be removed from the graph prior to printing final groups. For no node removal, set to 0 (default: 0)\n"
"   -i  Length (bp) of index sequence (default: 14)\n"
"   -v  Runs in verbose mode (optional, default: 0)\n";


ARCS::ArcsParams params;

static const char shortopts[] = "f:a:s:c:l:b:o:m:g:d:i:v";

enum { OPT_HELP = 1, OPT_VERSION};

static const struct option longopts[] = {
    {"file", required_argument, NULL, 'f'},
    {"fofName", required_argument, NULL, 'a'},
    {"seq_id", required_argument, NULL, 's'}, 
    {"min_reads", required_argument, NULL, 'c'},
    {"min_links", required_argument, NULL, 'l'},
    {"base_name", required_argument, NULL, 'b'},
    {"original_file", required_argument, NULL, 'o'},
    {"index_multiplicity", required_argument, NULL, 'm'},
    {"max_groupSize", required_argument, NULL, 'g'},
    {"max_degree", required_argument, NULL, 'd'},
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
 * Check if given SAM flag is one of the accepted ones.
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
/*
    std::smatch sm;
    std::string s = cigar;
    
    while (std::regex_search (s, sm, e)) {
        std::istringstream buffer(sm[1]);
        int value;
        buffer >> value;
        qalen += value;
        
        s = sm.suffix().str();
    }
*/
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

int getIntFromScafName(const std::string& scafName) {
    std::stringstream ss;
    for (auto i = scafName.begin(); i != scafName.end(); ++i) {
        if (isdigit(*i)) {
            ss << *i;
        }
    }
    
    //if (scafName.compare("Super-Scaffold_962476") == 0 || scafName.compare("Super-Scaffold_962661") == 0 || scafName.compare("Super-Scaffold_963983") == 0 || scafName.compare("Super-Scaffold_963663") == 0) {
    //    ss << 0;
    //}

    int numb;
    ss >> numb;
    return numb;
}

/* 
 * Read BAM file, if sequence identity greater than threashold
 * update indexMap. IndexMap also stores information about
 * contig number index algins with and counts.
 */
void readBAM(const std::string bamName, ARCS::IndexMap& imap) {

    /* Open BAM file */
    std::ifstream bamName_stream;
    bamName_stream.open(bamName.c_str());
    if (!bamName_stream) {
        std::cerr << "Could not open " << bamName << ". --fatal.\n";
        exit(EXIT_FAILURE);
    }

    std::string prevRN = "", readyToAddIndex = "";
    int prevSI = 0, prevFlag = 0, prevRef = 0, readyToAddRefName = 0;
    int ct = 1; 

    std::string line;
    int linecount = 0;
    //std::regex e("(\\d+)[M=XI]");
    while (getline(bamName_stream, line)) {
        /* Check to make sure it is not the header */
        if (line.substr(0, 1).compare("@") != 0) {
            linecount++;

            /* Read each line of the BAM file */
            std::stringstream ss(line);
            std::string readName, scafNameString, cigar, rnext, seq, qual;
            int flag, scafName, pos, mapq, pnext, tlen;

            ss >> readName >> flag >> scafNameString >> pos >> mapq >> cigar
                >> rnext >> pnext >> tlen >> seq >> qual;

            scafName = getIntFromScafName(scafNameString);

            /* Parse the index from the readName */
            std::string index = "", readNameSpl = "";
            std::size_t found = readName.find("_");
            if (found!=std::string::npos)
                readNameSpl = readName.substr(found+1);
            if (checkIndex(readNameSpl))
                index = readNameSpl;

            /* Calculate the sequence identity */
            int si = calcSequenceIdentity(line, cigar, seq);

            if (ct >= 3)
                ct = 1;
            if (ct == 1) {
                if (readName.compare(prevRN) != 0) {
                    prevRN = readName;
                    prevSI = si;
                    prevFlag = flag;
                    prevRef = scafName;

                    /* 
                     * Read names are different so we can add the previous index and scafName as
                     * long as there were only two mappings (one for each read)
                     */
                    if (!readyToAddIndex.empty() && readyToAddRefName != 0) {
                        //addToMap(imap, readyToAddIndex, readyToAddRefName);
                        imap[readyToAddIndex][readyToAddRefName]++;
                        readyToAddIndex = "";
                        readyToAddRefName = 0;
                    }
                } else {
                    ct = 0;
                    readyToAddIndex = "";
                    readyToAddRefName = 0;
                }
            } else if (ct == 2) {
                if (prevRN.compare(readName) != 0) {
                    std::cerr << "ERROR! BAM file should be sorted in order of read name. Exiting... \n Prev Read: " << prevRN << "; Curr Read: " << readName << "\n";
                    exit(EXIT_FAILURE);
                }

                if (!seq.empty() && checkFlag(flag) && checkFlag(prevFlag)
                        && si >= params.seq_id && prevSI >= params.seq_id) {
                    if ((prevRef == scafName) && scafName != 0 && !index.empty()) {
                        readyToAddIndex = index;
                        readyToAddRefName = scafName;
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
void readBAMS(const std::string& fofName, ARCS::IndexMap& imap) {

    std::ifstream fofName_stream(fofName.c_str());
    if (!fofName_stream) {
        std::cerr << "Could not open " << fofName << " ...\n";
        exit(EXIT_FAILURE);
    }

    std::string bamName;
    while (getline(fofName_stream, bamName)) {
        if (params.verbose)
            std::cout << "Reading bam " << bamName << std::endl;
        readBAM(bamName, imap);
        assert(fofName_stream);
    }
    fofName_stream.close();
}

/*
 * Calculate how many times given index is seen
 * across all alignments.
 */
int calcMultiplicity(const ARCS::ScafMap sm) {
    int mult = 0;

    ARCS::ScafMap::const_iterator it;
    for(it = sm.begin(); it != sm.end(); ++it) {
        mult += it->second;
    }
    return mult;
}

/* 
 * Iterate through IndexMap and for every pair of scaffolds
 * that align to the same index, store in PairMap. PairMap 
 * is a map with a key of pairs of saffold names, and value
 * of number of links between the pair. (Each link is one index).
 */
void pairContigs(const ARCS::IndexMap& imap, ARCS::PairMap& pmap) {

    /* Iterate through each index in IndexMap */
    for(auto it = imap.begin(); it != imap.end(); ++it) {

        int indexMult = calcMultiplicity(it->second);

        if (indexMult >= params.min_mult && indexMult <= params.max_mult) {

           /* Iterate through all the scafNames in ScafMap */ 
            for (auto o = it->second.begin(); o != it->second.end(); ++o) {
                for (auto p = it->second.begin(); p != it->second.end(); ++p) {

                    /* Only insert into PairMap if o->first less than p->first to avoid duplicates */
                    if (o->second >= params.min_reads && p->second >= params.min_reads 
                            && o->first < p->first) {
                        std::pair<int, int> pair (o->first, p->first);
                        pmap[pair]++;
                    }
                }
            }
        }
    }
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
        if (inserted)
            g[e].weight = it->second;
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
    dp.property("node_id", get(boost::vertex_index, g));
    boost::write_graphviz_dp(out, g, dp);
    assert(out);

    out.close();
}

void readGraph(const std::string& graphFile, ARCS::Graph& g) {
    std::ifstream in(graphFile.c_str());
    assert(in);

    boost::dynamic_properties dp;
    dp.property("id", boost::get(&ARCS::VertexProperties::id, g));
    dp.property("weight", boost::get(&ARCS::EdgeProperties::weight, g));
    dp.property("node_id", boost::get(boost::vertex_index, g));
    bool status = boost::read_graphviz(in, g, dp);
    if (!status) {
        std::cerr << "Error reading " << graphFile
            << "...fatal.\n";
        exit(EXIT_FAILURE);
    }

    std::cout << "      Finished reading graph file." << std::endl;
    in.close();
}

/*
 * Output a .gv file with all the nodes and edges in the graph.
 */
/*
void buildGroups(const ARCS::PairMap& pmap, const std::string graphFile_checkpoint) {
    std::ofstream graphFile_checkpoint_stream(graphFile_checkpoint);
    if (!graphFile_checkpoint_stream) {
        std::cerr << "Could not open " << graphFile_checkpoint 
            << " for writting...fatal.\n";
        exit(EXIT_FAILURE);
    }

    graphFile_checkpoint_stream << "graph scafGroups{\n" <<
        "\tlabel=\"Scaffold Groupings\"\n\trankdir=LR;\n\tnode [shape = circle];\n\toverlap=scale;\n";
    assert(graphFile_checkpoint_stream);

    ARCS::PairMap::const_iterator it;
    for (it = pmap.begin(); it != pmap.end(); ++it) { 
        int scaf1, scaf2;
        std::tie (scaf1, scaf2) = it->first;
        if (it->second >= params.min_links) {
            graphFile_checkpoint_stream << "\t" << scaf1 
                << " [style=filled, fillcolor=deepskyblue, color=deepskyblue]\n";
            graphFile_checkpoint_stream << "\t" << scaf2
                << "[style=filled, fillcolor=deepskyblue, color=deepskyblue]\n";
            graphFile_checkpoint_stream << "\t" << scaf1 << " -- " 
                << scaf2 << " [label= \"l=" << it->second 
                << " \", penwidth=2.0, color=deepskyblue]\n";
            assert(graphFile_checkpoint_stream);
        }
    }
    graphFile_checkpoint_stream << "}\n";
    graphFile_checkpoint_stream.close();
} 
*/

void removeDegreeNodes(ARCS::Graph& g, int max_degree, bool greater) {

    boost::graph_traits<ARCS::Graph>::vertex_iterator vi, vi_end, next;

    boost::tie(vi, vi_end) = boost::vertices(g);
    /*
    for (next = vi; vi != vi_end; vi = next) {
        ++next;
        if (static_cast<int>(boost::degree(*vi, g)) > params.max_degree) {
            boost::clear_vertex(*vi, g);
            //boost::remove_vertex_and_renumber_indices(vi, g);
            boost::remove_vertex(*vi, g);
        }
    }
    
    boost::renumber_indices(g);
    */
    std::vector<ARCS::VertexDes> dVertex;
    for (next = vi; vi != vi_end; vi = next) {
        ++next;
        if (greater) {
            if (static_cast<int>(boost::degree(*vi, g)) > max_degree) {
                dVertex.push_back(*vi);
            }
        } else {
            if (static_cast<int>(boost::degree(*vi, g)) < max_degree) {
                dVertex.push_back(*vi);
            }
        }
    }

    for (unsigned i = 0; i < dVertex.size(); i++) {
        boost::clear_vertex(dVertex[i], g);
        boost::remove_vertex(dVertex[i], g);
    }
    boost::renumber_indices(g);

}

void removeWeightEdges(ARCS::Graph& g, int min_links) {

    boost::graph_traits<ARCS::Graph>::edge_iterator ei, ei_end, next;

    boost::tie(ei, ei_end) = boost::edges(g);
    std::vector<ARCS::VertexDes> dVertex;
    for (next = ei; ei != ei_end; ei = next) {
        ++next;
        if (g[*ei].weight < min_links) {
            boost::remove_edge(*ei, g);
        }
    }
    boost::renumber_indices(g);
    removeDegreeNodes(g, 1, false);
}



void getConnectedComponents(ARCS::Graph& g, ARCS::VertexDescMap& compMap) {
    boost::associative_property_map<ARCS::VertexDescMap> componentMap(compMap);
    std::size_t compNum = boost::connected_components(g, componentMap);
    if (params.verbose)
        std::cout << "Number of connected components: " << compNum << std::endl; 
}


void writePostRemovalGraph(ARCS::Graph& g, const std::string postRemoval) {
    if (params.min_links != 0) {
        std::cout << "      Deleting edges with weight < " << params.min_links <<"... \n";
        removeWeightEdges(g, params.min_links);
    } else {
        std::cout << "      Min links (-l) set to: " << params.min_links << ". Will not delete any edges from graph.\n";
    }

    if (params.max_degree != 0) {
        std::cout << "      Deleting nodes with degree > " << params.max_degree <<"... \n";
        removeDegreeNodes(g, params.max_degree, true);
        //removeDegreeNodes(g, 1, false);
    } else {
        std::cout << "      Max Degree (-d) set to: " << params.max_degree << ". Will not delete any verticies from graph.\n";
    }

    std::cout << "      Writting graph file to " << postRemoval << "...\n";
    writeGraph(postRemoval, g);
}

/*
 * Output a fasta file with the group each scaffold belongs too
 * written in the header.
 */
void writeScafGroups(ARCS::Graph& g, const std::string file, const std::string scafOut) {

    std::ofstream scafOut_stream(scafOut);
    if (!scafOut_stream) {
        std::cerr << "Could not open " << scafOut
            << " for writting...fatal.\n";
        exit(EXIT_FAILURE);
    }

    ARCS::VertexDescMap compMap;
    getConnectedComponents(g, compMap);

    boost::graph_traits<ARCS::Graph>::vertices_size_type numVertex = boost::num_vertices(g);

    ARCS::VidVdesMap vmap;
    std::unordered_map<int, int> componentSize;

    for (int i = 0; i < int(numVertex); i++) {
        ARCS::VertexDes gv = boost::vertex(i, g);
        int vid = g[gv].id;
        vmap[vid] = gv;

        if (compMap.count(gv) != 0) {
            int compNum = compMap[gv];
            componentSize[compNum]++;
        }
    }
    
    FastaReader in(file.c_str(), FastaReader::FOLD_CASE);
    for (FastaRecord rec; in >> rec;) {

        //std::istringstream ss(rec.id);
        int vertex_id = getIntFromScafName(rec.id);
        //ss >> vertex_id;

        if (vmap.count(vertex_id) != 0) {
            ARCS::VertexDes vertex_des = vmap[vertex_id];
            if (compMap.count(vertex_des) != 0) {

                int compNum = compMap[vertex_des];
                if (componentSize[compNum] < params.max_grpSize && componentSize[compNum] > 1) {

                    /* output FASTQ record */
                    std::ostringstream id;
                    id << vertex_id << "_group" << compNum;
                    rec.id = id.str();
                    rec.comment = "";
                    scafOut_stream << rec;
                    assert(scafOut_stream);
                }
            }
        }
    }
}



void runArcs() {

    std::cout << "Running: " << PROGRAM << " " << VERSION 
        << "\n pid " << ::getpid()
        << "\n -f " << params.file 
        << "\n -s " << params.seq_id 
        << "\n -l " << params.min_links     
        << "\n -c " << params.min_reads 
        << "\n Min index multiplicity: " << params.min_mult 
        << "\n Max index multiplicity: " << params.max_mult 
        << "\n -d " << params.max_degree 
        << "\n -i " << params.indexLen << "\n";

    /* Setting output file names */
    std::string scafOut = params.base_name + "_scaffolds.fa";
    std::string graphFile_postRemoval = params.base_name + "_groups.gv";

    std::string graphFile_original = params.original_file;
    if (params.original_file.empty()) {
        std::ostringstream filename;
        filename << params.file << ".scaff" 
            << "_s" << params.seq_id 
            << "_c" << params.min_reads;
        graphFile_original = filename.str() + "_original.gv";
    }
    
    ARCS::IndexMap imap;
    ARCS::PairMap pmap;
    ARCS::Graph g;

    std::time_t rawtime;

    /* Check if graphFile_checkpoint exists. */
    std::ifstream graphFile_original_stream(graphFile_original.c_str());
    if (graphFile_original_stream.good()) {
        std::cout << "\n=> Graph file " << graphFile_original 
            << " found. Skipping reading BAM, pairing, writing graph file. \n";
        graphFile_original_stream.close();

        time(&rawtime);
        std::cout << "\n=>Reading graph file from " << graphFile_original << "... " << ctime(&rawtime);
        readGraph(graphFile_original, g);

    } else {
        graphFile_original_stream.close();

        time(&rawtime);
        std::cout << "\n=>Starting to read BAM files... " << ctime(&rawtime);
        readBAMS(params.fofName, imap);

        time(&rawtime);
        std::cout << "\n=>Starting pairing of scaffolds... " << ctime(&rawtime);
        pairContigs(imap, pmap);

        time(&rawtime);
        std::cout << "\n=>Starting to create graph... " << ctime(&rawtime);
        createGraph(pmap, g);

        time(&rawtime);
        std::cout << "\n=>Starting to write graph file... " << ctime(&rawtime) << "\n";
        writeGraph(graphFile_original, g);
        //buildGroups(pmap, graphFile_checkpoint);
    }

    time(&rawtime);
    std::cout << "\n=>Starting to create post removal graph file... " << ctime(&rawtime);
    writePostRemovalGraph(g, graphFile_postRemoval);

    time(&rawtime);
    std::cout << "\n=>Writting scaffold groups... " << ctime(&rawtime);
    writeScafGroups(g, params.file, scafOut); 

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
            case 'b':
                arg >> params.base_name; break;
            case 'o':
                arg >> params.original_file; break;
            case 'm': {
                std::string firstStr, secondStr;
                std::getline(arg, firstStr, '-');
                std::getline(arg, secondStr);
                std::stringstream ss;
                ss << firstStr << "\t" << secondStr;
                ss >> params.min_mult >> params.max_mult;
                }
                break;
            case 'g':
                arg >> params.max_grpSize; break;
            case 'd':
                arg >> params.max_degree; break;
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
            << "_l" << params.min_links 
            << "_s" << params.seq_id 
            << "_c" << params.min_reads
            << "_d" << params.max_degree 
            << "_pid" << ::getpid(); 
        params.base_name = filename.str();
    }

    runArcs();

    return 0;
}
