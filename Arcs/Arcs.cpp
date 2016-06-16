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
"   -d  Maximum degree of nodes in graph. All nodes with degree greater than this number will be removed from the graph prior to printing final groups. For no node removal, set to 0 (default: 0)\n"
"   -i  Length (bp) of index sequence (default: 14)\n"
"   -v  Runs in verbose mode (optional, default: 0)\n";


ARCS::ArcsParams params;

static const char shortopts[] = "f:a:s:c:l:b:m:g:d:i:v";

enum { OPT_HELP = 1, OPT_VERSION};

static const struct option longopts[] = {
    {"file", required_argument, NULL, 'f'},
    {"fofName", required_argument, NULL, 'a'},
    {"seq_id", required_argument, NULL, 's'}, 
    {"min_reads", required_argument, NULL, 'c'},
    {"min_links", required_argument, NULL, 'l'},
    {"base_name", required_argument, NULL, 'b'},
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
            std::string readName, cigar, rnext, seq, qual;
            int flag, scafName, pos, mapq, pnext, tlen;

            ss >> readName >> flag >> scafName >> pos >> mapq >> cigar
                >> rnext >> pnext >> tlen >> seq >> qual;

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
void pairContigs(const ARCS::IndexMap& imap, ARCS::PairMap& pmap, ARCS::LinkMap& lmap) {

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
                        lmap[pair].insert(it->first);
                        //lmap[o->first].insert(it->first);
                        //lmap[p->first].insert(it->first);
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
void createGraph(const ARCS::PairMap& pmap, ARCS::Graph& g, ARCS::LinkMap& lmap, ARCS::EdgeMap& emap) {

    ARCS::VidVdesMap vmap;

    ARCS::PairMap::const_iterator it;
    for(it = pmap.begin(); it != pmap.end(); ++it) {
        if (it->second >= params.min_links) {
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
            if (inserted) {
                g[e].weight = it->second;
                emap[e] = lmap[it->first];
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
    dp.property("node_id", get(boost::vertex_index, g));
    boost::write_graphviz_dp(out, g, dp);
    assert(out);

    out.close();
}

/*
void readGraph(const std::string& graphFile_dot, ARCS::Graph& g) {
    std::ifstream in(graphFile_dot.c_str());

    boost::dynamic_properties dp;
    dp.property("id", get(&ARCS::VertexProperties::id, g));
    dp.property("weight", get(&ARCS::EdgeProperties::weight, g));
    dp.property("node_id", get(boost::vertex_index, g));
    boost::read_graphviz(in, g, dp);
    assert(in);

    in.close();
}
*/
    

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

void removeDegreeNodes(ARCS::Graph& g) {

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
        if (static_cast<int>(boost::degree(*vi, g)) > params.max_degree) {
            dVertex.push_back(*vi);
        }
    }

    for (unsigned i = 0; i < dVertex.size(); i++) {
        boost::clear_vertex(dVertex[i], g);
        boost::remove_vertex(dVertex[i], g);
    }
    boost::renumber_indices(g);

}


void getConnectedComponents(ARCS::Graph& g, ARCS::VertexDescMap& compMap) {
    boost::associative_property_map<ARCS::VertexDescMap> componentMap(compMap);
    std::size_t compNum = boost::connected_components(g, componentMap);
    if (params.verbose)
        std::cout << "Number of connected components: " << compNum << std::endl; 
}

/*
 * Output a fasta file with the group each scaffold belongs too
 * written in the header.
 */
void writeScafGroups(ARCS::Graph& g, const std::string file, const std::string scafOut, const std::string postRemoval,  ARCS::LinkMap& lmap) {

    std::ofstream scafOut_stream(scafOut);
    if (!scafOut_stream) {
        std::cerr << "Could not open " << scafOut
            << " for writting...fatal.\n";
        exit(EXIT_FAILURE);
    }

    if (params.max_degree != 0) {
        std::cout << "=> Deleting nodes with degree: " << params.max_degree <<" ... \n";
        removeDegreeNodes(g);
        if (params.verbose)
            std::cout << "Finish removing nodes. Now writting graph..." << std::endl;
        writeGraph(postRemoval, g);
        if (params.verbose)
            std::cout << "Finished writing "<< postRemoval << std::endl;
    } else {
        if (params.verbose)
            std::cout << "Max Degree (-d) set to: " << params.max_degree << ". Will not delete any verticies from graph.\n";
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
    
    std::unordered_map<int, std::set<std::string>> indicies;
    
    FastaReader in(file.c_str(), FastaReader::FOLD_CASE);
    for (FastaRecord rec; in >> rec;) {

        std::istringstream ss(rec.id);
        int vertex_id;
        ss >> vertex_id;

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

                    ARCS::Graph::adjacency_iterator neighbourIt, neighbourEnd;
                    boost::tie(neighbourIt, neighbourEnd) = boost::adjacent_vertices(vertex_des, g);
                    for(; neighbourIt != neighbourEnd; ++neighbourIt) {
                        ARCS::VertexDes vNeighbour = *neighbourIt;
                        int vNeighbourId = g[vNeighbour].id; 
                        std::pair<int, int> pair;
                        if (vNeighbourId < vertex_id) {
                            pair = std::make_pair(vNeighbourId, vertex_id);
                        } else {
                            pair = std::make_pair(vertex_id, vNeighbourId);
                        }
                        indicies[compNum].insert(lmap[pair].begin(), lmap[pair].end());

                    }

                    //ARCS::Graph::out_edge_iterator outEdgeIt, outEdgeEnd;
                    //boost::tie(outEdgeIt, outEdgeEnd) = out_edges(vertex_des, g);
                    //for(; outEdgeIt != outEdgeEnd; ++outEdgeIt) {
                    //    ARCS::Graph::edge_descriptor e = *outEdgeIt;
                    //    indicies[compNum].insert(emap[e].begin(), emap[e].end());
                    //}
                }
            }
        }
    }



    std::ofstream indicies_stream("group_indicies.csv");
    if (!indicies_stream) {
        std::cerr << "Could not open " << "group_indicies.csv"
            << " for writting...fatal.\n";
        exit(EXIT_FAILURE);
    }

    std::unordered_map<int, std::set<std::string>>::iterator it;
    for(it = indicies.begin(); it != indicies.end(); ++it) {
        indicies_stream << it->first << '\t';
        for (std::set<std::string>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2){
            indicies_stream << *it2 << '\t';
        }
        indicies_stream << '\n';
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
    std::string scaffold = params.base_name + ".scaffold";
    std::string scafOut = params.base_name + "_d" + std::to_string(static_cast<long long>(params.max_degree)) + "_scaffolds.fa";
    std::string graphFile_checkpoint = params.base_name + "_groups.gv";
    std::string graphFile_postRemoval = params.base_name + "_groups-postNodeRemoval_d" + std::to_string(static_cast<long long>(params.max_degree)) + ".gv";

    ARCS::IndexMap imap;
    ARCS::PairMap pmap;
    ARCS::Graph g;

    ARCS::LinkMap lmap;
    ARCS::EdgeMap emap;

    std::time_t rawtime;

    /* Check if graphFile_checkpoint exists. */
    std::ifstream graphFile_checkpoint_stream(graphFile_checkpoint.c_str());
    if (graphFile_checkpoint_stream.good()) {
        std::cout << "\n Graph file " << graphFile_checkpoint 
            << " found. Skipping reading BAM, pairing, writing graph file. \n";
        graphFile_checkpoint_stream.close();

        //TODO read .gv file
        //readGraph(graphFile_checkpoint, g);

    } else {
        graphFile_checkpoint_stream.close();

        time(&rawtime);
        std::cout << "=>Starting to read BAM files... " << ctime(&rawtime) << "\n";
        readBAMS(params.fofName, imap);

        time(&rawtime);
        std::cout << "=>Starting pairing of scaffolds... " << ctime(&rawtime) << "\n";
        pairContigs(imap, pmap, lmap);

        time(&rawtime);
        std::cout << "=>Starting to create graph... " << ctime(&rawtime) << "\n";
        createGraph(pmap, g, lmap, emap);

        time(&rawtime);
        std::cout << "=>Starting to write graph file... " << ctime(&rawtime) << "\n";
        writeGraph(graphFile_checkpoint, g);
        //buildGroups(pmap, graphFile_checkpoint);
    }

    time(&rawtime);
    std::cout << "=>Writting scaffold groups... " << ctime(&rawtime) << "\n";
    writeScafGroups(g, params.file, scafOut, graphFile_postRemoval, lmap); 

    time(&rawtime);
    std::cout << "=>Done. " << ctime(&rawtime) << "\n";
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
            << "_pid" << ::getpid(); 
        params.base_name = filename.str();
    }

    runArcs();

    return 0;
}
