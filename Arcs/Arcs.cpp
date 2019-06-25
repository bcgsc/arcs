#include "config.h"
#include "Arcs.h"
#include "Arcs/DistanceEst.h"
#include "Common/ContigProperties.h"
#include "Common/Estimate.h"
#include "Common/SAM.h"
#include "Common/StringUtil.h"
#include "Graph/ContigGraph.h"
#include "Graph/DirectedGraph.h"
#include "Graph/DotIO.h"
#include <algorithm>
#include <cassert>
#include <string>
#include <utility>
#include "kseq.h"       // -- newly added
#include <omp.h>        // --newly added
#include <zlib.h>       // -- newly added for gzopen 

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
"Usage: arcs [OPTION]... ALIGNMENTS...\n"
"\n"
"ALIGNMENTS may be a SAM or BAM file.\n"
"The output of the aligner may be piped directly into ARCS by setting\n"
"ALIGNMENTS to /dev/stdin, in which case it must be in SAM format.\n"
"\n"
"Paired reads must occur consecutively (interleaved) in the SAM/BAM file.\n"
"The output of the aligner may either not be sorted,\n"
"or may be sorted by read name using samtools sort -n.\n"
"The SAM/BAM file must not be sorted by coordinate position.\n"
"\n"
"The barcode may be found in either the BX:Z:BARCODE SAM tag,\n"
"or in the read (query) name following an underscore, READNAME_BARCODE.\n"
"In the latter case the barcode must be compsed entirely of nucleotides.\n"
"\n"
"The contig sequence lengths must either be present in the SAM header\n"
"or provided by the -f option.\n"
"\n"
" Options:\n"
"\n"
"   -f, --file=FILE       FASTA file of contig sequences to scaffold [optional]\n"
"   -a, --fofName=FILE    text file listing input SAM/BAM filenames\n"
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

static const char shortopts[] = "f:a:B:s:c:Dl:z:b:g:m:d:e:r:v:t:x:p:u:q:w:i:o:k:h:j";   // --newly added part t:x:p:u:q:w:i:o:k:h:j

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
    // possibly new options must be put here as well.
static const struct option longopts[] = {
    {"file", required_argument, NULL, 'f'},
    {"fofName", required_argument, NULL, 'a'},
    {"bin_size", required_argument, NULL, 'B'},
    {"bx", no_argument, NULL, OPT_BX }, // ignored
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
    {"threads", required_argument, NULL, 't'},      //new part
    {"kmer_method", required_argument, NULL, 'x'},
    {"program", required_argument, NULL, 'p'},
    {"multfile", required_argument, NULL, 'u'},
    {"conrecfile", required_argument, NULL, 'q'},
    {"kmapfile", required_argument, NULL, 'w'},
    {"imapfile", required_argument, NULL, 'i'},
    {"checkpoint_outs", required_argument, NULL, 'o'},
    {"k_value", required_argument, NULL, 'k'},
    {"k_shift", required_argument, NULL, 'h'},
    {"j_index", required_argument, NULL, 'j'},
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

/** Added for program flow of ARKS(kmer) */
bool full = false;
bool alignc = false;
bool graph = false;

/** Added for program flow of ARKS(kmer) */
unsigned int s_numkmersmapped = 0, s_numkmercollisions = 0, s_numkmersremdup = 0,
		s_numbadkmers = 0, s_uniquedraftkmers = 0;

unsigned int s_totalnumckmers = 0, s_ckmersasdups = 0, s_numckmersfound = 0,
		s_numckmersrec = 0, s_numbadckmers = 0;

unsigned int s_numreadspassingjaccard = 0, s_numreadsfailjaccard = 0;

KSEQ_INIT(gzFile, gzread)       // --newly added

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

/** Check if SAM flag is one of the accepted ones. */
static inline bool checkFlag(int flag)
{
    flag &= ~0xc0; // clear READ1,READ2
    return flag == 19 // PAIRED,PROPER_PAIR,REVERSE
        || flag == 35; // PAIRED,PROPER_PAIR,MREVERSE
}

/*
 * Check if character is one of the accepted ones.
 */
bool checkChar(char c) {
    return (c == 'M' || c == '=' || c == 'X' || c == 'I');
}

//for multiple fastq chromium files     --- added from arks
vector<string> convertInputString(const string &inputString) {
	vector<string> currentInfoFile;
	string temp;
	stringstream converter(inputString);
	while (converter >> temp) {
		currentInfoFile.push_back(temp);
	}
	return currentInfoFile;
}

/* Calculate the jaccard index */
static inline double calcJacIndex(int smallCount, int overallCount) {
	return (double) smallCount / (double) overallCount;
}

/* Strip the trailing "/1" or "/2" from a FASTA ID, if such exists */
static inline void stripReadNum(std::string& readName)
{
	size_t pos = readName.rfind("/");
	if (pos == std::string::npos)
		return;
	if (pos == 0 || pos == readName.length() - 1)
		return;
	if (!std::isdigit(readName.at(pos + 1)))
		return;
	readName.resize(pos);
}

/* Track memory usage */
int memory_usage() {
	int mem = 0;
	ifstream proc("/proc/self/status");
	string s;
	while (getline(proc, s), !proc.fail()) {
		if (s.substr(0, 5) == "VmRSS") {
			stringstream convert(
					s.substr(s.find_last_of('\t'), s.find_last_of('k') - 1));
			if (!(convert >> mem)) {
				return 0;
			}
			return mem;
		}
	}
	return mem;
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
// -- newly added
/* Returns true if the contig sequence contains ATGC or IUPAC codes */
static inline bool checkContigSequence(std::string seq) {
	for (int i = 0; i < static_cast<int>(seq.length()); i++) {
		char c = toupper(seq[i]);
		if (c != 'A' && c != 'T' && c != 'G' && c != 'C' && c != 'N' && c != 'M'
				&& c != 'R' && c != 'W' && c != 'S' && c != 'Y' && c != 'K'
				&& c != 'V' && c != 'H' && c != 'D' && c != 'B') {
			std::cout << c << std::endl;
			return false;
		}
	}

	return true;
}

// /* Create ContigKmerMap */
// void createContigKmerMap(std::string kmaptsv, ARKS::ContigKMap &kmap) {

// 	std::ifstream kmaptsv_stream;
// 	kmaptsv_stream.open(kmaptsv.c_str());
// 	if (!kmaptsv_stream) {
// 		std::cerr << "Could not open " << kmaptsv << ". --fatal.\n";
// 		exit (EXIT_FAILURE);
// 	}

// 	std::string line;
// 	while (getline(kmaptsv_stream, line)) {
// 		std::stringstream sst(line);

// 		std::string kmer, contigreci_string;
// 		int contigreci;

// 		sst >> kmer >> contigreci_string;
// 		contigreci = std::stoi(contigreci_string);

// 		kmap[kmer] = contigreci;
// 	}
// 	kmaptsv_stream.close();
// }

// /* Create contigRecord vector */
// void createContigRecord(std::string contigrectsv, std::vector<ARKS::CI> &contigRecord) {

// 	std::ifstream contigrectsv_stream;
// 	contigrectsv_stream.open(contigrectsv.c_str());
// 	if (!contigrectsv_stream) {
// 		std::cerr << "Could not open " << contigrectsv << ". --fatal.\n";
// 		exit (EXIT_FAILURE);
// 	}

// 	std::string line;
// 	while (getline(contigrectsv_stream, line)) {
// 		std::stringstream sst(line);

// 		int contigreci;
// 		bool ht;
// 		std::string contigreci_string, contigname, headortail;
// 		sst >> contigreci_string >> contigname >> headortail;
// 		contigreci = std::stoi(contigreci_string);
// 		ht = HTtoBool(headortail);

// 		ARKS::CI contigID(contigname, ht);

// 		contigRecord[contigreci] = contigID;
// 	}
// 	contigrectsv_stream.close();
// }

/* Returns true if seqence only contains ATGC or N
 * 	Allows only if 0.98 of the read is not ambiguous (aka not N)
 *
 *	TODO: Allow ambiguity to be a user defined parameter
 */
static inline bool checkReadSequence(std::string seq) {

	double ambiguity = 0;

	for (int i = 0; i < static_cast<int>(seq.length()); i++) {
		char c = toupper(seq[i]);
		if (c != 'A' && c != 'T' && c != 'G' && c != 'C') {
			if (c == 'N') {
				ambiguity++;
			} else {
				return false;
			}
		}
	}

	double ar = ambiguity / (double) seq.length();

	if (ar > 0.02)
		return false;

	return true;
}

/* Reading TSV (or CSV for barcode mult) checkpoint files */

/* Create indexMultMap from Barcode Multiplicity File */
void createIndexMultMap(std::string multfile, std::unordered_map<std::string, int> &indexMultMap) {

	size_t numreadstotal=0;
	size_t numreadskept=0;
	size_t numbarcodes=0;

	// Decide if it is a tsv or csv file
	bool tsv = false;
	std::size_t found = multfile.find(".tsv");
	if (found!=std::string::npos)
		tsv = true;

	std::ifstream multfile_stream;
	multfile_stream.open(multfile.c_str());
	if (!multfile_stream) {
		std::cerr << "Could not open " << multfile << ". --fatal.\n";
		exit (EXIT_FAILURE);
	}

	std::string line;
	while(getline(multfile_stream, line)) {

		std::string barcode;
		std::string multiplicity_string;
		size_t multiplicity;

		if (tsv) {
			std::stringstream sst(line);
			sst >> barcode >> multiplicity_string;
			numbarcodes++;
			multiplicity = std::stoi(multiplicity_string);
		} else {
			std::istringstream iss(line);
			getline(iss, barcode, ',');
			iss >> multiplicity_string;
			numbarcodes++;
			multiplicity = std::stoi(multiplicity_string);
		}

		numreadstotal += multiplicity;

		if (!barcode.empty()) {
			numreadskept += multiplicity;
			indexMultMap[barcode] = multiplicity;
		} else {
			std::cout << "Please check your multiplicity file." << std::endl;
		}

		assert(multfile_stream);
	}
	multfile_stream.close();

	if (params.verbose) {
		std::cout << "Saw " << numbarcodes << " barcodes and keeping " << numreadskept << " read pairs out of " << numreadstotal << std::endl;
	}

}

/* Returns the size of the array for storing contigs */
size_t initContigArray(std::string contigfile) {

	size_t count = 0;

	gzFile fp;

	int l;
	const char* filename = contigfile.c_str();
	fp = gzopen(filename, "r");     // this method is copied from arks and not sure if dependencies of gzopne exist.
	kseq_t * seq = kseq_init(fp);

	while ((l = kseq_read(seq)) >= 0) {
		std::string sequence = seq->seq.s;
		unsigned sequence_length = sequence.length();
		if (checkContigSequence(sequence) && static_cast<int>(sequence_length) >= params.min_size)
			count++;
	}
	kseq_destroy(seq);
	gzclose(fp);
	
	//*
	// Not sure why it is multiplied by two.
	//*
	if (params.verbose) {
		cerr << "Number of contigs:" << count << "\nSize of Contig Array:"
				<< (count * 2) + 1 << endl;
	}

	//Let the first index represent the a null contig
	return (count * 2) + 1;
}


/* Get all scaffold sizes from FASTA file */
void getScaffSizes(std::string file, ARCS::ScaffSizeList& scaffSizes, ARCS::ContigToLength& contigToLength) {

    int counter = 0;
    FastaReader in(file.c_str(), FastaReader::FOLD_CASE);
    for (FastaRecord rec; in >> rec;) {
        counter++;
        std::string  scafName = rec.id;
        int size = rec.seq.length();
        scaffSizes.push_back(std::make_pair(rec.id, size));
        contigToLength[rec.id] = size;          //newly added
    }

    if (params.verbose)
        std::cout << "Saw " << counter << " sequences.\n";
}

/*
 * Read BAM file, if sequence identity greater than threashold
 * update indexMap. IndexMap also stores information about
 * contig number index algins with and counts.
 */
void readBAM(const std::string bamName, ARCS::IndexMap& imap, std::unordered_map<std::string, int>& indexMultMap,
        ARCS::ScaffSizeList& scaffSizeList, ARCS::ScaffSizeMap& sMap, ARCS::ContigToLength& contigToLengthMap)
{
    /* Open BAM file */
    std::ifstream bamName_stream;
    bamName_stream.open(bamName.c_str());
    assert_good(bamName_stream, bamName);
    if (bamName_stream.peek() == EOF) {
        std::cerr << "error: alignments file is empty: " << bamName << '\n';
        exit(EXIT_FAILURE);
    }

    std::string prevRN = "", readyToAddIndex = "", prevRef = "", readyToAddRefName = "";
    int prevSI = 0, prevFlag = 0, prevMapq = 0, prevPos = -1, readyToAddPos = -1;
    int ct = 1;

    std::string line;
    size_t linecount = 0;
    const unsigned int suppAligFlag = 0x00000800;
    const unsigned int notPrimaryAligFlag = 0x00000100;

    // Number of unpaired reads.
    size_t countUnpaired = 0;

    // Whether to add SAM SQ headers to sMap.
    //const bool addSAMSequenceLengths = sMap.empty();    // if scaffold file is not given the length and name of the scaffolds that reads are mapped will
                                                        // will be taken from the SAM header.
    const bool addSAMSequenceLengths = contigToLengthMap.empty();

    /* Read each line of the BAM file */
    while (getline(bamName_stream, line)) {
        if (line.empty())
            continue;
        if (line[0] == '@') {
            // Parse the SAM header.
            if (startsWith(line, "@SQ\t")) {
                std::stringstream ss(line);
                std::string name;
                size_t size = 0;
                ss >> expect("@SQ\tSN:") >> name >> expect("\tLN:") >> size;
                if (!ss) {
                    std::cerr << "error: parsing SAM header: " << line << '\n';
                    exit(EXIT_FAILURE);
                }
                if (addSAMSequenceLengths) {
                    //ARCS::ScaffSizeMap::value_type sq_ln(name, size);
                    ARCS::ContigToLength::value_type sq_ln(name, size);
                    scaffSizeList.push_back(sq_ln);
                    sMap.insert(sq_ln);
                    contigToLengthMap.insert(sq_ln);
                } else {
                    // auto it = sMap.find();
                    // if (it == sMap.end()) {
                    auto it = contigToLengthMap.find(name);  
                    if (it == contigToLengthMap.end()) { 
                        std::cerr << "error: unexpected sequence: " << name << " of size " << size;
                        exit(EXIT_FAILURE);
                    } else if (it->second != (int)size) {
                        std::cerr << "error: mismatched sequence lengths: sequence "
                            << name << ": " << it->second << " != " << size;
                        exit(EXIT_FAILURE);
                    }
                }
            }
        } else {
            linecount++;

            std::stringstream ss(line);
            std::string readName, scafName, cigar, rnext, seq, qual, tags;
            int flag, pos, mapq, pnext, tlen;

            ss >> readName >> flag >> scafName >> pos >> mapq >> cigar
                >> rnext >> pnext >> tlen >> seq >> qual >> std::ws;

            getline(ss, tags);

            /* Parse the index from the readName */
            std::string index = parseBXTag(tags);
            if (index.empty()) {
                std::size_t found = readName.rfind("_");
                if (found != std::string::npos) {
                    index = readName.substr(found + 1);
                    // Check that the barcode is composed of only ACGT.
                    if (index.find_first_not_of("ACGTacgt") != std::string::npos)
                      index.clear();
                }
            }
            
            /* Keep track of index multiplicity if alignment is not supplementary alignment and is not non-primary alignment*/
            if (!index.empty() && !((flag & suppAligFlag) || (flag & notPrimaryAligFlag)))
                indexMultMap[index]++;

            /* Calculate the sequence identity */
            int si = calcSequenceIdentity(line, cigar, seq);

            if (ct == 2 && readName != prevRN) {
                if (countUnpaired == 0)
                    std::cerr << "Warning: Skipping an unpaired read. Read pairs should be consecutive in the SAM/BAM file.\n"
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


                    /*
                     * Read names are different so we can add the previous index and scafName as
                     * long as there were only two mappings (one for each read)
                     */
                    if (!readyToAddIndex.empty() && !readyToAddRefName.empty() && readyToAddRefName.compare("*") != 0 && readyToAddPos != -1) {

                        //int size = sMap[readyToAddRefName];
                        int size = contigToLengthMap[readyToAddRefName];   
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
        std::cerr << "Warning: Skipped " << countUnpaired << " unpaired reads. Read pairs should be consecutive in the SAM/BAM file.\n";
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
    for (std::string singleFileName; getline(fin, singleFileName);)
        filenames.push_back(singleFileName);
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
void readBAMS(const std::vector<std::string> bamNames, ARCS::IndexMap& imap, std::unordered_map<std::string, int>& indexMultMap,
        ARCS::ScaffSizeList& scaffSizeList, ARCS::ScaffSizeMap& scaffSizeMap, ARCS::ContigToLength& contigToLengthMap)
{
    assert(!bamNames.empty());
    for (const auto& bamName : bamNames) {
        if (params.verbose)
            std::cout << "Reading alignments: " << bamName << std::endl;
        readBAM(bamName, imap, indexMultMap, scaffSizeList, scaffSizeMap, contigToLengthMap);
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

// --newly added
/* Shreds end sequence into kmers and inputs them one by one into the ContigKMap
 * 	std::pair<std::string, bool> 				specifies contigID and head/tail
 * 	std::string						the end sequence of the contig
 *	int							k-value specified by user
 *	ARKS::ContigKMap					ContigKMap for storage of kmers
 */
int mapKmers(std::string seqToKmerize, int k, int k_shift,
		ARCS::ContigKMap &kmap, ReadsProcessor &proc, int conreci) {

	int seqsize = seqToKmerize.length();

	// Checks if length of the subsequence is smaller than the k size
	// 	If the contig end is shorter than the k-value, then we ignore.
	if (seqsize < k) {
		std::string errmsg =
				"Warning: ends of contig is shorter than k-value for contigID (no k-mers added): ";
		std::cout << errmsg << conreci << std::endl;

		return 0;
	} else {
		//assert(seqsize >= k);
		int numKmers = 0;

		int i = 0;
		while (i <= seqsize - k) {
			//*
			// This proc.prepSeq() must have been returning the kmer as string.
			//
			const unsigned char* temp = proc.prepSeq(seqToKmerize, i); // prepSeq returns NULL if it contains an N
			// Ignore a NULL kmer
			if (temp != NULL) {
				std::string kmerseq = proc.getStr(temp);

				numKmers++;

				//Critical section to prevent multiple thread from accessing datastructure
				//only needed since we are writing to datastructure
				//ideally would want CAS operation but without special a datastructure this is difficult
				bool exists;
				int alreadyconreci;

//pragma omp critical(exists)
				exists = kmap.find(kmerseq) != kmap.end();
				//*
				// Wouldn't this line give an error if kmerseq is not in the
				// kmap?
				alreadyconreci = kmap[kmerseq];

				if (exists) {
					if (alreadyconreci != conreci) {
//#pragma omp atomic
						s_numkmersremdup++;
						if (alreadyconreci != 0) {
							s_uniquedraftkmers--;
//#pragma omp critical(duplicate)
							//*
							// Why setting it to zero if another
							// contig also has the kmer?
							//*
							kmap[kmerseq] = 0;
						}
					}
//#pragma omp atomic
					s_numkmercollisions++;
				} else {
//#pragma omp critical(insertkmerseq)
					kmap[kmerseq] = conreci;
					s_uniquedraftkmers++;
					s_numkmersmapped++;
				}
				i += k_shift;
			} else {
				i += k;
//#pragma omp atomic
				s_numbadkmers++;
			}
		}
		return numKmers;
	}
}

/* Returns best corresponding contig from read through kmers
 * 	ARKS::ContigKMap			tells me what kmers correspond to which contig
 *	std::string				read sequence
 *	int					size of k-mer
 *	int 					k_shift
 *      double j_index				Jaccard Index (default 0.5)
 *	ReadsProcessor				kmerizer
 */
int bestContig (ARCS::ContigKMap &kmap, std::string readseq, int k, int k_shift,
		double j_index, ReadsProcessor &proc) {

	// to keep track of what contig+H/T that the k-mer from barcode matches to
	// 	int					Index that corresponds to the contig in the contigRecord
	// 	count					# kmers found here
	std::map<int, int> ktrack;

	// k-merize readsequence
	int corrbestConReci = 0;
	int seqlen = readseq.length();

	int totalnumkmers = 0;
	int kmerdups = 0;
	int kmerfound = 0;
	int kmerstore = 0;

    //std::cout << "here! 10" << std::endl;

    //std::cout << "thread num: " << omp_get_thread_num() << "proc mem :" << &proc << std::endl;

	int i = 0;

	while (i <= seqlen - k) {
        //std::cout << "here! 11" << std::endl;
		const unsigned char* temp = proc.prepSeq(readseq, i);
        //std::cout << "here! 20" << std::endl;
#pragma omp atomic    //--commented
		totalnumkmers++;
		if (temp != NULL) {
            //std::cout << "here! 12" << std::endl;
			const std::string ckmerseq = proc.getStr(temp);
#pragma omp atomic    //--commented
			s_totalnumckmers++;

			// search for kmer in ContigKmerMap and only record if it is not the collisionmaker
			if (kmap.find(ckmerseq) != kmap.end()) {
                //std::cout << "here! 5" << std::endl;
				int corrConReci = kmap[ckmerseq];
                //std::cout << "here! 6" << std::endl;
				if (corrConReci != 0) {
                //std::cout << "here! 7" << std::endl;
					ktrack[corrConReci]++;
                //std::cout << "here! 8" << std::endl;
#pragma omp atomic    //--commented
					kmerstore++;
#pragma omp atomic    //--commented
					s_numckmersrec++;
				} else {
#pragma omp atomic    //--commented
					s_ckmersasdups++;
#pragma omp atomic    //--commented
					kmerdups++;
				}
#pragma omp atomic    //--commented
				kmerfound++;
#pragma omp atomic    //--commented
				s_numckmersfound++;

			}
		} else {
#pragma omp atomic    //--commented
			s_numbadckmers++;
		}
        //std::cout << "here! 9" << std::endl;
		i += k_shift;
        // //* put here for debugging purposes
        // if (i-k_shift / 50 != i / 50 )
        // {
        //     std::cout << "i: " << i << std::endl;
        // }
        if (s_numbadckmers % 100000 == 1)
        {
            std::cout << "s_numbadckmers: " << s_numbadckmers << std::endl;
        }
        
        
	}
    //std::cout << "here! 13" << std::endl;
	double maxjaccardindex = 0;
	// for the read, find the contig that it is most compatible with based on the jaccard index
	for (auto it = ktrack.begin(); it != ktrack.end(); ++it) {
		double currjaccardindex = calcJacIndex(it->second, totalnumkmers);
		if (maxjaccardindex < currjaccardindex) {
			maxjaccardindex = currjaccardindex;
			corrbestConReci = it->first;
		}
	}

	// default jaccard threshold is 0.5
	if (maxjaccardindex > j_index) {
		s_numreadspassingjaccard++;
		return corrbestConReci;
	} else {

		s_numreadsfailjaccard++;
		return 0;
	}
}

/// --newly added
/* Get the k-mers from the paired ends of the contigs and store them in map.
 * 	std::string file					FASTA (or later FASTQ) file
 *	std::sparse_hash_map<k-mer, pair<contidID, bool>> 	ContigKMap
 *	int k							k-value (specified by user)
 */
void getContigKmers(std::string contigfile, ARCS::ContigKMap &kmap,
	std::vector<ARCS::CI> &contigRecord, ARCS::ContigToLength& contigToLength)
{
	int totalNumContigs = 0;
	int skippedContigs = 0;
	int validContigs = 0;
	int totalKmers = 0;

	ARCS::CI collisionmarker("null contig", false);
	size_t conreci = 0; // 0 is the null contig so we will later increment before adding
	contigRecord[conreci] = collisionmarker;

	gzFile fp;

	int l;
	const char* filename = contigfile.c_str();
	fp = gzopen(filename, "r");
	kseq_t * seq = kseq_init(fp);

	//each thread gets a proc;
	//vector<ReadsProcessor*> procs(params.threads);

	int16_t k_proc = params.k_value;
	ReadsProcessor proc(k_proc);

//	for (unsigned i = 0; i < params.threads; ++i) {
//		procs[i] = new ReadsProcessor(params.k_value);
//	}

//#pragma omp parallel
	while ((l = kseq_read(seq)) >= 0) {
		totalNumContigs++;
		bool good = false;
		std::string contigID = "", sequence = "";
		size_t tempConreci1 = 0;
		size_t tempConreci2 = 0;
//#pragma omp critical(kseq)
		{
			//l = kseq_read(seq);
			//if (l >= 0) {
				contigID = seq->name.s;
				sequence = seq->seq.s;
				if (static_cast<int>(sequence.length()) >= params.min_size) {
					tempConreci1 = ++conreci;
					tempConreci2 = ++conreci;
					good = true;
				}
			//}
		}


		if (good) {
//			if (!checkContigSequence(sequence)) {
//				std::string errormsg =
//						"Error: Contig contains non-base characters. Please check your draft genome input file.";
//				if (params.verbose) {
//					std::cerr << contigID << ": " << errormsg << std::endl;
//				}
//#pragma omp atomic
//				skippedContigs++;
//			} else {
				// If the sequence is above minimum contig length, then will extract kmers from both ends
				// If not will ignore the contig
				int sequence_length = sequence.length();

				// record contig length for later use (distance estimation)
				contigToLength[contigID] = sequence_length;

				// If contig length is less than 2 x end_length, then we split the sequence
				// in half to decide head/tail (aka we changed the end_length)
				int cutOff = params.end_length;
				if (cutOff == 0 || sequence_length <= cutOff * 2)
					cutOff = sequence_length / 2;

				// Arbitrarily assign head or tail to ends of the contig
				ARCS::CI headside(contigID, true);
				ARCS::CI tailside(contigID, false);

				//get ends of the sequence and put k-mers into the map
				contigRecord[tempConreci1] = headside;		// in contigRecord head/tail of contig are kept as different contigs.
				std::string seqend;
				seqend = sequence.substr(0, cutOff);
				int num = mapKmers(seqend, params.k_value, params.k_shift,
						kmap, proc, tempConreci1);
//#pragma omp atomic
				totalKmers += num;

				contigRecord[tempConreci2] = tailside;
				seqend = sequence.substr(sequence_length - cutOff,
						sequence_length);
				num = mapKmers(seqend, params.k_value, params.k_shift, kmap,
						proc, tempConreci2);

//#pragma omp atomic
				totalKmers += num;
//#pragma omp atomic
				validContigs++;
		} else {
//#pragma omp atomic
			skippedContigs++;
		}
//			}
//		}
		// printprogress
		if (params.verbose) {
//#pragma omp critical(stdout)
			if (totalNumContigs % 1000 == 0) {

				printf("Finished %d Contigs...\n", totalNumContigs);
				// for memory tracking + debugging usage:
					//std::cout << "Cumulative memory usage: " << memory_usage() << std::endl;
					//std::cout << "Kmers so far: " << s_numkmersmapped << std::endl;
			}
		}
	}
	kseq_destroy(seq);
	gzclose(fp);

	// clean up
//	delete proc;
//	for (unsigned i = 0; i < params.threads; ++i) {
//		delete procs[i];
//	}

	if (params.verbose) {
		printf(
				"%s %u\n%s %u\n%s %u\n%s %u\n%s %u\n%s %u\n%s %u\n%s %u\n%s %u\n",
				"Total number of contigs in draft genome: ", totalNumContigs,
				"Total valid contigs: ", validContigs,
				"Total skipped contigs: ", skippedContigs,
				"Total number of Kmers: ", totalKmers, "Number Null Kmers: ",
				s_numbadkmers, "Number Kmers Recorded: ", s_numkmersmapped,
				"Number Kmer Collisions: ", s_numkmercollisions,
				"Number Times Kmers Removed (since duplicate in different contig): ",
				s_numkmersremdup, "Number of unique kmers (only one contig): ",
				s_uniquedraftkmers);
	}
}

/* Read through longranger basic chromium output fastq file */
void chromiumRead(std::string chromiumfile, ARCS::ContigKMap& kmap, ARCS::IndexMap& imap,
			const std::unordered_map<std::string, int> &indexMultMap,
			const std::vector<ARCS::CI> &contigRecord) {

	int stored_readpairs = 0;
	int skipped_unpaired = 0;
	int skipped_invalidreadpair = 0;
	int skipped_nogoodcontig = 0;
	int invalidbarcode = 0;
	int emptybarcode=0; 

	size_t count = 0;
	bool stop = false;

	//each thread gets a proc;
	vector<ReadsProcessor*> procs(params.threads);


	for (unsigned i = 0; i < params.threads; ++i) {
		procs[i] = new ReadsProcessor(params.k_value);
	}
    
	gzFile fp2;
	const char* filename = chromiumfile.c_str();
	fp2 = gzopen(filename, "r");
	if (fp2 == Z_NULL) {
		cerr << "File " << filename << " cannot be opened." << endl;
		exit(1);
	} else {
		cerr << "File " << filename << " opened." << endl;
	}
	kseq_t * seq2 = kseq_init(fp2);

#pragma omp parallel      //--commented
	while (!stop) {
		int l;
		std::string read1_name = "";
		std::string read2_name = "";
		std::string barcode1;
		std::string barcode2;
		std::string cread1 = "";
		std::string cread2 = "";
		std::string comment1 = "";
		std::string comment2 = "";
	        std::size_t foundTag;
	        std::size_t foundEnd;
		bool paired = false;
		int corrConReci1 = 0;
		int corrConReci2 = 0;
#pragma omp critical(checkread1or2) //chechkread1or2 doesn't exist. //--commented
		{
            //std::cout << "here! 14" << std::endl;
			l = kseq_read(seq2);
			if (l >= 0) {
				read1_name = seq2->name.s;
				if (seq2->comment.l) {
					comment1 = seq2->comment.s;
				}
				cread1 = seq2->seq.s;
				l = kseq_read(seq2);
				if (l >= 0) {
					read2_name = seq2->name.s;
					if (seq2->comment.l) {
						comment2 = seq2->comment.s;
					}
					cread2 = seq2->seq.s;
				} else {
					stop = true;
				}
			} else {
				stop = true;
			}
			stripReadNum(read1_name);
			stripReadNum(read2_name);
			if (read1_name == read2_name) {
				paired = true;
			} else {
				std::cout << "File contains unpaired reads: " << read1_name << " " << read2_name << std::endl;
				skipped_unpaired++;
			}
			count += 2;
			if (params.verbose) {
				if (count % 100000 == 0) {
					std::cout << "Processed " << count << " read pairs." << std::endl;
				}
			}
		}

		if (!stop) {
            //std::cout << "here! 15" << std::endl;
			barcode1.clear();
	    		//Find position of BX:Z:
			foundTag = comment1.find("BX:Z:");
	    		if(foundTag != std::string::npos){
				// End is space if there is another tag, newline otherwise
		                foundEnd = comment1.find(' ', foundTag);
				// Get substring from end of BX:Z: to space or end of string
                		if(foundEnd != std::string::npos){
		    			barcode1 = comment1.substr(foundTag + 5, foundEnd - foundTag - 5);
				}
				else {
		    			barcode1 = comment1.substr(foundTag + 5);
				}
	    		}
			
			barcode2.clear();
	    		//Find position of BX:Z:
			foundTag = comment2.find("BX:Z:");
	    		if(foundTag != std::string::npos){
				// End is space if there is another tag, newline otherwise
				foundEnd = comment2.find(' ', foundTag);
				// Get substring from end of BX:Z: to space or end of string
		                if(foundEnd != std::string::npos){
		    			barcode2 = comment2.substr(foundTag + 5, foundEnd - foundTag - 5);
				}else {
		    			barcode2 = comment2.substr(foundTag + 5);
				}
	    		}

			bool validbarcode = false; 
			if (barcode1.empty() || barcode2.empty()) {
				emptybarcode++; 
			} else {					// why taking the barcode mult file as input and checking again while reading the chromium file inside.
										// its not used anywhere before this readChrom function and can be created here.
				validbarcode = indexMultMap.find(barcode1) != indexMultMap.end();			
				if (!validbarcode) {
#pragma omp atomic    //--commented
					invalidbarcode++;
				}
			}

			if (paired && validbarcode && !barcode1.empty() && !barcode2.empty() && (barcode1==barcode2)) {
				const int indexMult = indexMultMap.at(barcode1);
				bool goodmult = indexMult > params.min_mult || indexMult < params.max_mult;
				if (goodmult && checkReadSequence(cread1) && checkReadSequence(cread2)) {
                    //std::cout << "Here! 3" << std::endl;
					corrConReci1 = bestContig(kmap, cread1, params.k_value, params.k_shift, params.j_index, *procs[omp_get_thread_num()]);
					corrConReci2 = bestContig(kmap, cread2, params.k_value, params.k_shift, params.j_index, *procs[omp_get_thread_num()]);
                    //std::cout << "Here! 4" << std::endl;
                } else {
#pragma omp atomic    //--commented
					skipped_invalidreadpair++;
				}
				// we only store barcode info in index map if read pairs have same contig + orientation
				// and if the corrContigId is not NULL (because it is above accuracy threshold)
				if (corrConReci1 != 0 && corrConReci1 == corrConReci2) {
					const ARCS::CI corrContigId = contigRecord[corrConReci1];
#pragma omp critical(imap)    //--commented
					{
                        //std::cout << "Here! 1" << std::endl;
						imap[barcode1][corrContigId]++;
                        //std::cout << "Here! 2" << std::endl;

					}

#pragma omp atomic    //--commented
					stored_readpairs++;

				} else {
#pragma omp atomic    //--commented
					skipped_nogoodcontig++;
				}
			}
		}
	}
	kseq_destroy(seq2);
	gzclose(fp2);

	// clean up
	for (unsigned i = 0; i < params.threads; ++i) {
		delete procs[i];
	}

	if (params.verbose) {
		printf(
				"Stored read pairs: %u\nSkipped invalid read pairs: %u\nSkipped unpaired reads: %u\nSkipped reads pairs without a good contig: %u\n",
				stored_readpairs, skipped_invalidreadpair, skipped_unpaired,
				skipped_nogoodcontig);
		printf(
				"Total valid kmers: %u\nNumber invalid kmers: %u\nNumber of kmers found in ContigKmap: %u\nNumber of kmers recorded in Ktrack: %u\nNumber of kmers found in ContigKmap but duplicate: %u\nNumber of reads passing jaccard threshold: %u\nNumber of reads failing jaccard threshold: %u\n",
				s_totalnumckmers, s_numbadckmers, s_numckmersfound, s_numckmersrec,
				s_ckmersasdups, s_numreadspassingjaccard, s_numreadsfailjaccard);
		if (emptybarcode>0) 
			printf("WARNING:: Your chromium read file has %d readpairs that have an empty barcode.", emptybarcode); 
		if (invalidbarcode > 0)
			printf("WARNING:: Your chromium read file has %d read pairs that have barcodes not in the barcode multiplicity file.", invalidbarcode);

	}
}

void readChroms(vector<string> inputFiles, ARCS::ContigKMap &kmap,
		ARCS::IndexMap &imap,
		const std::unordered_map<std::string, int> &indexMultMap,
		const std::vector<ARCS::CI> &contigRecord) {

	std::string chromFile;

	for (auto p = inputFiles.begin(); p != inputFiles.end(); ++p) {
		chromFile = *p;
		if (params.verbose)
			std::cout << "Reading chrom " << chromFile << std::endl;
		chromiumRead(chromFile, kmap, imap, indexMultMap, contigRecord);
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
    assert(!graphFile_dot.empty());

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
    assert(!graphFile.empty());

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
//void createAbyssGraph(const ARCS::ScaffSizeList& scaffSizes, const ARCS::Graph& gin, DistGraph& gout) {
void createAbyssGraph(const ARCS::ContigToLength& contigToLength, const ARCS::Graph& gin, DistGraph& gout) {
    // Add the vertices.
    for (const auto& it : contigToLength) {
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
    assert(!path.empty());
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
    assert(!tsvFile.empty());

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
    assert(!tsvFile.empty());

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
    std::cout << "\n\t=> Measuring intra-contig distances / shared barcodes... "
        << ctime(&rawtime);
    DistSampleMap distSamples;
    calcDistSamples(imap, contigToLength, indexMultMap, params, distSamples);

    time(&rawtime);
    std::cout << "\n\t=> Writing intra-contig distance samples to TSV... "
        << ctime(&rawtime);
    writeDistSamplesTSV(params.dist_samples_tsv, distSamples);

    time(&rawtime);
    std::cout << "\n\t=> Building Jaccard to distance map... "
        << ctime(&rawtime);
    JaccardToDist jaccardToDist;
    buildJaccardToDist(distSamples, jaccardToDist);

    time(&rawtime);
    std::cout << "\n\t=> Calculating barcode stats for scaffold pairs... "
        << ctime(&rawtime);
    PairToBarcodeStats pairToStats;
    buildPairToBarcodeStats(imap, indexMultMap, contigToLength, params, pairToStats);

    time(&rawtime);
    std::cout << "\n\t=> Adding edge distances... " << ctime(&rawtime);
    addEdgeDistances(pairToStats, jaccardToDist, params, g);

    if (!params.dist_tsv.empty()) {
        time(&rawtime);
        std::cout << "\n\t=> Writing distance estimates to TSV... "
            << ctime(&rawtime);
        writeDistTSV(params.dist_tsv, pairToStats, g);
    }
}

/** Run ARCS. */
void runArcs(const std::vector<std::string>& filenames) {

    std::cout << "Running: " << PROGRAM << " " << PACKAGE_VERSION
        << "\n this is add_arks branch and under modification"
        << "\n k-mer method: " << params.kmer_method
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

    ARCS::IndexMap imap;        // --same
    ARCS::PairMap pmap;         // --same except int<->unsigned difference
    ARCS::Graph g;              // --same
    std::unordered_map<std::string, int> indexMultMap;  // --same

    ARCS::ContigKMap kmap;      // --newly added
    kmap.set_deleted_key("");   // --newly added

    std::time_t rawtime;

    ARCS::ContigToLength contigToLength;    // --newly added normally it is declared in calcDistanceEstimates() function
    ARCS::ContigToLengthIt contigToLengthIt;

    ARCS::ScaffSizeList scaffSizeList;      // will be deleted at all.
    ARCS::ScaffSizeMap scaffSizeMap;        // will be converted to ContigToLength map structure.
    
    std::vector<ARCS::CI> contigRecord;

    if (!params.file.empty() && !params.kmer_method) {     // this is arcs method file type. scaffold file can be empty because same info can be gathered from SAM header file.(scaff id, length).
        time(&rawtime);
        std::cout << "\n=> Getting scaffold sizes... " << ctime(&rawtime);
        getScaffSizes(params.file, scaffSizeList, contigToLength);      // gets scaff ids and their sizes. Its said to be optional we'll see..
        scaffSizeMap.insert(scaffSizeList.begin(), scaffSizeList.end());
    }
    // -- results of above test are contigToLength works same as scaffSizeList after the methods called above. So no need for that.
    std::cout << "\n CTL Map size: " << contigToLength.size() << '\n';
    std::cout << "\n SSM Map size: " << scaffSizeMap.size() << '\n';
    
    // std::cout << "\n CTL Map key:468 value:" << contigToLength["10335"] << '\n';
    // std::cout << "\n SSM Map key:468 value:" << scaffSizeMap["10335"] << '\n';
    

    

    //std::unordered_map<std::string, int> indexMultMap;
    time(&rawtime);
    std::cout << "\n=> Reading alignment files... " << ctime(&rawtime);
    //std::vector<std::string> bamFiles = readFof(params.fofName);
    std::vector<std::string> fofFiles = readFof(params.fofName);    // readFof is generic for both arcs/arks.
    //std::copy(filenames.begin(), filenames.end(), std::back_inserter(bamFiles));
    std::copy(filenames.begin(), filenames.end(), std::back_inserter(fofFiles));
    std::cout << "\n Filenames have these files: \n";
    for(unsigned i = 0; i < fofFiles.size(); i++){
        std::cout << fofFiles.at(i) << '\n';
    }
    

    //readBAMS(bamFiles, imap, indexMultMap, scaffSizeList, scaffSizeMap);
    // this is a critical method as imap is constructed here, indexMultMap(barcodeMultiplicity map) constructed here. If scaffold file
    // is not given, then scaffSizeList and scaffSizeMap are also constructed here.
    if (!params.kmer_method)
    {
        readBAMS(fofFiles, imap, indexMultMap, scaffSizeList, scaffSizeMap, contigToLength);
    }else
    {
        std::cout << "\n----Full ARKS----\n" << std::endl;
        // if(!params.conrecfile.empty()){
        //     time(&rawtime);
		//     std::cout << "\n=>Detected ContigRecord file, making ContigRecord from checkpoint...\n" << ctime(&rawtime) << std::endl;
		//     createContigRecord(params.conrecfile, contigRecord);
        // }
        // if(!params.kmapfile.empty()){
		//     time(&rawtime);
		//     std::cout << "\n=>Detected ContigKmerMap file, making ContigKmerMap from checkpoint...\n" << ctime(&rawtime) << std::endl;
		//     createContigKmerMap(params.kmapfile, kmap);
        // }

        time(&rawtime);
        std::cout << "\n=>Preprocessing: Gathering barcode multiplicity information..." << ctime(&rawtime);
        createIndexMultMap(params.multfile, indexMultMap);
        
        std::cout << "barcode id: TACCACCCAAGGTCGA " << indexMultMap["TACCACCCAAGGTCGA-1"] << std::endl;

        time(&rawtime);
        std::cout << "\n=>Preprocessing: Gathering draft information..." << ctime(&rawtime) << "\n";
        size_t size = initContigArray(params.file); // why to read all of the file here for getting only the size.
        //std::vector<ARKS::CI> contigRecord(size);
        contigRecord.resize(size);

        //return;     // --- return border was set here and it worked!

    	time(&rawtime);
    	std::cout << "\n=>Storing Kmers from Contig ends... " << ctime(&rawtime) << std::endl;
    	getContigKmers(params.file, kmap, contigRecord, contigToLength);
        
        std::cout << "contigToLength cid:9 " << contigToLength["9"] << std::endl; //92
        std::cout << "contigToLength cid:10 " << contigToLength["10"] << std::endl; //674

        //return;     // --- return border is set here and working

        time(&rawtime);
  	    std::cout << "\n=>Reading Chromium FASTQ file(s)... " << ctime(&rawtime) << std::endl;
  	    //readChroms(inputFiles, kmap, imap, indexMultMap, contigRecord);
        readChroms(fofFiles, kmap, imap, indexMultMap, contigRecord);

        //return;

  	    std::cout << "Cumulative memory usage: " << memory_usage() << std::endl;
    }
    


    //return;     // --- return border is set here!
    //readBAMS(fofFiles, imap, indexMultMap, scaffSizeList, scaffSizeMap);
    



    size_t barcodeCount = countBarcodes(imap, indexMultMap);    // just for counting

    if (!params.barcode_counts_name.empty()) {
        time(&rawtime);
        std::cout << "\n=> Writing reads per barcode TSV file... " << ctime(&rawtime) << "\n";
        writeBarcodeCountsTSV(params.barcode_counts_name, indexMultMap);
    }

    time(&rawtime);
    std::cout << "\n=> Pairing scaffolds... " << ctime(&rawtime);
    pairContigs(imap, pmap, indexMultMap);

    time(&rawtime);
    std::cout << "\n=> Creating the graph... " << ctime(&rawtime);
    createGraph(pmap, g);

    if (params.dist_est) {
        std::cout << "\n=> Calculating distance estimates... " << ctime(&rawtime);
        calcDistanceEstimates(imap, indexMultMap, scaffSizeMap, g);
    }

    if (!params.base_name.empty()) {
        time(&rawtime);
        std::cout << "\n=> Writing graph file... " << ctime(&rawtime) << "\n";
        std::string graphFile = params.base_name + "_original.gv";
        writePostRemovalGraph(g, graphFile);
    }

    if (!params.dist_graph_name.empty()) {
        time(&rawtime);
        std::cout << "\n=> Creating the ABySS graph... " << ctime(&rawtime);
        DistGraph gdist;
        //createAbyssGraph(scaffSizeList, g, gdist);
        createAbyssGraph(contigToLength, g, gdist);

        time(&rawtime);
        std::cout << "\n=> Writing the ABySS graph file... " << ctime(&rawtime) << "\n";
        writeAbyssGraph(params.dist_graph_name, gdist);
    }

    if (!params.tsv_name.empty()) {
        time(&rawtime);
        std::cout << "\n=> Writing TSV file... " << ctime(&rawtime) << "\n";
        writeTSV(params.tsv_name, imap, pmap, barcodeCount);
    }

    time(&rawtime);
    std::cout << "\n=> Done. " << ctime(&rawtime);
}

int main(int argc, char** argv)
{
    
    printf("Reading user inputs...\n");
    // Couldn't understand how rawInputFiles string works.
	std::string rawInputFiles = "";

    // namesspace opt is declared in Graph/Options.h and used by DotIO.
    
    // trimMasked is not used any where in this script.(used only in arcs)
    opt::trimMasked = false;
    // bool die used in both.
    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            //*
            // Options only arks have are:
            
            // @ full, align and graph options are presented to user in arks with program param.
            // case 'p':
			//    arg >> params.program;
			//    break;

            // @ a parameter is used for multiplicity file declaration in arks. Can be tsv or gsv but required to create indexMultMap.
            // case 'a':
			//     arg >> params.multfile;
			//     break;

            // @ from this file contig record id(conreci), contig name and head/tail is read and vector is created in methods from it.
            // Not needed in full pipeline as it can be built in getContigKmers() method.
            // case 'q':
			//     arg >> params.conrecfile;
			//     break;

            // @ Gets a file to create kmap hash table. Not needed in full program.
            // case 'w':
			//     arg >> params.kmapfile;
			//     break;

            // @ Used only in graph program flow and specifies the input of index mappings and used to create imap.
            //case 'i':
			//    arg >> params.imapfile;
			//    break;

            // @ can get {0,1,2,3} which decides which of the checkpoint tsv output files will be created.
            // case 'o':
			//     arg >> params.checkpoint_outs;
			//     break;

            // @ size of k-mer.
            // case 'k':
			//     arg >> params.k_value;
			//     break;
		    
            // kmer shift step size.
            // case 'g':
			//     arg >> params.k_shift;
			//     break;

            // @ jaccard threshold used in determining the best contig according to kmer hits.
            // case 'j':
			//    arg >> params.j_index;
			//    break;

            // @ number fo threads running in the parallelized parts of chromiumRead() and atomic precautions taken in its 'calling' function
            // bestContig().
            // case 't':
			//     arg >> params.threads;
			//     break;
            

            // case 's':
			//     arg >> params.intra_contig_tsv;
			//     break;
		    
            // @ written distance/barcode data in TSV as output.
            // case 'S':
			//     arg >> params.inter_contig_tsv;
			//     break;

            case 'p':
			    arg >> params.program; break;
            case 'u':
			    arg >> params.multfile; break;
            case 'q':
			    arg >> params.conrecfile; break;
            case 'w':
			    arg >> params.kmapfile; break;
            case 'i':
		        arg >> params.imapfile; break;
            case 'o':
			    arg >> params.checkpoint_outs; break;
            case 'k':
			    arg >> params.k_value; break;
            case 'h':
		        arg >> params.k_shift; break;
            case 'j':
			    arg >> params.j_index; break;
            case 't':
			    arg >> params.threads; break;
            case 'x':
                arg >> params.kmer_method; break;   // -- newly added if kmer_method is 1 kmer arks method will be used.
                                                    // --- newly added params are above.

            case '?':
                die = true; break; //same
            case 'f':
                arg >> params.file; break; // same
            case 'a':
                arg >> params.fofName; break;   // -a is multfile in arks which is barcode multiplicity file having barcode and occurences(?).
                                                // fof file will be implemented to arks as well for chromium reads.
            case 'B':
                arg >> params.dist_bin_size; break; //same
            case 's':
                arg >> params.seq_id; break; // -s is intra_contig_tsv in arks which is a output, seq_id min sequence identity for read alignments in arcs.
            case 'c':
                arg >> params.min_reads; break; // same
            case 'D':
                params.dist_est = true; break; // same
            case 'l':
                arg >> params.min_links; break; // same
            case 'z':
                arg >> params.min_size; break; // same
            case 'b':
                arg >> params.base_name; break; // same --output file prefix.
            case 'g':
                arg >> params.dist_graph_name; break; // -g is k_shift in arks and in arcs an abyss graph is created additionally.
            case OPT_TSV:
                arg >> params.tsv_name; break;  // apart from a graph arcs outputs a detailed tsv file about the barcodes/scaffolds/links... where arks only
                                                // create checkpoint tsv files only. 
            case OPT_GAP:
                arg >> params.gap; break;       // fixed gap size for ABySS dist.gv file, naturally arks doesn't have it.
            case OPT_BARCODE_COUNTS:
                arg >> params.barcode_counts_name; break; // reads per barcode tsv file.
            case OPT_SAMPLES_TSV:
                arg >> params.dist_samples_tsv; break; // only name difference, intra_contig_tsv.
            case OPT_DIST_TSV:
                arg >> params.dist_tsv; break;  // only name difference, inter_contig_tsv.
            case OPT_NO_DIST_EST:
                params.dist_est = false; break; // same
            case OPT_DIST_MEDIAN:
                params.dist_mode = ARCS::DIST_MEDIAN; break;    // this is distance estimation method option not included in arks
            case OPT_DIST_UPPER:
                params.dist_mode = ARCS::DIST_UPPER; break;     // second option of the upper one.
            case 'm': {                                         // same
                std::string firstStr, secondStr;
                std::getline(arg, firstStr, '-');
                std::getline(arg, secondStr);
                std::stringstream ss;
                ss << firstStr << "\t" << secondStr;
                ss >> params.min_mult >> params.max_mult;
                }
                break;
            case 'd':
                arg >> params.max_degree; break;    // same
            case 'e':
                arg >> params.end_length; break;    // same
            case 'r':
                arg >> params.error_percent; break; // same
            case 'v':
                ++params.verbose; break;            // same
            case OPT_HELP:                          // same after here -->
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

    std::vector<std::string> filenames(argv + optind, argv + argc); // now filenames have specified SAM/BAM files or chromuim reads. 
    // later edit the readfof function that is compatible with reading of chromium reads.
    if (params.fofName.empty() && filenames.empty()) {      //user whether needs to put SAM/BAM files as parameter without flag or in a fof file.
            cerr << PROGRAM ": error: specify input (SAM/BAM file(s) || chromium reads) or a list of files with -a option\n";
            die = true;
    }

    std::ifstream g(params.file.c_str());       // this part can not be used for both as file is optional for arcs.
	if (!g.good() && params.kmer_method) {
		std::cerr << "Cannot find -f [required for kmer]" << params.file << ". Exiting... \n";
		die = true;
	}

    if (params.dist_graph_name.empty() && !params.base_name.empty())
            params.dist_graph_name = params.base_name + ".dist.gv";
    
    if (params.base_name.empty()                // think about this later whether it is necessary or not.
        && params.dist_graph_name.empty()
        && params.tsv_name.empty()) {
            cerr << PROGRAM ": error: specify an output file using one or more of -b, -g, or --tsv\n";
            die = true;
    }

    if (!params.file.empty())
            assert_readable(params.file);
    if (!params.fofName.empty())
        assert_readable(params.fofName);
    for (const auto& filename : filenames)
        assert_readable(filename);

    //std::vector<std::string> inputFiles;

    // this is not a good implementation of parameter processing in this case. But improving this part is left as later work.
    if(params.kmer_method){
        omp_set_num_threads(params.threads);    //omp method called only if kmer_method is using as openmp only implemented in arks.
        //gather all the chromium files
	    //inputFiles = convertInputString(rawInputFiles); // --convertInputString() is copied from arks.
                                                        // those elements in input file is just filename for chromium reads that iz opened by gzopen.

        // if (optind == argc) {   // optind is the number of arguments processed above. So the remaining must be chromium files.
		//     std::cerr << "No Chromium read files are specified. Exiting... \n"; 
		//     die=true;
	    // } else {
		//     while (optind < argc) {
		// 	    inputFiles.push_back(argv[optind]);
		// 	    optind++;
		//     }
	    // }
        // std::ifstream g(params.file.c_str());       // this part can not be used for both as file is optional for arcs.
	    // if (!g.good()) {
		//     std::cerr << "Cannot find -f " << params.file << ". Exiting... \n";
		//     die = true;
	    // }

	    if (params.program == "full") {         //those booleans specified here, are declared on the top of the script as global variables.
		    full = true;
	    } else if (params.program == "align") {
		    alignc = true;
	    } else if (params.program == "graph") {
		    graph = true;
	    }else {
		    std::cerr << "You must specify where you want ARKS(kmer) to start. Exiting... \n";
		    die = true;
	    }

	    // if (die) {              // this part can be used for both
		//     std::cerr << "Try " << PROGRAM << " --help for more information.\n";
		//     exit(EXIT_FAILURE);
	    // }   

	    /* Setting base name if not previously set */
	    if (params.base_name.empty()) {         // check here if you can add components by easily checking a boolean
		    std::ostringstream filename;
		    filename << params.file << ".scaff" << "k-method" << "_c"
				<< params.min_reads << "_k" << params.k_value << "_g"
				<< params.k_shift << "_j" << params.j_index << "_l"
				<< params.min_links << "_d" << params.max_degree << "_e"
				<< params.end_length << "_r" << params.error_percent;
		    params.base_name = filename.str();
	    }

	    // printf("%s\n", "Finished reading user inputs...entering runArks()...");

	    //runArks(inputFiles);
    }else{
        // Set base name if not previously set.
        if (params.base_name.empty() && !params.file.empty()) {
            std::ostringstream filename;
            filename << params.file << ".scaff" << "align-method"
                << "_s" << params.seq_id
                << "_c" << params.min_reads
                << "_l" << params.min_links
                << "_d" << params.max_degree
                << "_e" << params.end_length
                << "_r" << params.error_percent;
            params.base_name = filename.str();
        }

        // Set dist graph name unless specified.
        // if (params.dist_graph_name.empty() && !params.base_name.empty())
        //     params.dist_graph_name = params.base_name + ".dist.gv";

        // if (params.base_name.empty()
        // && params.dist_graph_name.empty()
        // && params.tsv_name.empty()) {
        //     cerr << PROGRAM ": error: specify an output file using one or more of -b, -g, or --tsv\n";
        //     die = true;
        // }

        //filenames(argv + optind, argv + argc);
        // if (params.fofName.empty() && filenames.empty()) {      //user whether needs to put SAM/BAM files as parameter without flag or in a fof file.
        //     cerr << PROGRAM ": error: specify input SAM/BAM file(s) or a list of files with -a option\n";
        //     die = true;
        // }

        // if (die) {
        //     std::cerr << "Try " PROGRAM " --help for more information.\n";
        //     exit(EXIT_FAILURE);
        // }

        // if (!params.file.empty())
        //     assert_readable(params.file);
        // if (!params.fofName.empty())
        //     assert_readable(params.fofName);
        // for (const auto& filename : filenames)
        //     assert_readable(filename);

        //runArcs(filenames);
    }

    
    if (die) {              // this part can be used for both
		    std::cerr << "Try " << PROGRAM << " --help for more information.\n";
		    exit(EXIT_FAILURE);
	}  

    printf("%s\n", "Finished reading user inputs...entering runArqs()...");
    runArcs(filenames);

    return 0;
}


//*
// Changes being made.
// Notes: Current script is arcs and being changed to be able to work with arks inputs and mechanism.
// 1. Starting from main function differences, mainly on parameters are investigated. Same and different parameters are idetified. Same parameters will
// remain unchanged. Additional parameters for arks will be added directly if the same parameter flag is not being used.
// These parameter changes are put in main function and param struct in header file but not in other methods like 'printing'. 

// The flags that are already in use are listed here: (---flag already taken)
// -a barcodeMultiplicity file in arks -a fofName in arcs. So it will be added as -u.
// -g k_shift file in arks -g dist_graph_name in arcs which is an abyss graph created additionally. So it will be added as -h.

// (---flag empty)
// -p program flow options for kmer mapping method. {full,align,graph}.
// -q conrecfile readed to create contig records not included in full flow.
// -w kmap input file not included in full flow.
// -i imap file
// -o checkpoint_out options {0,1,2,3}
// -k kvalue for kmer mapping
// -j threshold used in determining the best contig according to kmer hits. (j_index)
// -t thread number

// (---same flag same function different variable)
// -OPT_SAMPLES_TSV dist_samples_tsv(arcs) -> intra_contig_tsv(arks)
// OPT_DIST_TSV dist_tsv(arcs) -> inter_contig_tsv(arks)
// for those upper options arcs editions will be used for ease for now.

// BIG WARNING: Libraries didn't changed to be compatible as well.


// Questions to Warren a distance graph is outputted as called ABySS graph should I put it in new version? Example can be found on 
// /projects/btl_scratch/tgoktas/myDemo/k80_tigmint/arcs_results/elegans_outputs-scaffolds.tigmint_c5_m50-10000_s98_r0.05_e30000_z3000.dist.gv
// the main output graph of arks can be found at /projects/btl_scratch/tgoktas/newArks/arks/Examples/arks_test-demo/test_scaffolds_c5_m50-6000_k30_r0.05_e30000_z500_original.gv
// Firstly arks part will be merged as it is in full program.