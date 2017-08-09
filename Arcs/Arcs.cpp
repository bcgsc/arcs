#include "Arcs.h"
// for the inputted Fasta or Fastq file of contigs (if it is compressed with gzip)
#include <zlib.h>
#include "kseq.h"
#include <cassert>
#include <omp.h>

KSEQ_INIT(gzFile, gzread)

#define PROGRAM "arcs"
#define VERSION "1.0.1"

static const char VERSION_MESSAGE[] =
		"VERSION: " PROGRAM " " VERSION "\n"
		"\n"
		"http://www.bcgsc.ca/platform/bioinfo/software/links \n"
		"We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca.\n"
		"If you use LINKS, ARCS code or ideas, please cite our work. \n"
		"\n"
		"LINKS and ARCS Copyright (c) 2014-2016 Canada's Michael Smith Genome Science Centre.  All rights reserved. \n";

static const char USAGE_MESSAGE[] =
		"Usage: [" PROGRAM " " VERSION "]\n"
		"    arcs [Options] <chrom file list>\n"
		"Options:\n"
		"	=> TYPE OPTIONS: <=\n"
		"	   -p can be one of:\n"
		"			1) full    uses the full ARKS process (kmerize draft, kmerize and align chromium reads, scaffold).\n"
		"			2) align   skips kmerizing of draft and starts with kmerizing and aligning chromium reads.\n"
		"			3) graph   skips kmerizing draft and kmerizing/aligning chromium reads and only scaffolds.\n"
		"	=> INPUT OPTIONS: <=\n"
		"	    A) Always required (specific type 'full'):\n"
		"   		-f  Using kseq parser, these are the contig sequences to further scaffold and can be in either FASTA or FASTQ format. (required)\n"
		"   		-a  tsv or csv file for barcode multiplicities. Can be acquired by running CalcBarcodeMultiplicity script included. (required)\n"
		"	    B) If you want to skip kmerizing contigs you will need (specified type 'align'):\n" 
		"   		-q  tsv file for ContigRecord (a record of all the contigs + h/t). \n"
		"			--> Format of file should be: <contig record index number> <contig name> <H/T>\n"
		"   		-w  tsv file for the ContigKmerMap (a record of all the kmers to contig ID index (corresponds to contigrecord tsv).\n"
		"			--> Format of file should be: <kmer> <contig record index number>\n"
		"	    C) If you want to skip both the full kmer alignment based step you will need (specific type 'graph'):\n" 
		"			**Note that you are using ARKS as a graphing application**\n"
		"		-i  tsv file for the IndexMap.\n"
		"			--> Format of file should be: <barcode> <contig name> <H/T> <count>\n"
		"	=> EXTRA OUTPUT OPTIONS: <= \n"
		"		-o   can be one of: \n"
		"			0    no checkpoint files (default)\n"
		"			1    outputs of kmerizing the draft (ContigRecord + ContigKmerMap only)\n"
		"			2    output of aligning chromium to draft (IndexMap only)\n"
		"			3    all checkpoint files (ContigRecord, ContigKmerMap, and IndexMap)\n"
		"   -c  Minimum number of mapping read pairs/Index required before creating edge in graph. (default: 5)\n"
		"   -k  k-value for the size of a k-mer. (default: 30) (required)\n"
		"   -g  shift between k-mers (default: 1)\n"
		"   -j  Minimum Jaccard Index for a read to be associated with a contigId. (default: 0.55)\n"
		"   -l  Minimum number of links to create edge in graph (default: 0)\n"
		"   -z  Minimum contig length to consider for scaffolding (default: 500)\n"
		"   -b  Base name for your output files (optional)\n"
		"   -m  Range (in the format min-max) of index multiplicity (only reads with indices in this multiplicity range will be included in graph) (default: 50-10000)\n"
		"   -d  Maximum degree of nodes in graph. All nodes with degree greater than this number will be removed from the graph prior to printing final graph. For no node removal, set to 0 (default: 0)\n"
		"   -e  End length (bp) of sequences to consider (default: 30000)\n"
		"   -r  Maximum p-value for H/T assignment and link orientation determination. Lower is more stringent (default: 0.05)\n"
		"   -t 	Number of threads.(default: 1)\n"
		"   -v  Runs in verbose mode (optional, default: 0)\n";

/* ARCS PREPARATION AKA GLOBAL VARIABLES: */

ARCS::ArcsParams params;

static const char shortopts[] = "p:f:a:q:w:i:o:c:k:g:j:l:z:b:m:d:e:r:vt:";

enum { OPT_HELP = 1, OPT_VERSION};

static const struct option longopts[] = {
    {"program", required_argument, NULL, 'p'}, 
    {"file", required_argument, NULL, 'f'},
    {"multfile", required_argument, NULL, 'a'}, 
    {"conrecfile", required_argument, NULL, 'q'}, 
    {"kmapfile", required_argument, NULL, 'w'}, 
    {"imapfile", required_argument, NULL, 'i'},
    {"checkpoint_outs", required_argument, NULL, 'o'}, 
    {"min_reads", required_argument, NULL, 'c'},
    {"k_value", required_argument, NULL, 'k'}, 
    {"k_shift", required_argument, NULL, 'g'}, 
    {"j_index", required_argument, NULL, 'j'}, 
    {"min_links", required_argument, NULL, 'l'},
    {"min_size", required_argument, NULL, 'z'},
    {"base_name", required_argument, NULL, 'b'},
    {"index_multiplicity", required_argument, NULL, 'm'},
    {"max_degree", required_argument, NULL, 'd'},
    {"end_length", required_argument, NULL, 'e'},
    {"error_percent", required_argument, NULL, 'r'},
    {"run_verbose", no_argument, NULL, 'v'},
    {"threads", required_argument, NULL, 't'},
    {"version", no_argument, NULL, OPT_VERSION},
    {"help", no_argument, NULL, OPT_HELP},
    { NULL, 0, NULL, 0 }
};

//TODO: PLEASE ADD THE s_ prefex (eg. s_numkmercollisions) to denote that these are static variables during the runtime of the program

unsigned int numkmersmapped = 0, numkmercollisions = 0, numkmersremdup = 0,
		numbadkmers = 0, uniquedraftkmers = 0;

unsigned int totalnumckmers = 0, ckmersasdups = 0, numckmersfound = 0,
		numckmersrec = 0, numbadckmers = 0;

unsigned int numreadspassingjaccard = 0, numreadsfailjaccard = 0;

bool full = false; 
bool alignc = false; 
bool graph = false;

/* HELPERS FOR CHECKING AND PRINTING: */

std::string HeadOrTail(bool orientation) {
	if (orientation) {
		return "H";
	} else {
		return "T";
	}
}

/* Precondition: input must either be "H" or "T" */
bool HTtoBool(std::string ht) {
	if (ht == "H") {
		return true; 
	} else {
		return false; 
	}
}

/* Returns true if the barcode only contains ATGC */
static inline bool checkIndex(std::string seq) {
	for (int i = 0; i < static_cast<int>(seq.length()); i++) {
		char c = toupper(seq[i]);
		if (c != 'A' && c != 'T' && c != 'G' && c != 'C')
			return false;

	}
	//return (static_cast<int>(seq.length()) == params.indexLen);
	return true;
}

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

//for multiple fastq chromium files
vector<string> convertInputString(const string &inputString) {
	vector<string> currentInfoFile; 
	string temp; 
	stringstream converter(inputString); 
	while (converter >> temp) {
		currentInfoFile.push_back(temp); 
	}
	return currentInfoFile; 
}

/* Write TSV checkpoint files */

/* writes contigRecord table to TSV */
void writeContigRecord(std::vector<ARCS::CI> contigRecord) {

	std::string outputfilename = params.base_name + "_contigrec.tsv"; 

	FILE* fout = fopen(outputfilename.c_str(), "w"); 

	size_t conreci = 0; 
	for (std::vector<ARCS::CI>::iterator it = contigRecord.begin(); it != contigRecord.end(); ++it) {
		std::string contigname = it->first; 
		std::string ht = HeadOrTail(it->second); 
		fprintf(fout, "%zu\t%s\t%s\n", conreci, contigname.c_str(), ht.c_str()); 
		conreci++; 
	}
	fclose(fout); 
} 

/* writes contigKMap to TSV */
void writeContigKmerMap(ARCS::ContigKMap kmap) {
	
	std::string outputfilename = params.base_name + "_kmercontigrec.tsv"; 

	FILE* fout = fopen(outputfilename.c_str(), "w"); 

	for (auto it=kmap.begin(); it != kmap.end(); ++it) {
		std::string kmer = it->first; 
		size_t contigreci = it->second; 
		fprintf(fout, "%s\t%zu\n", kmer.c_str(), contigreci); 
	}
	fclose(fout); 
}

/* write IndexMap to TSV */
void writeIndexMap(ARCS::IndexMap imap) {

	std::string barcode = ""; 
	std::string contigname = ""; 
	std::string orientation = ""; 
	int count = 0; 

	std::string outputfilename = params.base_name + "_imap.tsv"; 

	FILE* fout = fopen(outputfilename.c_str(), "w"); 

	for (auto it = imap.begin(); it != imap.end(); ++it) {
		barcode = it->first; 
		ARCS::ScafMap smap = it->second; 
		for (auto j = smap.begin(); j != smap.end(); ++j) {
			contigname = j->first.first; 
			orientation = HeadOrTail(j->first.second); 
			count = j->second; 
			fprintf(fout, "%s\t%s\t%s\t%d\n", barcode.c_str(), contigname.c_str(), orientation.c_str(), count); 
		}
	}
	fclose(fout); 
}
			
/* Reading TSV (or CSV for barcode mult) checkpoint files */			

/* Create indexMultMap from Barcode Multiplicity File */
void createIndexMultMap(std::string multfile, std::unordered_map<std::string, int> &indexMultMap) {

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

		std::stringstream sst(line); 

		std::string barcode; 
		std::string multiplicity_string;
		int multiplicity; 

		if (tsv) { 
			sst >> barcode >> multiplicity_string; 
			multiplicity = std::stoi(multiplicity_string);  
		} else {
			getline(getline(sst, barcode, ','), multiplicity_string, ','); 
			multiplicity = std::stoi(multiplicity_string); 
		}

		if (checkIndex(barcode)) {
			indexMultMap[barcode] = multiplicity; 
		} else {
			std::cout << "Please check your multiplicity file." << std::endl; 
		}

		assert(multfile_stream); 
	}	
	multfile_stream.close(); 		
}

/* Create contigRecord vector */
void createContigRecord(std::string contigrectsv, std::vector<ARCS::CI> contigRecord) {
	
	std::ifstream contigrectsv_stream; 
	contigrectsv_stream.open(contigrectsv.c_str()); 
	if (!contigrectsv_stream) {
		std::cerr << "Could not open " << contigrectsv << ". --fatal.\n"; 
		exit (EXIT_FAILURE); 
	}

	std::string line; 
	while (getline(contigrectsv_stream, line)) {
		std::stringstream sst(line); 

		int contigreci;
		bool ht;
		std::string contigreci_string, contigname, headortail; 
		sst >> contigreci_string >> contigname >> headortail; 
		contigreci = std::stoi(contigreci_string); 
		ht = HTtoBool(headortail); 

		ARCS::CI contigID(contigname, ht); 
	
		contigRecord[contigreci] = contigID; 
	}
	contigrectsv_stream.close(); 
}

/* Create ContigKmerMap */
void createContigKmerMap(std::string kmaptsv, ARCS::ContigKMap kmap) {

	std::ifstream kmaptsv_stream; 
	kmaptsv_stream.open(kmaptsv.c_str()); 
	if (!kmaptsv_stream) {
		std::cerr << "Could not open " << kmaptsv << ". --fatal.\n"; 
		exit (EXIT_FAILURE); 
	} 

	std::string line; 
	while (getline(kmaptsv_stream, line)) {
		std::stringstream sst(line); 

		std::string kmer, contigreci_string; 
		int contigreci; 

		sst >> kmer >> contigreci_string; 
		contigreci = std::stoi(contigreci_string); 

		kmap[kmer] = contigreci; 
	}
	kmaptsv_stream.close(); 
}

/* Create IndexMap */
void createIndexMap(std::string imaptsv, ARCS::IndexMap imap) {

	std::ifstream imaptsv_stream; 
	imaptsv_stream.open(imaptsv.c_str()); 
	if (!imaptsv_stream) {
		std::cerr << "Could not open " << imaptsv << ". --fatal.\n"; 
		exit (EXIT_FAILURE); 
	}

	std::string line; 
	while (getline(imaptsv_stream, line)) {
		std::stringstream sst(line); 

		std::string barcode, contigname, ht_string, count_string; 
		bool ht; 
		int count; 
		
		sst >> barcode >> contigname >> ht_string >> count_string; 
		ht = HTtoBool(ht_string); 
		count = std::stoi(count_string); 

		ARCS::CI contigID(contigname, ht); 

		imap[barcode][contigID] = count; 
	}
	imaptsv_stream.close(); 
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

/* ARCS PROCESSES FUNCTIONS */

/* Returns the size of the array for storing contigs */
size_t initContigArray(std::string contigfile) {

	size_t count = 0;

	gzFile fp;

	int l;
	const char* filename = contigfile.c_str();
	fp = gzopen(filename, "r");
	kseq_t * seq = kseq_init(fp);

	while ((l = kseq_read(seq)) >= 0) {
		std::string sequence = seq->seq.s;
		unsigned sequence_length = sequence.length();
		if (checkContigSequence(sequence) && sequence_length >= params.min_size)
			count++;
	}
	kseq_destroy(seq);
	gzclose(fp);

	if (params.verbose) {
		cerr << "Number of contigs:" << count << "\nSize of Contig Array:"
				<< (count * 2) + 1 << endl;
	}

	//Let the first index represent the a null contig
	return (count * 2) + 1;
}

/* Shreds end sequence into kmers and inputs them one by one into the ContigKMap 
 * 	std::pair<std::string, bool> 				specifies contigID and head/tail 
 * 	std::string						the end sequence of the contig 
 *	int							k-value specified by user 
 *	ARCS::ContigKMap					ContigKMap for storage of kmers
 */
int mapKmers(std::string seqToKmerize, int k, int k_shift,
		ARCS::ContigKMap& kmap, ReadsProcessor &proc, int conreci) {
	/*
	 //For debugging purposes:
	 std::string kmaplist;*/

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
			const unsigned char* temp = proc.prepSeq(seqToKmerize, i); // prepSeq returns NULL if it contains an N
			// Ignore a NULL kmer
			if (temp != NULL) {
				std::string kmerseq = proc.getStr(temp);
				/*
				 kmaplist += proc.getBases(temp);
				 kmaplist += "\t\t\t";
				 kmaplist += proc.getBinary(temp);
				 kmaplist += "\t\t\t";
				 kmaplist += kmerseq;
				 */

				numKmers++;

				//Critical section to prevent multiple thread from accessing datastructure
				//only needed since we are writing to datastructure
				//ideally would want CAS operation but without special a datastructure this is difficult

#pragma omp critical(kmerseq)
				if (kmap.find(kmerseq) != kmap.end()) {
//#pragma omp atomic
					numkmercollisions++;
					if (kmap[kmerseq] != conreci) {
//#pragma omp atomic
						numkmersremdup++;
						if (kmap[kmerseq] != 0) {
							//kmaplist += "\t\t\tremoved"; 
//#pragma omp atomic
							uniquedraftkmers--;
							kmap[kmerseq] = 0;
						} else {
							kmap[kmerseq] = 0;
						}
					}
				} else {
					assert(kmap.count(kmerseq) == 0);
					kmap[kmerseq] = conreci;
//#pragma omp atomic
					uniquedraftkmers++;
//#pragma omp atomic
					numkmersmapped++;
				}
				i += k_shift;
			} else {
				i += k;
#pragma omp atomic
				numbadkmers++;
			}
			//kmaplist += "\n"; 
		}
		//writeToKmapDebugLogFile(kmaplist); 
		return numKmers;
	}
}

/* Get the k-mers from the paired ends of the contigs and store them in map. 
 * 	std::string file					FASTA (or later FASTQ) file 
 *	std::sparse_hash_map<k-mer, pair<contidID, bool>> 	ContigKMap
 *	int k							k-value (specified by user)
 */
void getContigKmers(std::string contigfile, ARCS::ContigKMap& kmap,
		std::vector<ARCS::CI> &contigRecord) {

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
	vector<ReadsProcessor*> procs(params.threads);

	for (unsigned i = 0; i < params.threads; ++i) {
		procs[i] = new ReadsProcessor(params.k_value);
	}

#pragma omp parallel
	while (l >= 0) {
		bool good = false;
		std::string contigID = "", sequence = "";
		size_t tempConreci1 = 0;
		size_t tempConreci2 = 0;
#pragma omp critical(kseq)
		{
			l = kseq_read(seq);
			if (l >= 0) {
				contigID = seq->name.s;
				sequence = seq->seq.s;
				if (sequence.length() >= params.min_size) {
					tempConreci1 = ++conreci;
					tempConreci2 = ++conreci;
					good = true;
				}
			}
		}
#pragma omp atomic
		totalNumContigs++;

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
//				if (sequence_length >= params.min_size) {

					// If contig length is less than 2 x end_length, then we split the sequence
					// in half to decide head/tail (aka we changed the end_length)
					int cutOff = params.end_length;
					if (cutOff == 0 || sequence_length <= cutOff * 2)
						cutOff = sequence_length / 2;

					// Arbitrarily assign head or tail to ends of the contig
					ARCS::CI headside(contigID, true);
					ARCS::CI tailside(contigID, false);

					//get ends of the sequence and put k-mers into the map
					contigRecord[tempConreci1] = headside;
					std::string seqend;
					seqend = sequence.substr(0, cutOff);
					int num = mapKmers(seqend, params.k_value, params.k_shift,
							kmap, *procs[omp_get_thread_num()], tempConreci1);
#pragma omp atomic
					totalKmers += num;

					contigRecord[tempConreci2] = tailside;
					seqend = sequence.substr(sequence_length - cutOff,
							sequence_length);
					num = mapKmers(seqend, params.k_value, params.k_shift, kmap,
							*procs[omp_get_thread_num()], tempConreci2);
#pragma omp atomic
					totalKmers += num;
#pragma omp atomic
					validContigs++;
		} else {
#pragma omp atomic
			skippedContigs++;
		}
//			}
//		}

		// printprogress
		if (params.verbose) {
#pragma omp critical(stdout)
			if (totalNumContigs % 1000 == 0) {
				printf("Finished %d Contigs...\n", totalNumContigs);
				//std::cout << "Cumulative memory usage: " << memory_usage()
				//		<< std::endl;
				//std::cout << "Kmers so far: " << numkmersmapped << std::endl;
			}
		}
	}
	kseq_destroy(seq);
	gzclose(fp);

	// clean up
	for (unsigned i = 0; i < params.threads; ++i) {
		delete procs[i];
	}

	if (params.verbose) {
		printf(
				"%s %u\n%s %u\n%s %u\n%s %u\n%s %u\n%s %u\n%s %u\n%s %u\n%s %u\n",
				"Total number of contigs in draft genome: ", totalNumContigs,
				"Total valid contigs: ", validContigs,
				"Total skipped contigs: ", skippedContigs,
				"Total number of Kmers: ", totalKmers, "Number Null Kmers: ",
				numbadkmers, "Number Kmers Recorded: ", numkmersmapped,
				"Number Kmer Collisions: ", numkmercollisions,
				"Number Times Kmers Removed (since duplicate in different contig): ",
				numkmersremdup, "Number of unique kmers (only one contig): ",
				uniquedraftkmers);
	}
}


/* Calculate the jaccard index */
static inline double calcJacIndex(int smallCount, int overallCount) {
	return (double) smallCount / (double) overallCount;
}

/* Returns best corresponding contig from read through kmers
 * 	ARCS::ContigKMap			tells me what kmers correspond to which contig
 *	std::string				read sequence
 *	int					size of k-mer
 *	int 					k_shift
 *      double j_index				Jaccard Index (default 0.55)
 *	ReadsProcessor				kmerizer
 */
int bestContig(ARCS::ContigKMap &kmap, std::string readseq, int k, int k_shift,
		double j_index, ReadsProcessor &proc) {

	// to keep track of what contig+H/T that the k-mer from barcode matches to
	// 	int					Index that corresponds to the contig in the contigRecord
	// 	count					# kmers found here
	std::map<int, int> ktrack;

	// k-merize readsequence
	int corrConReci = 0;
	int seqlen = readseq.length();

	int totalnumkmers = 0;
	int kmerdups = 0;
	int kmerfound = 0;
	int kmerstore = 0;
	
	 //For debugging purposes
	 //std::string ckmerlist;

	int i = 0;

	while (i <= seqlen - k) {
		const unsigned char* temp = proc.prepSeq(readseq, i);
		totalnumkmers++;
		if (temp != NULL) {
			std::string ckmerseq = proc.getStr(temp);
			totalnumckmers++;

#pragma omp critical(ckmerseq)
			// search for kmer in ContigKmerMap and only record if it is not the collisionmaker
			if (kmap.find(ckmerseq) != kmap.end()) {
				kmerfound++;
				numckmersfound++;
				corrConReci = kmap[ckmerseq];
				if (corrConReci != 0) {
					ktrack[corrConReci]++;
					kmerstore++;
					numckmersrec++;
				} else {
					ckmersasdups++;
					kmerdups++;
				}
			}
			i += k_shift; 
		} else {
			i += k_shift;
#pragma omp atomic 
			numbadckmers++;
		}
	}

	double maxjaccardindex = 0;
	// for the read, find the contig that it is most compatible with based on the jaccard index
	for (auto it = ktrack.begin(); it != ktrack.end(); ++it) {
		double currjaccardindex = calcJacIndex(it->second, totalnumkmers);
		if (maxjaccardindex < currjaccardindex) {
			maxjaccardindex = currjaccardindex;
			corrConReci = it->first;
		}
	}

	// default accuracythreshold is 0.55)
	if (maxjaccardindex > j_index) {
		numreadspassingjaccard++;
		return corrConReci;
	} else {
		numreadsfailjaccard++;
		return 0;
	}
}


/* Read through longranger basic chromium output fastq file */
void chromiumRead(std::string chromiumfile, ARCS::ContigKMap& kmap, ARCS::IndexMap& imap,
			std::unordered_map<std::string, int> &indexMultMap,
			std::vector<ARCS::CI> &contigRecord) {

	int stored_readpairs = 0;
	int skipped_unpaired = 0;
	int skipped_invalidbarcode = 0;
	int skipped_nogoodcontig = 0;

	int count = 0;

	//each thread gets a proc;
	vector<ReadsProcessor*> procs(params.threads);

	for (unsigned i = 0; i < params.threads; ++i) {
		procs[i] = new ReadsProcessor(params.k_value);
	}

	gzFile fp2;
	int l;
	const char* filename = chromiumfile.c_str();
	fp2 = gzopen(filename, "r");
	if (fp2 == Z_NULL) {
		cerr << "File " << filename << " cannot be opened." << endl;
		exit(1);
	} else {
		cerr << "File " << filename << " opened." << endl;
	}
	kseq_t * seq2 = kseq_init(fp2);

	
#pragma omp parallel 
	while (l>=0) {
		std::string read1_name = ""; 
		std::string read2_name = ""; 
		std::string barcode1 = ""; 
		std::string barcode2 = ""; 
		std::string cread1 = ""; 
		std::string cread2 = ""; 
		std::string comment1 = ""; 
		std::string comment2 = ""; 
		bool paired = false; 
		int indexMult = 0; 
		int corrConReci1 = 0;
		int corrConReci2 = 0; 
		ARCS::CI corrContigId;
#pragma omp atomic
		count++; 
#pragma omp critical(checkread1or2) 
		{	
			l = kseq_read(seq2);
			if (l >= 0) {
				read1_name = seq2->name.s; 
				if (seq2->comment.l) {
					comment1 = seq2->comment.s; 
				}
				cread1 = seq2->seq.s; 
			}
			l = kseq_read(seq2); 
			if (l >= 0) {
				read2_name = seq2->name.s; 
				if (seq2->comment.l) {
					comment2 = seq2->comment.s;
				}
				cread2 = seq2->seq.s;
			}
			if (read1_name == read2_name) {
				paired = true; 
			} else {
				std::cout << "File contains unpaired reads: " << read1_name << " " << read2_name << std::endl; 
				skipped_unpaired++; 		
			}
		}

		for (std::string::iterator i = comment1.begin(); i != comment1.end();
				i++) {
			if (*i != 'B' && *i != 'X' && *i != ':' && *i != 'Z'
					&& *i != '-' && *i != '1' && *i != '\n') {
				barcode1 += *i;
			}
		}
		for (std::string::iterator i = comment2.begin(); i != comment2.end();
				i++) {
			if (*i != 'B' && *i != 'X' && *i != ':' && *i != 'Z'
					&& *i != '-' && *i != '1' && *i != '\n') {
				barcode2 += *i;
			}
		}

		bool good = paired && !barcode1.empty() && !barcode2.empty() && (barcode1==barcode2); 
		if (good) {			
			indexMult = indexMultMap[barcode1]; 
			bool goodmult = indexMult > params.min_mult || indexMult < params.max_mult;
			if (goodmult) {
				corrConReci1 = bestContig(kmap, cread1, params.k_value, params.k_shift, params.j_index, *procs[omp_get_thread_num()]); 
				corrConReci2 = bestContig(kmap, cread2, params.k_value, params.k_shift, params.j_index, *procs[omp_get_thread_num()]);
			} else {
#pragma omp atomic
				skipped_invalidbarcode++; 
			}
			// we only store barcode info in index map if read pairs have same contig + orientation	
			// and if the corrContigId is not NULL (because it is above accuracy threshold)
			if (corrConReci1 != 0 && corrConReci2 != 0 && corrConReci1 == corrConReci2) {
#pragma omp critical(imap)
				{
					corrContigId = contigRecord[corrConReci1];
					imap[barcode1][corrContigId]++;
				}

#pragma omp atomic
				stored_readpairs++;

			} else {
#pragma omp atomic
				skipped_nogoodcontig++;
			}
		}

		if (params.verbose) {
#pragma omp critical(readstdout)
			if (count % 2000000 == 0) {
				std::cout << "Processed " << count << " read pairs." << std::endl; 
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
				"Stored read pairs: %u\nSkipped reads with invalid barcodes: %u\nSkipped unpaired reads: %u\nSkipped reads pairs without a good contig: %u\n",
				stored_readpairs, skipped_invalidbarcode, skipped_unpaired,
				skipped_nogoodcontig);
		printf(
				"Total valid kmers: %u\nNumber invalid kmers: %u\nNumber of kmers found in ContigKmap: %u\nNumber of kmers recorded in Ktrack: %u\nNumber of kmers found in ContigKmap but duplicate: %u\nNumber of reads passing jaccard threshold: %u\nNumber of reads failing jaccard threshold: %u\n",
				totalnumckmers, numbadckmers, numckmersfound, numckmersrec,
				ckmersasdups, numreadspassingjaccard, numreadsfailjaccard);
	}
} 


void readChroms(vector<string> inputFiles, ARCS::ContigKMap &kmap,
		ARCS::IndexMap &imap,
		std::unordered_map<std::string, int> &indexMultMap,
		std::vector<ARCS::CI> &contigRecord) {

	std::string chromFile;

	for (auto p = inputFiles.begin(); p != inputFiles.end(); ++p) {
		chromFile = *p; 
		if (params.verbose)
			std::cout << "Reading chrom " << chromFile << std::endl;
		chromiumRead(chromFile, kmap, imap, indexMultMap, contigRecord);
	}
}

		
/*
 * Check if SAM flag is one of the accepted ones.
 */
static inline bool checkFlag(int flag) {
	return (flag == 99 || flag == 163 || flag == 83 || flag == 147);
}

/*
 * Check if character is one of the accepted ones.
 */
static inline bool checkChar(char c) {
	return (c == 'M' || c == '=' || c == 'X' || c == 'I');
}

/* Normal approximation to the binomial distribution */
static inline float normalEstimation(int x, float p, int n) {
	float mean = n * p;
	float sd = std::sqrt(n * p * (1 - p));
	return 0.5 * (1 + std::erf((x - mean) / (sd * std::sqrt(2))));
}

/*
 * Based on number of read pairs that align to the 
 * head or tail of scaffold, determine if is significantly 
 * different from a uniform distribution (p=0.5)
 */
static inline std::pair<bool, bool> headOrTail(int head, int tail) {
	int max = std::max(head, tail);
	int sum = head + tail;
	if (sum < params.min_reads) {
		return std::pair<bool, bool>(false, false);
	}
	float normalCdf = normalEstimation(max, 0.5, sum);
	if (1 - normalCdf < params.error_percent) {
		bool isHead = (max == head);
		return std::pair<bool, bool>(true, isHead);
	} else {
		return std::pair<bool, bool>(false, false);
	}
}


/* 
 * Iterate through IndexMap and for every pair of scaffolds
 * that align to the same index, store in PairMap. PairMap 
 * is a map with a key of pairs of saffold names, and value
 * of number of links between the pair. (Each link is one index).
 */
void pairContigs(ARCS::IndexMap& imap, ARCS::PairMap& pmap,
		std::unordered_map<std::string, int>& indexMultMap) {

	/* Iterate through each index in IndexMap */
	for (auto it = imap.begin(); it != imap.end(); ++it) {

		/* Get index multiplicity from indexMultMap */
		std::string index = it->first;
		int indexMult = indexMultMap[index];

		if (indexMult >= params.min_mult && indexMult <= params.max_mult) {

			/* Iterate through all the scafNames in ScafMap */
			for (auto o = it->second.begin(); o != it->second.end(); ++o) {
				for (auto p = it->second.begin(); p != it->second.end(); ++p) {
					std::string scafA, scafB;
					bool scafAflag, scafBflag;
					std::tie(scafA, scafAflag) = o->first;
					std::tie(scafB, scafBflag) = p->first;

					/* Only insert into pmap if scafA < scafB to avoid duplicates */
					if (scafA != scafB) {
						bool validA, validB, scafAhead, scafBhead;

						std::tie(validA, scafAhead) = headOrTail(
								it->second[std::pair<std::string, bool>(scafA,
										true)],
								it->second[std::pair<std::string, bool>(scafA,
										false)]);
						std::tie(validB, scafBhead) = headOrTail(
								it->second[std::pair<std::string, bool>(scafB,
										true)],
								it->second[std::pair<std::string, bool>(scafB,
										false)]);

						if (validA && validB) {
							std::pair<std::string, std::string> pair(scafA,
									scafB);
							std::pair<std::string, std::string> reversepair(
									scafB, scafA);

							if (pmap.count(reversepair) != 0) {
								pair = reversepair;
								bool temp = scafAhead;
								scafAhead = scafBhead;
								scafBhead = temp;
							}

							if (pmap.count(pair) == 0) {
								std::vector<int> init(4, 0);
								pmap[pair] = init;
							}

							// Head - Head
							if (scafAhead && scafBhead) {
								pmap[pair][0]++;
								// Head - Tail
							} else if (scafAhead && !scafBhead) {
								pmap[pair][1]++;
								// Tail - Head
							} else if (!scafAhead && scafBhead) {
								pmap[pair][2]++;
								// Tail - Tail
							} else if (!scafAhead && !scafBhead) {
								pmap[pair][3]++;
							}
						}
					}
				}
			}
		}
	}
	//writePairMapToLog(pmap);
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
static inline bool checkSignificance(int max, int second) {
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
	for (it = pmap.begin(); it != pmap.end(); ++it) {
		std::string scaf1, scaf2;
		std::tie(scaf1, scaf2) = it->first;

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
			std::tie(e, inserted) = boost::add_edge(vmap[scaf1], vmap[scaf2],
					g);
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
		std::cout << "      Deleting nodes with degree > " << params.max_degree
				<< "... \n";
		removeDegreeNodes(g, params.max_degree);
	} else {
		std::cout << "      Max Degree (-d) set to: " << params.max_degree
				<< ". Will not delete any verticies from graph.\n";
	}

	std::cout << "      Writting graph file to " << graphFile << "...\n";
	writeGraph(graphFile, g);
}

void runArcs(vector<string> inputFiles) {
    std::cout << "Entered runArcs()..." << std::endl; 

    std::cout << "Running: " << PROGRAM << " " << VERSION 
        << "\n pid " << ::getpid()
	<< "\n -p " << params.program
        << "\n -f " << params.file
	<< "\n -a " << params.multfile
	<< "\n -q " << params.conrecfile
	<< "\n -w " << params.kmapfile
	<< "\n -i " << params.imapfile
	<< "\n -o " << params.checkpoint_outs
        << "\n -c " << params.min_reads 
	<< "\n -k " << params.k_value
	<< "\n -g " << params.k_shift
	<< "\n -j " << params.j_index
        << "\n -l " << params.min_links     
        << "\n -z " << params.min_size
        << "\n -b " << params.base_name
        << "\n Min index multiplicity: " << params.min_mult 
        << "\n Max index multiplicity: " << params.max_mult 
        << "\n -d " << params.max_degree 
        << "\n -e " << params.end_length
        << "\n -r " << params.error_percent
	<< "\n -t " << params.threads
        << "\n -v " << params.verbose << "\n";

    std::string graphFile = params.base_name + "_original.gv";

    ARCS::ContigKMap kmap; 
    kmap.set_deleted_key(""); 

    ARCS::IndexMap imap;
    ARCS::PairMap pmap;
    ARCS::Graph g;
    std::unordered_map<std::string, int> indexMultMap;

    std::time_t rawtime;

    std::cout << "\n---We are using KMER method.---\n" << std::endl; 	

    time(&rawtime); 
    std::cout << "\n=>Preprocessing: Gathering barcode multiplicity information..." << ctime(&rawtime); 
    createIndexMultMap(params.multfile, indexMultMap);

    time(&rawtime); 
    std::cout << "\n=>Preprocessing: Gathering draft information..." << ctime(&rawtime) << "\n"; 
    size_t size = initContigArray(params.file); 
    std::vector<ARCS::CI> contigRecord(size); 

    if (full) {
	assert(!graph);
	assert(!alignc);

	std::cout << "\n----Full ARKS----\n" << std::endl;  
		
    	// Read contig file, shred sequences into k-mers, and then map them 
    	time(&rawtime); 
    	std::cout << "\n=>Storing Kmers from Contig ends... " << ctime(&rawtime) << "\n"; 
    	getContigKmers(params.file, kmap, contigRecord);
    } 

    if (full || alignc) {
	  assert(!graph);

	  if (!full && alignc) {
		std::cout << "\n----Kmer Align ARKS----\n" << std::endl; 

		time(&rawtime); 
		std::cout << "\n=>Detected ContigRecord file, making ContigRecord from checkpoint...\n" << ctime(&rawtime) << "\n"; 
		createContigRecord(params.conrecfile, contigRecord); 

		time(&rawtime); 
		std::cout << "\n=>Detected ContigKmerMap file, making ContigKmerMap from checkpoint...\n" << ctime(&rawtime) << "\n";
		createContigKmerMap(params.kmapfile, kmap); 
	  } 

  	  time(&rawtime); 
  	  std::cout << "\n=>Reading Chromium FASTQ file(s)... " << ctime(&rawtime) << "\n"; 
  	  readChroms(inputFiles, kmap, imap, indexMultMap, contigRecord); 

  	  std::cout << "Cumulative memory usage: " << memory_usage() << std::endl;
    }

    if (graph) {
	assert(!full); 
	assert(!alignc); 

	std::cout << "\n----Graph ARKS----\n" << std::endl; 

	time(&rawtime); 
	std::cout << "\n=>Detected IndexMap file, making IndexMap from checkpoint...\n" << ctime(&rawtime) << "\n"; 
	createIndexMap(params.imapfile, imap); 
    }
	 
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
    std::cout << "\n=>Outputting desired checkpoint files... " << ctime(&rawtime) << "\n"; 
    int o = params.checkpoint_outs; 
    switch (o) {
	case 3:
		writeContigRecord(contigRecord); 
		writeContigKmerMap(kmap); 
		writeIndexMap(imap); 
		break; 
	case 2: 
		writeIndexMap(imap); 
		break; 
	case 1: 
		writeContigRecord(contigRecord); 
		writeContigKmerMap(kmap); 
		break; 
	case 0:
	default: 
		break; 
    }

    time(&rawtime);
    std::cout << "\n=>Done. " << ctime(&rawtime);
}

int main(int argc, char** argv) {

	printf("Reading user inputs...\n");

	std::string rawInputFiles = "";  

	bool die = false;
	for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;
			) {
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case 'p': 
			arg >> params.program; 
			break; 
		case '?':
			die = true;
			break;
		case 'f':
			arg >> params.file;
			break;
		case 'a': 
			arg >> params.multfile; 
			break; 
		case 'q':
			arg >> params.conrecfile; 
			break; 
		case 'w': 
			arg >> params.kmapfile; 
			break; 
		case 'i': 
			arg >> params.imapfile; 
			break; 
		case 'o': 
			arg >> params.checkpoint_outs; 
			break; 
		case 'c':
			arg >> params.min_reads;
			break;
		case 'k':
			arg >> params.k_value;
			break;
		case 'g':
			arg >> params.k_shift;
			break;
		case 'j':
			arg >> params.j_index;
			break;
		case 'l':
			arg >> params.min_links;
			break;
		case 'z':
			arg >> params.min_size;
			break;
		case 'b':
			arg >> params.base_name;
			break;
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
			arg >> params.max_degree;
			break;
		case 'e':
			arg >> params.end_length;
			break;
		case 'r':
			arg >> params.error_percent;
			break;
		case 'v':
			++params.verbose;
			break;
		case 't':
			arg >> params.threads;
			break;
		case OPT_HELP:
			std::cout << USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		case OPT_VERSION:
			std::cout << VERSION_MESSAGE;
			exit(EXIT_SUCCESS);
		}
		if (optarg != NULL && (!arg.eof() || arg.fail())) {
			std::cerr << PROGRAM ": invalid option: `-" << (char) c << optarg
					<< "'\n";
			exit(EXIT_FAILURE);
		}
	}
	omp_set_num_threads(params.threads);

	//gather all the chromium files
	vector<string> inputFiles = convertInputString(rawInputFiles); 

	while (optind < argc) {
		inputFiles.push_back(argv[optind]); 
		optind++; 
	}

	std::ifstream g(params.file.c_str());
	if (!g.good()) {
		std::cerr << "Cannot find -f " << params.file << ". Exiting... \n";
		die = true;
	}

	if (params.program == "full") {
		full = true; 
	} else if (params.program == "align") {
		alignc = true; 
	} else if (params.program == "graph") {
		graph = true; 
	} else {
		std::cerr << "You must specify where you want ARKS to start. Exiting... \n"; 
		die = true;  
	}

	if (die) {
		std::cerr << "Try " << PROGRAM << " --help for more information.\n";
		exit(EXIT_FAILURE);
	}

	/* Setting base name if not previously set; name depends on type of arcs used */
	if (params.base_name.empty()) {
		std::ostringstream filename;
		filename << params.file << ".scaff" << "k-method" << "_c"
				<< params.min_reads << "_k" << params.k_value << "_g"
				<< params.k_shift << "_j" << params.j_index << "_l"
				<< params.min_links << "_d" << params.max_degree << "_e"
				<< params.end_length << "_r" << params.error_percent;
		params.base_name = filename.str();
	}

	printf("%s\n", "Finished reading user inputs...entering runArcs()...");

	runArcs(inputFiles);

	return 0;
}
