#include "Arcs_work.h"
// for the inputted Fasta or Fastq file of contigs (if it is compressed with gzip)
#include <zlib.h>
#include "kseq.h"
#include <cassert>

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
// -p tells us if we are using read alignment method (preprocessed input with bwa) or we are using kmer (no preprocessing)
//"   -p  If 0 then use read alignment; If 1 then use kmer (so default is 0)\n"
"   -f  Using kseq parser, these are the contig sequences to further scaffold and can be in either FASTA or FASTQ format.\n"
"		If using read alignment option, then must be Multi-Fasta format\n"
/*"   -a  File of File Names listing all input BAM alignment files (required if you are using read alignment option). \n"
"       NOTE: alignments must be sorted in order of name\n"
"             index must be included in read name in the format read1_indexA\n" */
// use -h if we are using the k-merizing method rather than read alignment
"   -h  Chromium read file (output from longranger) (required if using k-mer option)\n"
//"   -s  Minimum sequence identity (min. required to include the read's scaffold alignment in the graph file, default: 98)\n"
"   -c  Minimum number of mapping read pairs/Index required before creating edge in graph. (default: 5)\n"
// add k-value required by the user
"   -k  k-value for the size of a k-mer. (default: 30) (required)\n"
// add shift between k-mers (optional supplied by user)
"   -g  shift between k-mers (default: 1)\n"
"   -j  Minimum Jaccard Index for a read to be associated with a contigId. (default: 0.55)\n"
"   -l  Minimum number of links to create edge in graph (default: 0)\n"
"   -z  Minimum contig length to consider for scaffolding (default: 500)\n"
"   -b  Base name for your output files (optional)\n"
"   -m  Range (in the format min-max) of index multiplicity (only reads with indices in this multiplicity range will be included in graph) (default: 50-10000)\n"
"   -d  Maximum degree of nodes in graph. All nodes with degree greater than this number will be removed from the graph prior to printing final graph. For no node removal, set to 0 (default: 0)\n"
"   -e  End length (bp) of sequences to consider (default: 30000)\n"
"   -r  Maximum p-value for H/T assignment and link orientation determination. Lower is more stringent (default: 0.05)\n"
"   -v  Runs in verbose mode (optional, default: 0)\n";


ARCS::ArcsParams params;

static const char shortopts[] = "f:h:s:c:k:g:j:l:z:b:m:d:e:r:v";

enum { OPT_HELP = 1, OPT_VERSION};

static const struct option longopts[] = {
    //{"arcs_type", required_argument, NULL, 'p'},    
    {"file", required_argument, NULL, 'f'},
    //{"fofName", required_argument, NULL, 'a'},
    {"c_input", required_argument, NULL, 'h'}, 
    //{"seq_id", required_argument, NULL, 's'}, 
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
    {"run_verbose", required_argument, NULL, 'v'},
    {"version", no_argument, NULL, OPT_VERSION},
    {"help", no_argument, NULL, OPT_HELP},
    { NULL, 0, NULL, 0 }
};

/* Writes things to a log file */
void writeToLogFile(const std::string & text) {
	std::string logfilename = params.base_name + ".log_file.txt";
	std::ofstream log_file(
		logfilename, std::ios_base::out | std::ios_base::app); 
	log_file << text << std::endl; 
}

/*
void writeToKmerLogFile(const std::string & text) {
	std::ofstream log_file(
		"Arcs_log_file.txt", std::ios_base::out | std::ios_base::app); 
	log_file << text << std::endl; 
}
*/

void writeIndexMultToLog(std::unordered_map<std::string, int> indexMultMap) {
	writeToLogFile("=>Index Multiplicity Map: "); 
	for (auto it = indexMultMap.begin(); it != indexMultMap.end(); ++it) {
		writeToLogFile(it->first + "  " + std::to_string(it->second)); 
	}
}

std::string HeadOrTail(bool orientation) {
	if (orientation) {
		return " Head"; 
	} else {
		return " Tail"; 
	}
}

//typedef std::unordered_map<std::string, ScafMap> IndexMap; 
//typedef std::map<std::pair<std::string, bool>, int> ScafMap;
void writeImapToLog(ARCS::IndexMap IndexMap) {
	writeToLogFile("=>Index Map: ");
	for (auto it = IndexMap.begin(); it != IndexMap.end(); ++it) {
		std::string barcode = it->first; 
		writeToLogFile(barcode);
		ARCS::ScafMap smap = it->second; 
		for (auto j = smap.begin(); j != smap.end(); ++j) {
			std::string contigname = j->first.first;
			std::string orientation = HeadOrTail(j->first.second); 
			writeToLogFile("    " + contigname + orientation + "     " + std::to_string(j->second)); 
		}
	}
}

void writePairMapToLog(ARCS::PairMap pmap) {
	writeToLogFile("=>Pair Map: "); 
	for (auto it = pmap.begin(); it != pmap.end(); ++it) {
		writeToLogFile("\t" + it->first.first + " and " + it->first.second); 
		writeToLogFile("\t\t\tHead-Head " + std::to_string(it->second[0])); 
		writeToLogFile("\t\t\tHead-Tail " + std::to_string(it->second[1])); 
		writeToLogFile("\t\t\tTail-Head " + std::to_string(it->second[2])); 
		writeToLogFile("\t\t\tTail-Tail " + std::to_string(it->second[3]));
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
        if (c != 'A' && c != 'T' && c != 'G' && c != 'C' && c != 'N' &&
		c != 'M' && c != 'R' && c != 'W' && c != 'S' && c != 'Y' && c != 'K' && 
		c != 'V' && c != 'H' && c != 'D' && c != 'B') {
		std::cout << c << std::endl; 
	        return false;
	}
    }

    return true;
}

/* Returns true if seqence only contains ATGC or N
 * 	Allows only if 0.98 of the read is not ambiguous (aka not N) 
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

/* Returns the size of the array for storing contigs */
int initContigArray(std::string contigfile) {

	int count = 0; 

	gzFile fp;  

	int l; 
	const char* filename = contigfile.c_str(); 
	fp = gzopen(filename, "r"); 
	kseq_t * seq = kseq_init(fp); 

	while((l= kseq_read(seq)) >= 0) {
		std::string sequence = seq->seq.s;
		int sequence_length = sequence.length();
		if (checkContigSequence(sequence) && sequence_length >= params.min_size )
			count++;
	}
	kseq_destroy(seq); 
	gzclose(fp); 

	//Let the first index represent the collision
	return (count * 2) + 1; 
}

int numkmersmapped = 0, numkmersduplicate = 0; 

/* Shreds end sequence into kmers and inputs them one by one into the ContigKMap 
 * 	std::pair<std::string, bool> 				specifies contigID and head/tail 
 * 	std::string						the end sequence of the contig 
 *	int							k-value specified by user 
 *	ARCS::ContigKMap					ContigKMap for storage of kmers
 */
int mapKmers(std::string seqToKmerize, int k, int k_shift, ARCS::ContigKMap& kmap,
		 ReadsProcessor &proc, int conreci) {

	int seqsize = seqToKmerize.length();  

	// Checks if length of the subsequence is smaller than the k size
	// 	If the contig end is shorter than the k-value, then we ignore it for now. 
	if (seqsize < k) {
		std::string errmsg = "Warning: ends of contig is shorter than k-value for contigID (no k-mers added): ";
		std::cout << errmsg << conreci << std::endl; 

		return 0; 
	} else {
		//assert(seqsize >= k); 
		int numKmers = 0; 

		int i = 0; 
		while (i <= seqsize - k) {
			const unsigned char* temp = proc.prepSeq(seqToKmerize, i);  // prepSeq returns NULL if it contains an N
			// Ignore a NULL kmer
			if (temp != NULL) {
				std::string kmerseq = proc.getBases(temp);

				//writeToKmerLogFile(proc.getBases(temp)); 

				numKmers++; 

				if (kmap.count(kmerseq) == 1) {
					numkmersduplicate++; 
					if (kmap[kmerseq] == conreci) {
						//std::string warningmsg = "kmer duplicate and not in same contigID: " + kmerseq;
						//writeToLogFile(warningmsg); 
						kmap[kmerseq] = 0; 
	
					}
				} else {
					assert(kmap.count(kmerseq) == 0); 
					kmap[kmerseq] = conreci; 
					numkmersmapped++; 
				}
				i += k_shift; 	
			} else {
				i += k; 
			}
		}
		return numKmers; 
	}
}

/* Get the k-mers from the paired ends of the contigs and store them in map. 
 * 	std::string file					FASTA (or later FASTQ) file 
 *	std::sparse_hash_map<k-mer, pair<contidID, bool>> 	ContigKMap
 *	int k							k-value (specified by user)
 */ 
void getContigKmers(std::string contigfile, ARCS::ContigKMap& kmap, ReadsProcessor &proc,
			std::vector<std::pair<std::string, bool>> &contigRecord){
	
	int totalNumContigs = 0; 
	int skippedContigs = 0; 
	int validContigs = 0; 
	int totalKmers = 0; 

	std::pair<std::string, bool> collisionmarker("null contig", false); 
	int conreci = 0; 
	contigRecord[conreci] = collisionmarker;

	gzFile fp;  

	int l; 
	const char* filename = contigfile.c_str(); 
	fp = gzopen(filename, "r"); 
	kseq_t * seq = kseq_init(fp); 

	while((l= kseq_read(seq)) >= 0) {
		std::string contigID = seq->name.s; 
		std::string sequence = seq->seq.s; 

		totalNumContigs++; 

		if (!checkContigSequence(sequence)) {
			std::string errormsg = "Error: Contig contains non-base characters. Please check your draft genome input file.";
			if (params.verbose) {
				std::cerr << contigID << ": " << errormsg << std::endl; 
			}
			writeToLogFile(contigID + ": " + errormsg); 
			skippedContigs++; 
		} else {

			// If the sequence is above minimum contig length, then will extract kmers from both ends 
			// If not will ignore the contig
			int sequence_length = sequence.length();

			if (sequence_length >= params.min_size) {
			
				// If contig length is less than 2 x end_length, then we split the sequence 
				// in half to decide head/tail (aka we changed the end_length)
				int cutOff = params.end_length; 
				if (cutOff == 0 || sequence_length <= cutOff * 2)
					cutOff = sequence_length / 2; 
	
				// Arbitrarily assign head or tail to ends of the contig
				std::pair<std::string, bool> headside(contigID, true); 
				std::pair<std::string, bool> tailside(contigID, false); 

				//get ends of the sequence and put k-mers into the map
				conreci++; 
				contigRecord[conreci] = headside;
				std::string seqend; 
				seqend = sequence.substr(0, cutOff); 
				//writeToLogFile(contigID + " head comments:"); 
				int num = mapKmers(seqend, params.k_value, params.k_shift, kmap, proc, conreci); 
				totalKmers += num; 

				conreci++; 
				contigRecord[conreci] = tailside; 
				seqend = sequence.substr(sequence_length - cutOff, sequence_length);
				//writeToLogFile(contigID + " tail comments:"); 
				num = mapKmers(seqend, params.k_value, params.k_shift, kmap, proc, conreci); 
				totalKmers += num; 
	
				validContigs++; 
			} else {
				skippedContigs++; 
			}
		}

		// printprogress
		if (params.verbose) {
			if (totalNumContigs % 1000 == 0)
				printf("Finished %d Contigs...\n", totalNumContigs); 
		} 
	}	
	kseq_destroy(seq); 
	gzclose(fp); 

	if (params.verbose) {
		printf("%s %d\n%s %d\n%s %d\n%s %d\n%s %d\n%s %d\n", 
			"Total number of contigs in draft genome: ", totalNumContigs,
			"Total valid contigs: ", validContigs,  
			"Total skipped contigs: ", skippedContigs, 
			"Total number of Kmers: ", totalKmers, 
			"Number Kmers Recorded: ", numkmersmapped, 
			"Number Duplicate Kmers (recorded by appearance): ", numkmersduplicate); 
	}

	writeToLogFile("Total number of contigs in draft genome: " + std::to_string(totalNumContigs)
			 + "\nTotal valid contigs: " + std::to_string(validContigs)
			 + "\nTotal skipped contigs: " + std::to_string(skippedContigs)
			 + "\nTotal number of Kmers: " + std::to_string(totalKmers)
			 + "\nNumber Kmers Recorded: " + std::to_string(numkmersmapped)
			 + "\nNumber Duplicate Kmers (recorded by appearance): " + std::to_string(numkmersduplicate)); 

}

void filterChromiumFile(std::string chromfile, std::unordered_map<std::string, int>& indexMultMap) {
	gzFile fp3;  
	int l; 
	std::string chromiumfile = chromfile;
	const char* filename = chromiumfile.c_str(); 
	fp3 = gzopen(filename, "r"); 
	kseq_t * seq3 = kseq_init(fp3); 
	
	int numreads = 0; 
	int numreadskept = 0; 
	while((l= kseq_read(seq3)) >= 0) {
		numreads++;  
		std::string comment = ""; 
		std::string barcode = ""; 
		if (seq3->comment.l) {
			comment = seq3->comment.s; 
			// Extract the barcode which is in the comment under this format: 
			// 		BX:Z:<barcodeseq>-1
			barcode.clear(); 
			for (std::string::iterator i = comment.begin(); i != comment.end(); i++) {
				if (*i != 'B' && *i != 'X' && *i != ':' && *i != 'Z'
					&& *i != '-' && *i != '1' && *i != '\n') {
					barcode += *i;
				}
			}	
		}
		
		if (!barcode.empty() && checkIndex(barcode)) {
			numreadskept++; 
			indexMultMap[barcode]++; 
		}
	}

	kseq_destroy(seq3); 
	gzclose(fp3); 

	std::string msg = "Total Number of Reads in Chromium File: " + std::to_string(numreads) + "\nNumber of Reads that have Valid Barcodes: " + std::to_string(numreadskept); 
	writeToLogFile(msg); 

	if (params.verbose) {
		printf(msg.c_str()); 
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
int bestContig(ARCS::ContigKMap &kmap, std::string readseq,
			 int k, int k_shift, double j_index, ReadsProcessor &proc) {

	// to keep track of what contig+H/T that the k-mer from barcode matches to
	// 	int					Index that corresponds to the contig in the contigRecord
	// 	count					# kmers found here
	std::map<int, int> ktrack; 

	// k-merize readsequence
	int corrConReci = 0;
	int totalnumkmers = 0; 
	int numkmersfound = 0;
	int seqlen = readseq.length(); 

	int i = 0; 
	while (i <= seqlen-k) {
		const unsigned char* temp = proc.prepSeq(readseq, i); 
		if (temp != NULL) {

			std::string ckmerseq = proc.getBases(temp); 
 
			// search for kmer in ContigKmerMap and only record if it is not the collisionmaker */
			if (kmap.count(ckmerseq) == 1) {
				corrConReci = kmap[ckmerseq]; 
				if (corrConReci != 0)
					ktrack[corrConReci]++; 
				numkmersfound++; 
			}
		}
		totalnumkmers++;
		i += k_shift; 
	}

	std::string msg = "Total Number of Kmers from read: " + std::to_string(totalnumkmers)
			 + "\nNumber of kmers found in draft genome: " + std::to_string(numkmersfound); 
	writeToLogFile(msg); 

	/*if (params.verbose) {
		printf("%s\n", msg.c_str()); 
	}*/
	
	double maxjaccardindex = 0; 
	// for the read, find the contig that it is most compatible with based on the jaccard index
	for (auto it = ktrack.begin(); it != ktrack.end(); ++it) {
		double currjaccardindex = calcJacIndex(it->second, totalnumkmers); 
		if (maxjaccardindex < currjaccardindex){
			maxjaccardindex = currjaccardindex; 
			corrConReci = it->first; 
		}
	}

	// default accuracythreshold is 0.55)
	if (maxjaccardindex >= j_index) {
		/*if (params.verbose) {
			std::string ht = HeadOrTail(corrContigId.second); 
			printf("%s %s %s\n%s %f\n", "Corresponding contig: ", corrContigId.first.c_str(), ht.c_str(), 
				"Jaccard Index higher than threshold: ", maxjaccardindex);
		}*/
		//writeToLogFile("Best jaccard index = " + std::to_string(maxjaccardindex)); 
		return corrConReci; 
	} else {
		/*if (params.verbose) {
			printf("%s %f\n", "Jaccard Index lower than threshold: ", maxjaccardindex); 
		}
		writeToLogFile("Low jaccard index = " + std::to_string(maxjaccardindex)); */
		return 0; 
	}
}

/* Read through longranger basic chromium output fastq file */ 
void chromiumRead(std::string chromiumfile, ARCS::ContigKMap& kmap, ARCS::IndexMap& imap, 
			     std::unordered_map<std::string, int> indexMultMap,
			     ReadsProcessor &proc, 
			     std::vector<std::pair<std::string, bool>> &contigRecord) {

	int ctpername = 1;
	int skipped_unpaired = 0; 
	int skipped_invalidbarcode = 0; 
	int skipped_nogoodcontig = 0; 
	std::string prevname; 
	int prevConReci = 0; 

	int count = 0; 

	gzFile fp2;  
	int l; 
	const char* filename = chromiumfile.c_str(); 
	fp2 = gzopen(filename, "r"); 
	kseq_t * seq2 = kseq_init(fp2); 
	while((l= kseq_read(seq2)) >= 0) {
		count++; 
		if (params.verbose) {
			if (count % 1000000 == 0) 
				std::cout << "Processed " << count << " reads." << std::endl; 
		}

		std::string name = seq2->name.s; 
		/*if (params.verbose) {
			std::cout << name << std::endl; 
		}*/
		std::string comment = ""; 
		std::string barcode = ""; 
		if (seq2->comment.l) {
			comment = seq2->comment.s; 
			// Extract the barcode which is in the comment under this format: 
			// 		BX:Z:<barcodeseq>-1
			barcode.clear(); 
			for (std::string::iterator i = comment.begin(); i != comment.end(); i++) {
				if (*i != 'B' && *i != 'X' && *i != ':' && *i != 'Z'
					&& *i != '-' && *i != '1' && *i != '\n'){
					barcode += *i;
				}
			}	
		}
	
		/*if (params.verbose) {
			std::cout << "\t\t" << barcode << std::endl; 
		}*/

		writeToLogFile(name + "\t" + barcode); 
         	int indexMult = indexMultMap[barcode]; 
       		if (indexMult < params.min_mult || indexMult > params.max_mult) {
			writeToLogFile("Skipped Kmerization of barcode: " + barcode + " because not a good barcode."); 
			skipped_invalidbarcode++; 
			ctpername = 1; 
		} else {		
		
			std::string cread = seq2->seq.s; 

			if (!checkReadSequence(cread)) {
				/*std::string errmsg = "Warning: Skipped poor quality read: "
							+ name ;
				std::cerr << errmsg << std::endl; 
				writeToLogFile(errmsg); */
			} else {
				// Checks that it is a read pair
				if (ctpername == 2 && name != prevname) {
					skipped_unpaired ++; 
					/*std::string warningmsg = "Warning: Skipping an unpaired read. File should be sorted in order of read name.\n\tPrev read: " + prevname + "\n\tCurr read: " + name; 
					writeToLogFile(warningmsg); */				
				
					// reset count
					ctpername = 1; 
				}
	
				int corrConReci = 0;
				std::pair<std::string, bool> corrContigId;
				if (ctpername >=3) 
					ctpername = 1; 
				if (ctpername == 1) {
					if (name.compare(prevname) != 0) {
						prevname = name; 
	
						// get the contig that the read matches
						corrConReci = bestContig(kmap, cread, params.k_value, params.k_shift,
										params.j_index, proc); 

						/*std::string ht = HeadOrTail(corrContigId.second); 
						//std::cout << ht << std::endl; 
						writeToLogFile("-->Best contig for read1: " + corrContigId.first + ht); */
						
						// store the corrContigId so that we can match it with its pair	
						// 	if corrContigId was not chosen, it will still set prevContig to NULL
						prevConReci = corrConReci; 
					} else {
						prevname = ""; 
						ctpername = 0; 
					}
				} else if (ctpername == 2) {
					assert(name == prevname);
				
					// get the contig that the read matches
					corrConReci = bestContig(kmap, cread, params.k_value, params.k_shift, 
									params.j_index, proc); 

					// we only store barcode info in index map if read pairs have same contig + orientation	
					// and if the corrContigId is not NULL (because it is above accuracy threshold)
					if (corrConReci != 0 && corrConReci == prevConReci) {
						corrContigId = contigRecord[corrConReci]; 
						imap[barcode][corrContigId]++; 
						std::string ht = HeadOrTail(corrContigId.second);
						writeToLogFile("barcode: " + barcode + "\tContigID: "
								+ corrContigId.first + ht + "	"
								+ std::to_string(imap[barcode][corrContigId])); 
					} else {
						writeToLogFile("No contig match for barcode."); 
						skipped_nogoodcontig++; 
					}
		
					// reset previous records
					prevname = ""; 
					prevConReci = 0; 
				}
				ctpername++;
			} 
		}
		writeToLogFile("\n"); 

	}
	kseq_destroy(seq2); 
	gzclose(fp2); 


	if (skipped_unpaired > 0) {
		std::cerr << "Warning:: Skipped: " << skipped_unpaired << " unpaired reads. Check your input file." << std::endl; 
	}

	if (skipped_invalidbarcode > 0) {
		std::cerr << "Warning:: Skipped: " << skipped_invalidbarcode << " invalid barcode reads. Check your input file." << std::endl; 
	}

	if (params.verbose) {
		if (skipped_nogoodcontig > 0)
			printf("%s %d\n", "Skipped reads without good match: ", skipped_nogoodcontig);
	}

	writeToLogFile("Skipped: " + std::to_string(skipped_unpaired + skipped_invalidbarcode) + " unpaired reads." + 
			"\nSkipped: " + std::to_string(skipped_invalidbarcode) + " invalid barcode reads." + 
			"\nSkipped: " + std::to_string(skipped_nogoodcontig) + " reads without good match"); 

	//fordebugging
	//writeIndexMultToLog(indexMultMap);
	//writeImapToLog(imap); 

}

/* Track Memory Usage */

			
		
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
    return 0.5 * (1 + std::erf((x - mean)/(sd * std::sqrt(2))));
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
                    if (scafA != scafB) {
                        bool validA, validB, scafAhead, scafBhead;

                        std::tie(validA, scafAhead) = headOrTail(it->second[std::pair<std::string, bool>(scafA, true)], it->second[std::pair<std::string, bool>(scafA, false)]);
                        std::tie(validB, scafBhead) = headOrTail(it->second[std::pair<std::string, bool>(scafB, true)], it->second[std::pair<std::string, bool>(scafB, false)]);


                        if (validA && validB) {
                            std::pair<std::string, std::string> pair (scafA, scafB);
			    std::pair<std::string, std::string> reversepair (scafB, scafA); 

		    	    if (pmap.count(reversepair) != 0) {
				pair = reversepair; 
				bool temp = scafAhead; 
				scafAhead = scafBhead;
				scafBhead = temp; 
			    }
		    
	
                            if (pmap.count(pair) == 0) {
                                std::vector<int> init(4,0); 
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
    writePairMapToLog(pmap); 
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
    for(it = pmap.begin(); it != pmap.end(); ++it) {
        std::string scaf1, scaf2;
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
    std::cout << "Entered runArcs()..." << std::endl; 

    std::cout << "Running: " << PROGRAM << " " << VERSION 
        << "\n pid " << ::getpid()
	//<< "\n -p " << params.arcs_type
        << "\n -f " << params.file
        //<< "\n -a " << params.fofName
	<< "\n -h " << params.c_input
        //<< "\n -s " << params.seq_id 
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
        << "\n -v " << params.verbose << "\n";

    std::string parametertext = "contig file: " + params.file + "\nchromium read file: " 
			+ params.c_input + "\nminimum reads per contig: " 
			+ std::to_string(params.min_reads) + "\nk-value: " 
			+ std::to_string(params.k_value) + "\nk_shift: " + std::to_string(params.k_shift)
			+ "\nminimum links per pair: " + std::to_string(params.min_links) 
			+ "\nminimum contig size: " + std::to_string(params.min_size) 
			+ "\nfile base name: " + params.base_name + "\nminimum index multiplicity: " 
			+ std::to_string(params.min_mult) + "\nmaximum index multiplicity: " 
			+ std::to_string(params.max_mult) + "\nmax degree: " 
			+ std::to_string(params.max_degree) + "\ncontig end length: " 
			+ std::to_string(params.end_length) + "\nerror percentage: " + std::to_string(params.error_percent); 
    writeToLogFile(parametertext); 

    std::string graphFile = params.base_name + "_original.gv";

    ARCS::ContigKMap kmap; 
    kmap.set_deleted_key(""); 

    ARCS::IndexMap imap;
    ARCS::PairMap pmap;
    ARCS::Graph g;
    std::unordered_map<std::string, int> indexMultMap;

    std::time_t rawtime;


    std::cout << "\n---We are using KMER method.---\n" << std::endl; 	

    int16_t k_proc = params.k_value; 
    ReadsProcessor proc(k_proc);


    // Initialize the contigRecord
    time(&rawtime); 
    std::cout << "\n=>Allocating the Contig Record... " << ctime(&rawtime); 
    writeToLogFile("\n=>Allocating the Contig Record... "); 
    int size = initContigArray(params.file); 
    std::vector<std::pair<std::string, bool>> contigRecord (size); 
	
    // Read contig file, shred sequences into k-mers, and then map them 
    time(&rawtime); 
    std::cout << "\n=>Storing Kmers from Contig ends... " << ctime(&rawtime); 
    writeToLogFile("\n=>Storing Kmers from Contig ends... ");
    getContigKmers(params.file, kmap, proc, contigRecord); 
	 

     /* Attempt to filter through chromium reads and remove not good barcodes first */
    time(&rawtime); 
    std::cout << "\n=>Filtering Chromium FastQ file... " << ctime(&rawtime); 
    writeToLogFile("\n=>Filtering Chromium FastQ file... "); 
    filterChromiumFile(params.c_input, indexMultMap); 

    time(&rawtime); 
    std::cout << "\n=>Reading Chromium FASTQ file... " << ctime(&rawtime); 
    writeToLogFile("\n=>Reading Chromium FASTQ file... ");
    chromiumRead(params.c_input, kmap, imap, indexMultMap, proc, contigRecord); 


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

    printf("Reading user inputs... (in main)\n"); 

    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case '?':
                die = true; break;
	   /* case 'p': 
		arg >> params.arcs_type; break; */
            case 'f':
                arg >> params.file; break;
           /* case 'a':
                arg >> params.fofName; break; */
	    case 'h': 
		arg >> params.c_input; break; 
            /*case 's':
                arg >> params.seq_id; break; */
            case 'c':
                arg >> params.min_reads; break;
	    case 'k':
		arg >> params.k_value; break; 
	    case 'g':
		arg >> params.k_shift; break; 
	    case 'j': 
		arg >> params.j_index; break; 
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

    // if we have selected to use read alignment based arcs, then check to find correct files
    /*if (params.arcs_type == 0){
    	std::ifstream f(params.fofName.c_str());
    	if (!f.good()) {
    	    std::cerr << "Cannot find -a " << params.fofName << ". Exiting... \n";
    	    die = true;
    	}
    }*/
    // if we have selected to use k-mer based arcs, then check to find correct file
   // if (params.arcs_type == 1){
	std::ifstream f(params.c_input.c_str()); 
	if (!f.good()) {
		std::cerr << "Cannot find -h " << params.c_input << ". Exiting... \n"; 
		die = true; 
	}
   // }

    std::ifstream g(params.file.c_str());
    if (!g.good()) {
        std::cerr << "Cannot find -f " << params.file << ". Exiting... \n";
        die = true;
    }

    if (die) {
        std::cerr << "Try " << PROGRAM << " --help for more information.\n";
        exit(EXIT_FAILURE);
    }

    /* Setting base name if not previously set; name depends on type of arcs used */
    if (params.base_name.empty()) {
        std::ostringstream filename;
	//if (params.arcs_type) {
        	filename << params.file << ".scaff" 
          	//<< "_s" << params.seq_id 
		<< "k-method"
            	<< "_c" << params.min_reads
	   	<< "_k" << params.k_value     
		<< "_g" << params.k_shift
		<< "_j" << params.j_index  
           	<< "_l" << params.min_links 
           	<< "_d" << params.max_degree 
            	<< "_e" << params.end_length
           	<< "_r" << params.error_percent;
	/*} else {
        	filename << params.file << ".scaff" 
            	<< "_s" << params.seq_id 
		<< "alignment method"
            	<< "_c" << params.min_reads
            	<< "_l" << params.min_links 
            	<< "_d" << params.max_degree 
            	<< "_e" << params.end_length
            	<< "_r" << params.error_percent;
	}*/
        params.base_name = filename.str();
    }

    printf("%s\n", "Finished reading user inputs...entering runArcs()..."); 

    runArcs();

    return 0;
}
