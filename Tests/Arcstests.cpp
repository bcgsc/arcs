/* 
 *  Arcs.cpp test for kmer arcs
 *
 */

#include "Arcs.h"
#include <assert.h>
#include <string>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include "Common/kseq.h"
#include <cassert>

KSEQ_INIT(gzFile,gzread)

using namespace std; 

int numkmersduplicate, numkmersmapped; 

int main(int argc, char ** argv) {
	
/* Tests for the Helpers */

	// Testing HeadOrTail: 
	string r0 = HeadOrTail(true); 
	assert(r0==" Head"); 
	result = HeadOrTail(false); 
	assert(r0==" Tail"); 

	// Testing checkIndex:
	bool r1 = checkIndex("ATGC"); 
	assert(r1); 
	r1 = checkIndex("atgc"); 
	assert(r1); 
	r1 = checkIndex("AaTtGgCc"); 
	assert(r1); 
	r1 = checkIndex("ABCD"); 
	assert(!r1);
	r1 = checkIndex("1234"); 
	assert(!r1); 

	//Testing checkContigSequence: 
	r1 = checkContigSequence("ATGCNMRWSYKVHDB"); 
	assert(r1); 
	r1 = checkContigSequence("AaTtGgCcNnMmRrWwSsYyKkVvHhDdBb"); 
	assert(r1); 
	r1 = checkContigSequence("ABCDslaksjeo13409"); 
	assert(!r1); 
	r1 = checkContigSequence("123(*"); 
	assert(!r1); 
	r1 = checkContigSequence("@"); 
	assert(!r1); 

	//Testing checkReadSequence: 
	// TODO: implement allowable ambiguity (future)
	r1 = checkReadSequence("ATGCN"); 
	assert(r1); 
	r1 = checkReadSequence("AaTtGgCcNn"); 
	assert(r1); 
	r1 = checkReadSequence("ABCdsla"); 
	assert(!r1); 
	r1 = checkReadSequence("12349"); 
	assert(!r1); 


/* Tests for nonfile dependent ARCS Process Functions */

	// Setup: 
	ARCS::ContigKMap& kmap; 
	kmap.set_deleted_key(""); 
	ARCS::IndexMapimap; 
	ARCS::PairMap pmap; 
	unordered_map<string, int> indexMultMap; 

	
	//Testing mapKmers: 
	int16_t k_proc = 5;
	ReadsProcessor proc(k_proc); 
	int r2 = mapKmers("AAAAAAAAAA", 5, 1, kmap, proc, 2); 
	const unsigned char* r3 = proc.prepSeq("AAAAAAAAAA", 0); 
	r0 = proc.getStr(r3); 
	assert(kmap[r0] == 0); 
	assert(r2 == 6); 
	assert(numkmersduplicate == 5); 
	assert(numkmersmapped == 1); 

	r2 = mapKmers("ATGCGCT", 5, 1, kmap, proc, 2); 
	r0 = proc.getStr(proc.prepSeq("ATGCGCT",0));		//ATGCG CGCAT
	assert(kmap[r0] == 2); 			
	r0 = proc.getStr(proc.prepSeq("ATGCGCT",1));		//TGCGC GCGCA
	assert(kmap[r0] == 2);  
	r0 = proc.getStr(proc.prepSeq("ATGCGCT",2));		//GCGCT AGCGC
	assert(kmap[r0] == 2); 
	assert(r2 == 3); 
	assert(numkmersduplicate == 5); 
	assert(numkmersmapped == 4); 
	
	r2 = mapKmers("NACGCT", 5, 1, kmap, proc, 2); 
	r0 = proc.getStr(proc.prepSeq("NACGCT",1));		//ACGCT AGCGT
	assert(kmap[r0] == 2); 
	assert(r2 == 1); 
	//TODO: have it print out the entire map and manually check

	

	//Testing calcJacIndex:
	double r4 = calcJacIndex(1, 2); 
	assert(r4 == 0.5); 
	r4 = calcJacIndex(55, 100); 
	assert(r4 == 0.55); 
	r4 = calcJacIndex(0, 100); 
	assert(r4 == 0); 
	r4 = calcJacIndex(0, 0); 
	assert(r4 == 0); 

	//Testing bestContig: 
	r2 = bestContig(kmap, "AAAAA", 5, 1, 0.55, kmap); 
	assert(r2 == 0); 
	r2 = bestContig(kmap, "ATGCG", 5, 1, 0.55, kmap); 
	assert(r2 == 2); 
	r2 = bestContig(kmap, "CGCAT", 5, 1, 0.55, kmap); 
	assert(r2 == 2); 
	r2 = bestContig(kmap, "TGCGC", 5, 1, 0.55, kmap); 
	assert(r2 == 2); 
	r2 = bestContig(kmap, "GCGCA", 5, 1, 0.55, kmap); 
	assert(r2 == 2); 
	r2 = bestContig(kmap, "GCGCT", 5, 1, 0.55, kmap); 
	assert(r2 == 2); 
	r2 = bestContig(kmap, "AGCGC", 5, 1, 0.55, kmap); 
	assert(r2 == 2); 
	r2 = bestContig(kmap, "ACGCT", 5, 1, 0.55, kmap); 
	assert(r2 == 2); 
	r2 = bestContig(kmap, "AGCGT", 5, 1, 0.55, kmap); 
	assert(r2 == 2); 

	return 0; 
}













	
