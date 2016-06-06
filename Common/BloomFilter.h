/*
 *
 * BloomFilter.h
 *
 *  Created on: Aug 10, 2012
 *      Author: cjustin
 */

#ifndef BLOOMFILTER_H_
#define BLOOMFILTER_H_
#include <string>
#include <vector>
#include <stdint.h>
#include "city.h"
#include <math.h>
#if _OPENMP
# include <omp.h>
#endif

using namespace std;

static const uint8_t bitsPerChar = 0x08;
static const unsigned char bitMask[0x08] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20,
		0x40, 0x80 };

//TODO Work out better way to deal with kmerSize since conversion to kmerSize in bytes is needed

/*
 * For precomputing hash values. kmerSize is the number of bytes of the string used.
 */
static inline vector<size_t> multiHash(const unsigned char* kmer, size_t num,
		unsigned kmerSize) {
	vector<size_t> tempHashValues(num);
	//use raw kmer number as first hash value
	size_t kmerSizeInBytes = (kmerSize + 4 - 1) / 4;

	for (size_t i = 0; i < num; ++i) {
		tempHashValues[i] = CityHash64WithSeed(
				reinterpret_cast<const char*>(kmer), kmerSizeInBytes, i);
	}
	return tempHashValues;
}

class BloomFilter {
public:
	//for generating a new filter
	explicit BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize);
	void insert(vector<size_t> const &precomputed);
	void insert(const unsigned char* kmer);
	bool contains(vector<size_t> const &precomputed) const;
	bool contains(const unsigned char* kmer) const;

	unsigned getHashNum() const;
	unsigned getKmerSize() const;

	//for storing/restoring the filter
	void storeFilter(string const &filterFilePath) const;
	explicit BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize,
			string const &filterFilePath);

	virtual ~BloomFilter();
private:
	BloomFilter(const BloomFilter& that); //to prevent copy construction
	void initSize(size_t size);
	uint8_t* m_filter;
	size_t m_size;
	size_t m_sizeInBytes;
	unsigned m_hashNum;
	unsigned m_kmerSize;
	unsigned m_kmerSizeInBytes;
};

#endif /* BLOOMFILTER_H_ */
