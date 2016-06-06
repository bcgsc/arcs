/*
 * BloomFilterInfo.h
 * Intended to be used to collect/store and compute information about bloom filters
 * Can output information in IEE format text file
 *
 *  Created on: Aug 20, 2012
 *      Author: cjustin
 */

#ifndef BLOOMFILTERINFO_H_
#define BLOOMFILTERINFO_H_
#include <string>
#include <vector>
#include <boost/unordered/unordered_map.hpp>

using namespace std;

// For Calculating aspects of bloom filter
// todo: Tweak calculations as they are approximations and may not be 100% optimal
// see http://en.wikipedia.org/wiki/Bloom_filter

/*
 * Calculation assumes optimal ratio of bytes per entry given a fpr
 */
//Note: Rounded down because in practice you want to calculate as few hash values as possible
//static uint8_t calcOptiHashNum(float fpr) {
//	return uint8_t(-log(fpr) / log(2));
//}
class BloomFilterInfo {
public:
	explicit BloomFilterInfo(string const &filterID, unsigned kmerSize,
			unsigned hashNum, double desiredFPR, size_t expectedSize,
			const vector<string> &seqSrc);
	explicit BloomFilterInfo(string const &fileName);
	void addHashFunction(const string &fnName, size_t seed);
	void setRedundancy(size_t redunSeq);
	void setTotalNum(size_t totalNum);

	void printInfoFile(const string &fileName) const;
	virtual ~BloomFilterInfo();

	//getters
	unsigned getKmerSize() const;
	unsigned getHashNum() const;
	size_t getCalcuatedFilterSize() const;
	const string &getFilterID() const;
	const string &getPresetType() const;
	double getRedundancyFPR() const;
	double getFPR() const;

private:
	//user specified input
	string m_filterID;
	unsigned m_kmerSize;
	double m_desiredFPR;
	vector<string> m_seqSrcs;
	unsigned m_hashNum;
	size_t m_expectedNumEntries;

	//determined at run time
	struct runtime {
		size_t size;
		size_t numEntries;
		double FPR;
		size_t redundantSequences;
		double redundantFPR;
	};

	runtime m_runInfo;

	const vector<string> convertSeqSrcString(const string &seqSrcStr) const;
	double calcApproxFPR(size_t size, size_t numEntr,
			unsigned hashFunctNum) const;
	double calcRedunancyFPR(size_t size, size_t numEntr,
			unsigned hashFunctNum) const;
	size_t calcOptimalSize(size_t entries, float fpr) const;
	size_t calcOptimalSize(size_t entries, float fpr,
			unsigned hashNum) const;
};

#endif /* BLOOMFILTERINFO_H_ */
