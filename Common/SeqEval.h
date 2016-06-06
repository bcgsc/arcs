/*
 * SeqEval.h
 *
 * Algorithms for bloomfilter based evaluation on sequences
 *
 *  Created on: Mar 10, 2015
 *      Author: cjustin
 *
 * Todo: try to expand and transfer methods from BioBloomClassifier
 */

#ifndef SEQEVAL_H_
#define SEQEVAL_H_

#include <string>
#include <cmath>
#include <cassert>
#include "boost/unordered/unordered_map.hpp"
#include "DataLayer/FastaReader.h"
#include "Common/Options.h"
#include "Common/ReadsProcessor.h"
#include "Common/BloomFilter.h"

using namespace std;
using namespace boost;

namespace SeqEval {

inline double denormalizeScore(double score, unsigned kmerSize, size_t seqLen)
{
	assert(score >= 0 && score <= 1);
	return score * (seqLen - kmerSize + 1);
}

inline double normalizeScore(double score, unsigned kmerSize, size_t seqLen)
{
	return score / (seqLen - kmerSize + 1);
}

/*
 * Evaluation algorithm with hashValue storage (minimize redundant work)
 */
inline bool evalSingle(const FastqRecord &rec, unsigned kmerSize, const BloomFilter &filter,
		double threshold, double antiThreshold, unsigned hashNum,
		vector<vector<size_t> > *hashValues, const BloomFilter *subtract)
{
	threshold = denormalizeScore(threshold, kmerSize, rec.seq.length());
	antiThreshold = floor(denormalizeScore(antiThreshold, kmerSize, rec.seq.length()));

	ReadsProcessor proc(kmerSize);

	size_t currentLoc = 0;
	double score = 0;
	unsigned antiScore = 0;
	unsigned streak = 0;
	while (rec.seq.length() >= currentLoc + kmerSize) {
		const unsigned char* currentSeq = proc.prepSeq(rec.seq, currentLoc);
		if (streak == 0) {
			if (currentSeq != NULL) {
				vector<size_t> hash = multiHash(currentSeq, hashNum, kmerSize);
				if (hashValues != NULL)
					(*hashValues)[currentLoc] = hash;
				if ((subtract == NULL || !subtract->contains(hash))
						&& filter.contains(hash)) {
					score += 0.5;
					++streak;
					if (threshold <= score) {
						return true;
					}
				}
				else if (antiThreshold <= ++antiScore) {
					return false;
				}
				++currentLoc;
			} else {
				if (currentLoc > kmerSize) {
					currentLoc += kmerSize + 1;
					antiScore += kmerSize + 1;
				} else {
					++antiScore;
					++currentLoc;
				}
				if (antiThreshold <= antiScore) {
					return false;
				}
			}
		} else {
			if (currentSeq != NULL) {
				vector<size_t> hash = multiHash(currentSeq, hashNum, kmerSize);
				if (hashValues != NULL)
					(*hashValues)[currentLoc] = hash;
				if ((subtract == NULL || !subtract->contains(hash))
						&& filter.contains(hash)) {
					++streak;
					score += 1 - 1 / (2 * streak);
					++currentLoc;

					if (threshold <= score) {
						return true;
					}
					continue;
				}
				else if (antiThreshold <= ++antiScore) {
					return false;
				}
			} else {
				currentLoc += kmerSize + 1;
				antiScore += kmerSize + 1;
			}
			if (streak < opt::streakThreshold) {
				++currentLoc;
			} else {
				currentLoc += kmerSize;
				antiScore += kmerSize;
			}
			if (antiThreshold <= antiScore) {
				return false;
			}
			streak = 0;
		}
	}
	return false;
}

/*
 * Evaluation algorithm with hashValue storage (minimize redundant work)
 */
inline bool evalSingle(const FastqRecord &rec, unsigned kmerSize, const BloomFilter &filter,
		double threshold, double antiThreshold, unsigned hashNum,
		vector<vector<size_t> > &hashValues)
{
	return evalSingle(rec, kmerSize, filter, threshold, antiThreshold, hashNum,
			&hashValues, NULL);
}

/*
 * Evaluation algorithm with no hashValue storage (optimize speed for single queries)
 */
inline bool evalSingle(const FastqRecord &rec, unsigned kmerSize, const BloomFilter &filter,
		double threshold, size_t antiThreshold)
{
	return evalSingle(rec, kmerSize, filter, threshold, antiThreshold, filter.getHashNum(),
			NULL, NULL);
}

/*
 * Evaluation algorithm with hashValue storage (minimize redundant work)
 */
inline bool evalSingle(const FastqRecord &rec, unsigned kmerSize, const BloomFilter &filter,
		double threshold, double antiThreshold, unsigned hashNum,
		vector<vector<size_t> > &hashValues, const BloomFilter &subtract)
{
	return evalSingle(rec, kmerSize, filter, threshold, antiThreshold, hashNum,
			&hashValues, &subtract);
}

/*
 * Evaluation algorithm based on minimum number of contiguous matching bases.
 */
inline bool evalMinMatchLen(const FastqRecord &rec, unsigned kmerSize, const BloomFilter &filter,
		unsigned minMatchLen, unsigned hashNum, vector<vector<size_t> > *hashValues,
		const BloomFilter *subtract)
{
	ReadsProcessor proc(kmerSize);
	// number of contiguous k-mers matched
	unsigned matchLen = 0;
	size_t l = rec.seq.length();

	for (size_t i = 0; i < l + kmerSize - 1; ++i) {
		// quit early if there is no hope
		if (l - i + matchLen < minMatchLen)
			return false;
		// get 2-bit encoding of k-mer at index i
		const unsigned char* kmer = proc.prepSeq(rec.seq, i);
		if (kmer == NULL) {
			matchLen = 0;
			continue;
		}
		// compute Bloom filter hash functions
		vector<size_t> hash = multiHash(kmer, hashNum, kmerSize);
		// cache hash values for future use
		if (hashValues != NULL)
			(*hashValues)[i] = hash;
		// ignore k-mers in subtract filter
		if (subtract != NULL && subtract->contains(hash)) {
			matchLen = 0;
			continue;
		}
		if (filter.contains(hash)) {
			if (matchLen == 0)
				matchLen = kmerSize;
			else
				++matchLen;
		}
		else {
			matchLen = 0;
		}
		// if min match length reached
		if (matchLen >= minMatchLen)
			return true;
	}
	return false;
}

enum EvalMode { EVAL_STANDARD, EVAL_MIN_MATCH_LEN };

inline bool evalRead(const FastqRecord &rec, unsigned kmerSize, const BloomFilter &filter,
		double threshold, double antiThreshold, unsigned hashNum,
		vector<vector<size_t> > *hashValues, const BloomFilter *subtract, EvalMode mode)
{
	// compute enough hash values to check both the main Bloom filter
	// and the subtractive Bloom filter
	if (subtract != NULL) {
		if (subtract->getHashNum() > hashNum) {
			hashNum = subtract->getHashNum();
		}
	}

	switch(mode) {
	case EVAL_MIN_MATCH_LEN:
		return evalMinMatchLen(rec, kmerSize, filter, (unsigned)round(threshold),
			hashNum, hashValues, subtract);
	case EVAL_STANDARD:
	default:
		return evalSingle(rec, kmerSize, filter, threshold, antiThreshold,
			hashNum, hashValues, subtract);
	}
	assert(false);
}

inline bool evalRead(const FastqRecord &rec, unsigned kmerSize, const BloomFilter &filter,
		double threshold, double antiThreshold, unsigned hashNum,
		vector<vector<size_t> > &hashValues, const BloomFilter &subtract, EvalMode mode)
{
	return evalRead(rec, kmerSize, filter, threshold, antiThreshold, hashNum,
		&hashValues, &subtract, mode);
}

inline bool evalRead(const FastqRecord &rec, unsigned kmerSize, const BloomFilter &filter,
		double threshold, double antiThreshold, unsigned hashNum,
		vector<vector<size_t> > &hashValues, EvalMode mode)
{
	return evalRead(rec, kmerSize, filter, threshold, antiThreshold, hashNum,
		&hashValues, NULL, mode);
}

inline bool evalRead(const FastqRecord &rec, unsigned kmerSize, const BloomFilter &filter,
		double threshold, double antiThreshold, EvalMode mode)
{
	return evalRead(rec, kmerSize, filter, threshold, antiThreshold,
		filter.getHashNum(), NULL, NULL, mode);
}

/*
 * Evaluation algorithm with no hashValue storage (optimize speed for single queries)
 * Returns score and does not have a stopping threshold
 */
inline double evalSingleExhaust(const FastqRecord &rec, unsigned kmerSize,
		const BloomFilter &filter)
{
	ReadsProcessor proc(kmerSize);
	size_t currentLoc = 0;
	double score = 0;
	unsigned streak = 0;
	while (rec.seq.length() >= currentLoc + kmerSize) {
		const unsigned char* currentKmer = proc.prepSeq(rec.seq, currentLoc);
		if (streak == 0) {
			if (currentKmer != NULL) {
				if (filter.contains(currentKmer)) {
					score += 0.5;
					++streak;
				}
				++currentLoc;
			} else {
				currentLoc += kmerSize + 1;
			}
		} else {
			if (currentKmer != NULL) {
				if (filter.contains(currentKmer)) {
					++streak;
					score += 1 - 1 / (2 * streak);
					++currentLoc;
					continue;
				}
			} else {
				currentLoc += kmerSize + 1;
			}
			if (streak < opt::streakThreshold) {
				++currentLoc;
			} else {
				currentLoc += kmerSize;
			}
			streak = 0;
		}
	}
	return score;
}

///*
// * Evaluation algorithm with no hashValue storage (optimize speed for single queries)
// * Returns length when stopping threshold is met
// * Returns the length of the read where threshold was met
// */
//inline unsigned evalSingleLength(const FastqRecord &rec, unsigned kmerSize,
//		const BloomFilter &filter, double threshold, double antiThreshold)
//{
//	ReadsProcessor proc(kmerSize);
//	size_t currentLoc = 0;
//	double score = 0;
//	unsigned antiScore = 0;
//	unsigned streak = 0;
//	while (rec.seq.length() >= currentLoc + kmerSize) {
//		const unsigned char* currentKmer = proc.prepSeq(rec.seq, currentLoc);
//		if (streak == 0) {
//			if (currentKmer != NULL) {
//				if (filter.contains(currentKmer)) {
//					score += 0.5;
//					++streak;
//					if (threshold <= score) {
//						return currentLoc;
//					}
//				} else if (antiThreshold <= ++antiScore) {
//					return rec.seq.length() - kmerSize;
//				}
//				++currentLoc;
//			} else {
//				if (currentLoc > kmerSize) {
//					currentLoc += kmerSize + 1;
//					antiScore += kmerSize + 1;
//				} else {
//					++antiScore;
//					++currentLoc;
//				}
//				if (antiThreshold <= antiScore) {
//					return rec.seq.length() - kmerSize;
//				}
//			}
//		} else {
//			if (currentKmer != NULL) {
//				if (filter.contains(currentKmer)) {
//					++streak;
//					score += 1 - 1 / (2 * streak);
//
//					if (threshold <= score) {
//						return currentLoc;
//					}
//					++currentLoc;
//					continue;
//				} else if (antiThreshold <= ++antiScore) {
//					return rec.seq.length() - kmerSize;
//				}
//			} else {
//				currentLoc += kmerSize + 1;
//				antiScore += kmerSize + 1;
//			}
//			if (streak < opt::streakThreshold) {
//				++currentLoc;
//			} else {
//				currentLoc += kmerSize;
//				antiScore += kmerSize;
//			}
//			if (antiThreshold <= antiScore) {
//				return rec.seq.length() - kmerSize;
//			}
//			streak = 0;
//		}
//	}
//	return rec.seq.length() - kmerSize;
//}

/*
 * Core evaluation algorithm, with ability start evaluating sequence midway
 * Evaluation algorithm with hashValue storage (minimize redundant work)
 * Also stores if position has already been visited to minimize work
 * Takes in last position visited and score and updates them accordingly
 */
inline bool eval(const FastqRecord &rec, unsigned kmerSize,
		const BloomFilter &filter, double threshold, double antiThreshold,
		vector<bool> &visited, vector<vector<size_t> > &hashValues,
		unsigned &currentLoc, double &score, ReadsProcessor &proc)
{
	threshold = denormalizeScore(threshold, kmerSize, rec.seq.length());
	antiThreshold = denormalizeScore(antiThreshold, kmerSize, rec.seq.length());
	score = denormalizeScore(score, kmerSize, rec.seq.length());

	unsigned antiScore = 0;
	unsigned streak = 0;
	bool hit = false;

	while (rec.seq.length() >= currentLoc + kmerSize) {

		//prepare hash values for filter

		//check if hash value is already generated
		if (hashValues[currentLoc].size() == 0) {
			if (!visited[currentLoc]) {
				const unsigned char* currentSeq = proc.prepSeq(rec.seq,
						currentLoc);
				if (currentSeq != NULL) {
					hashValues[currentLoc] = multiHash(currentSeq, filter.getHashNum(),
							kmerSize);
				}
				visited[currentLoc] = true;
			}
		}

		if (streak == 0) {
			if (hashValues[currentLoc].size() > 0) {
				if (filter.contains(hashValues[currentLoc])) {
					score += 0.5;
					++streak;
					if (threshold <= score) {
						++currentLoc;
						hit = true;
						break;
					}
				}
				else if (antiThreshold <= ++antiScore) {
					++currentLoc;
					hit = false;
					break;
				}
				++currentLoc;
			} else {
				if (currentLoc > kmerSize) {
					currentLoc += kmerSize + 1;
					antiScore += kmerSize + 1;
				} else {
					++antiScore;
					++currentLoc;
				}
				if (antiThreshold <= antiScore) {
					hit = false;
					break;
				}
			}
		} else {
			if (hashValues[currentLoc].size() > 0) {
				if (filter.contains(hashValues[currentLoc])) {
					++streak;
					score += 1 - 1 / (2 * streak);
					++currentLoc;

					if (threshold <= score) {
						hit = true;
						break;
					}
					continue;
				}
				else if (antiThreshold <= ++antiScore) {
					++currentLoc;
					hit = false;
					break;
				}
			} else {
				//if has non atcg character
				currentLoc += kmerSize + 1;
				antiScore += kmerSize + 1;
			}
			if (streak < opt::streakThreshold) {
				++currentLoc;
			} else {
				currentLoc += kmerSize;
				antiScore += kmerSize;
			}
			if (antiThreshold <= antiScore) {
				hit = false;
				break;
			}
			streak = 0;
		}
	}
	score = normalizeScore(score, kmerSize, rec.seq.length());
	return hit;
}

}
;

#endif /* SEQEVAL_H_ */
