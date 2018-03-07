#ifndef _SCAFF_SIZES_H_
#define _SCAFF_SIZES_H_ 1

#include <string>
#include <utility>
#include <unordered_map>
#include <vector>

/**
 * a list of the input scaffolds and their lengths, in the order
 * that they appear in the input contigs FASTA file
 */
typedef std::pair<std::string, int> ScaffSizeRec;
typedef std::vector<ScaffSizeRec> ScaffSizeList;
typedef typename ScaffSizeList::const_iterator ScaffSizeConstIt;

/**
 * a list of the input scaffolds and their lengths, in the order
 * that they appear in the input contigs FASTA file
 */
typedef std::unordered_map<std::string, int> ScaffSizeMap;

#endif
