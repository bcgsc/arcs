#ifndef _CONTIG_H_
#define _CONTIG_H_ 1

#include <string>
#include <utility>

typedef std::string ContigName;

enum Sense { SENSE_FORWARD = 0, SENSE_REVERSE, NUM_SENSE };
typedef std::pair<ContigName, Sense> ContigSense;


/** Hash functor for std::pair of two arbitrary types */
struct ContigSenseHash {
	std::size_t operator()(const ContigSense& contig) const {
		return std::hash<ContigName>()(contig.first)
			^ std::hash<int>()((int)contig.second);
	}
};

static inline std::string ContigSenseToString(const ContigSense& sense)
{
	std::string s(sense.first);
	s += sense.second == SENSE_FORWARD ? "+" : "-";
	return s;
}

#endif
