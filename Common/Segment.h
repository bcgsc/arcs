#ifndef _SEGMENT_H_
#define _SEGMENT_H_ 1

#include <cassert>

/** One segment of a contig. */
typedef std::pair<std::string, unsigned> Segment;

/** Hash a Segment. */
struct HashSegment {
    size_t operator()(const Segment& key) const {
        return std::hash<std::string>()(key.first)
            ^ std::hash<unsigned>()(key.second);
    }
};

/**
 * Return the index of the contig segment containing the given
 * (one-based) sequence position. If the sequence length
 * is not evenly divisible by the region length,
 * the "remainder" region placed in the middle of the sequence.
 * Return std::numeric_limits<unsigned>::max() if the given position
 * falls within the middle remainder region.
 */
static inline unsigned posToSegmentIndex(
    unsigned pos, unsigned seqLen, unsigned segLen)
{
    /* input should be one-based position, as in SAM format */

    assert(pos > 0);
    assert(pos <= seqLen);

    /* translate to zero-based pos */
    pos--;

    if (seqLen % segLen == 0) {
        /* sequence length is perfectly divisible by segment length */
        return pos / segLen;
    }

    unsigned index;
    unsigned segsPerHalf = seqLen / 2 / segLen;
    if (pos < seqLen / 2) {
        /* pos is in left half of seq */
        index = pos / segLen;
        if (index >= segsPerHalf) {
            /* pos is in remainder middle seg */
            index = std::numeric_limits<unsigned>::max();
        }
    } else {
        /* pos is in right half of seq */
        unsigned middleLen = seqLen - segLen * segsPerHalf * 2;
        index = (pos - middleLen) / segLen;
        if (index < segsPerHalf) {
            /* pos is in middle remainder seg */
            index = std::numeric_limits<unsigned>::max();
        }
    }
    return index;
}

#endif
