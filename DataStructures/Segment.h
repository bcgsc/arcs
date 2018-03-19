#ifndef _SEGMENT_H_
#define _SEGMENT_H_ 1

#include "Common/SetUtil.h"
#include "DataStructures/Barcode.h"
#include "DataStructures/Contig.h"

#include <cassert>
#include <iterator>
#include <limits>
#include <string>
#include <utility>

#define MIDDLE_SEGMENT unsigned(-1)

typedef unsigned SegmentIndex;
typedef unsigned Position;

/** One segment of a contig. */
typedef std::pair<ContigName, SegmentIndex> Segment;

/** A list of segments */
typedef std::vector<Segment> SegmentList;
typedef typename SegmentList::const_iterator SegmentListConstIt;

/** Hash a Segment. */
struct HashSegment {
    size_t operator()(const Segment& key) const {
        return std::hash<ContigName>()(key.first)
            ^ std::hash<SegmentIndex>()(key.second);
    }
};

/** Perform calculations related to contig segments */
class SegmentCalc
{
public:

    SegmentCalc() : m_segmentSize(0) {}
    SegmentCalc(unsigned segmentSize) : m_segmentSize(segmentSize) {}

    /**
    * Return the index of the contig segment containing the given
    * (one-based) sequence position. If the sequence length
    * is not evenly divisible by the region length,
    * the "remainder" region placed in the middle of the sequence.
    * Return std::numeric_limits<unsigned>::max() if the given position
    * falls within the middle remainder region.
    */
    unsigned index(unsigned pos, unsigned l) const
    {
        /* input should be one-based position, as in SAM format */
        assert(pos > 0);
        assert(pos <= l);
        assert(m_segmentSize <= l / 2);

        /* translate to zero-based pos */
        pos--;

        if (l % m_segmentSize == 0) {
            /* sequence length is perfectly divisible by segment length */
            return pos / m_segmentSize;
        }

        unsigned index;
        if (pos < l / 2) {
            /* pos is in left half of seq */
            index = pos / m_segmentSize;
            if (index >= segmentsPerHalf(l)) {
                /* pos is in remainder middle seg */
                return MIDDLE_SEGMENT;
            }
        } else {
            /* pos is in right half of seq */
            index = (pos - remainder(l)) / m_segmentSize;
            if (index < segmentsPerHalf(l)) {
                /* pos is in middle remainder seg */
                return MIDDLE_SEGMENT;
            }
        }
        return index;
    }

    /**
     * Return the segment index range (first and last segment indices)
     * containing the giving (1-based) alignment coordinate range.
     */
    std::pair<std::pair<SegmentIndex, SegmentIndex>, bool>
        indexRange(Position start, Position end, unsigned l)
    {
        if (l < m_segmentSize * 2) {
            return std::make_pair(std::make_pair(0, 0), false);
        }

        std::pair<SegmentIndex, SegmentIndex> range =
            std::make_pair(index(start, l), index(end, l));

        if (range.first == MIDDLE_SEGMENT && range.second == MIDDLE_SEGMENT) {
            return std::make_pair(range, false);
        } else if (range.first == MIDDLE_SEGMENT) {
            range.first = segmentsPerHalf(l);
            return std::make_pair(range, true);
        } else if (range.second == MIDDLE_SEGMENT) {
            range.second = segmentsPerHalf(l) - 1;
        }

        return std::make_pair(range, true);
    }

    /**
     * Return the number of segments contained in each
     * half of the sequence.  Precondition: Sequence length
     * is not evenly divisible by segment length.
     */
    unsigned segmentsPerHalf(unsigned l) const
    {
        assert(m_segmentSize <= l / 2);
        assert(l % m_segmentSize > 0);
        return l / 2 / m_segmentSize;
    }

    /** Return the number of segments in a sequence of the given length */
    unsigned segments(unsigned l) const
    {
        assert(m_segmentSize <= l / 2);
        if (l % m_segmentSize == 0) {
            /* sequence length is perfectly divisible by segment length */
            return l / m_segmentSize;
        }
        return segmentsPerHalf(l) * 2;
    }

    /**
     * Return one-based start position of the segment
     * with the given index
     */
    unsigned start(unsigned l, unsigned index, bool rc=false) const
    {
        assert(m_segmentSize <= l / 2);

        unsigned _start;
        if (l % m_segmentSize == 0) {
            _start = index * m_segmentSize + 1;
        } else {
            unsigned segsPerHalf = segmentsPerHalf(l);
            if (index < segsPerHalf) {
                _start = index * m_segmentSize + 1;
            } else {
                unsigned rightIndex = index - segsPerHalf;
                _start = segsPerHalf * m_segmentSize
                    + remainder(l) + rightIndex * m_segmentSize + 1;
            }
        }

        if (rc)
            _start = l - _start - m_segmentSize + 2;

        return _start;
    }

    /**
     * Return length of remainder area in middle of sequence,
     * which is not coverage by any segment
     */
    unsigned remainder(unsigned l) const
    {
        assert(m_segmentSize <= l / 2);
        if (l % m_segmentSize == 0) {
            return 0;
        }
        assert(l > m_segmentSize * segmentsPerHalf(l) * 2);
        return l - m_segmentSize * segmentsPerHalf(l) * 2;
    }

protected:

    /** length of a contig segment in base pairs */
    unsigned m_segmentSize;
};


typedef std::pair<Segment, Segment> SegmentPair;

/**
 * Iterate over all segments pairs whose distances are integer
 * multiples of the segment length.
 */
class SegmentPairIterator
    : public std::iterator<std::input_iterator_tag, const SegmentPair>
{
public:

    SegmentPairIterator()
        : m_length(0), m_segmentSize(0), m_divisible(false),
        m_pair(std::make_pair("", std::numeric_limits<unsigned>::max()),
            std::make_pair("", std::numeric_limits<unsigned>::max()))
    {}

    SegmentPairIterator(const std::string& id,
        unsigned length, unsigned segmentSize)
    	: m_length(length), m_segmentSize(segmentSize),
        m_divisible(false), m_calc(segmentSize),
        m_pair(std::make_pair(id, std::numeric_limits<unsigned>::max()),
            std::make_pair(id, std::numeric_limits<unsigned>::max()))
    {
        if (m_length < 2 * m_segmentSize) {
            return;
        }

        unsigned& i = m_pair.first.second;
        unsigned& j = m_pair.second.second;

        m_divisible = m_length % m_segmentSize == 0;

        if (m_divisible) {
            i = 0;
            j = 1;
            return;
        }

        if (m_calc.segmentsPerHalf(m_length) < 2) {
            return;
        }

        i = 0;
        j = 1;
    }

    const SegmentPair& operator *() const
    {
        return m_pair;
    }

    const SegmentPair* operator ->() const
    {
        return &m_pair;
    }

    bool operator ==(const SegmentPairIterator& it) const
    {
        const unsigned& i1 = m_pair.first.second;
        const unsigned& j1 = m_pair.second.second;
        const unsigned& i2 = it.m_pair.first.second;
        const unsigned& j2 = it.m_pair.second.second;
        return i1 == i2 && j1 == j2;
    }

    bool operator !=(const SegmentPairIterator& it) const
    {
        return !this->operator==(it);
    }

    SegmentPairIterator& operator ++()
    {
        next();
        return *this;
    }

    SegmentPairIterator operator ++(int)
    {
        SegmentPairIterator it = *this;
        next();
        return it;
    }

protected:

    /** advance to the next segment pair */
    void next()
    {
        unsigned& i = m_pair.first.second;
        unsigned& j = m_pair.second.second;

        if (m_divisible) {
            /*
             * special case: sequence length is perfectly
             * divisible by segment length
             */
            unsigned segments = m_calc.segments(m_length);
            for (; i < segments - 1; ++i, j=i) {
                for (++j; j < segments;) {
                    return;
                }
            }
        } else {
            unsigned segsPerHalf = m_calc.segmentsPerHalf(m_length);
            /* segment pairs in left half of seq */
            for (; i < segsPerHalf - 1; ++i, j=i) {
                for (++j; j < segsPerHalf;) {
                    return;
                }
            }
            if (i < segsPerHalf) {
                ++i;
                j=i;
            }
            /* segment pairs in right half of seq */
            for (; i < 2 * segsPerHalf - 1; ++i, j=i) {
                for (++j; j < 2 * segsPerHalf;) {
                    return;
                }
            }
        }

        i = std::numeric_limits<unsigned>::max();
        j = std::numeric_limits<unsigned>::max();
    }

    /** sequence length */
    const unsigned m_length;
    /** segment length */
    const unsigned m_segmentSize;
    /** true if sequence is perfectly divisible by segment length */
    bool m_divisible;
    /** performs segment-related calculations */
    SegmentCalc m_calc;
    /** current segment pair */
    SegmentPair m_pair;
};

#endif
