#ifndef _ARCS_FISHER_H_
#define _ARCS_FISHER_H_ 1

#include "Common/ContigNode.h"
#include "DataStructures/Contig.h"
#include "DataStructures/Segment.h"
#include "DataStructures/SegmentToBarcode.h"
#include "DistanceEst/DistanceModel.h"

#include <boost/math/special_functions/gamma.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

enum SegmentPairResultCode {
	SPRC_COMPUTED_P_VALUE=0,
	SPRC_SEGMENTS_TOO_CLOSE,
	SPRC_NOT_ENOUGH_SAMPLES,
	SPRC_OVERFLOW
};

template <class Graph>
struct SegmentPairResult
{
    ContigNode contig1;
	ContigNode contig2;
	const Graph& g;
	SegmentIndex segmentIndex1;
	SegmentIndex segmentIndex2;
	unsigned segmentSize;
	int distance;
	float jaccard;
	float p;
	SegmentPairResultCode code;

	SegmentPairResult(const ContigNode& contig1, const ContigNode& contig2,
		const Graph& g, SegmentIndex segmentIndex1, SegmentIndex segmentIndex2,
		unsigned segmentSize, int distance, float jaccard, float p,
		SegmentPairResultCode code)
		: contig1(contig1), contig2(contig2), g(g),
		segmentIndex1(segmentIndex1), segmentIndex2(segmentIndex2),
		segmentSize(segmentSize), distance(distance), jaccard(jaccard),
		p(p), code(code)
		{}

	friend std::ostream& operator <<(std::ostream& out, const SegmentPairResult& o)
	{
		unsigned length1 = o.g[o.contig1].length;
		unsigned length2 = o.g[o.contig2].length;

		std::stringstream name1;
		name1 << get(vertex_name, o.g, o.contig1);
		assert(name1);
		std::stringstream name2;
		name2 << get(vertex_name, o.g, o.contig2);
		assert(name2);

		/*
		 * calc start/end coordinates of the segments w.r.t. to the given
		 * contig orientations
		 */

		SegmentCalc calc(o.segmentSize);
		unsigned start1 = calc.start(length1, o.segmentIndex1, o.contig1.sense());
		unsigned start2 = calc.start(length2, o.segmentIndex2, o.contig2.sense());
		unsigned end1 = start1 + o.segmentSize - 1;
		unsigned end2 = start2 + o.segmentSize - 1;

		/* calc segment distance rounded to a multiple of segment size */

		unsigned roundedDist =
			(unsigned)std::round((double)o.distance / o.segmentSize)
			* o.segmentSize;

		/* calc end-to-end distance between the two contigs */

		int endDist = o.distance - start2 - (length1 - start1);
		assert(endDist > -int(length1));

		out << name1.str() << '\t' /* id1 */
			<< length1 << '\t'     /* length1 */
			<< name2.str() << '\t' /* id2 */
			<< length2 << '\t'     /* l2 */
			<< endDist << '\t'     /* end-to-end distance between sequences */
			<< o.segmentIndex1 << '\t' /* segment index w.r.t. forward orientation */
			<< start1 << '\t'      /* start1 */
			<< end1 << '\t'        /* end1 */
			<< o.segmentIndex2 << '\t' /* segment index w.r.t. forward orientation */
			<< start2 << '\t'      /* start2 */
			<< end2 << '\t'        /* end2 */
			<< o.distance << '\t'  /* dist */
			<< roundedDist << '\t' /* roundedDist */
			<< o.jaccard << '\t'   /* jaccard */
			<< o.p << '\n';        /* p */

		return out;
	}

};

enum FisherResultCode {
	FRC_COMPUTED_P_VALUE=0,
	FRC_OVERFLOW
};

template <class Graph>
struct FisherResult
{
	ContigNode contig1;
	ContigNode contig2;
	int distance;
	float p;
	FisherResultCode code;

	std::vector< SegmentPairResult<Graph> > segmentPairResults;

	FisherResult() : distance(0), p(0.0f), code(FRC_OVERFLOW) {}
};

/**
 * Compute a p-value that two contigs are located next
 * to each other with the given orientations and distance.
 * The p-value is calculated using a Fisher combined
 * probability test for all segment pairs across
 * contig1 and contig2.
 */
template <class Graph>
static inline FisherResult<Graph> fisher(
	const ContigNode& contig1, const ContigNode& contig2, const Graph& g,
	const SegmentToBarcode& segmentToBarcode,
	const DistanceModel& distModel, int dist)
{
	FisherResult<Graph> result;
	result.contig1 = contig1;
	result.contig2 = contig2;
	result.distance = dist;
	result.code = FRC_COMPUTED_P_VALUE;

	unsigned length1 = g[contig1].length;
	unsigned length2 = g[contig2].length;

	const unsigned segmentSize = distModel.segmentSize();

	SegmentCalc calc(segmentSize);

	unsigned segments1 = calc.segments(length1);
	unsigned segments2 = calc.segments(length2);
	assert(segments1 >= 1 && segments2 >= 1);

	/* compute p-values for segment pairs */

	unsigned MIN_SAMPLES = 300;
	/* running sum for Fisher combined probability test */
	double logsum = 0.0;
	unsigned segmentPairs = 0;
	bool overflow = false;
	for (unsigned i = 0; i < segments1 && !overflow; ++i) {
		for (unsigned j = 0; j < segments2 && !overflow; ++j) {

			/* compute distance between segments */

			unsigned start1 = calc.start(length1, i, contig1.sense());
			unsigned start2 = calc.start(length2, j, contig2.sense());

			assert(dist > -int(length1));
			unsigned segmentDist = length1 + dist + start2 - start1;

			/* round distance to a multiple of segment size */

			unsigned roundedSegmentDist =
				(unsigned)std::round((double)segmentDist / segmentSize)
				* segmentSize;

			if (roundedSegmentDist < segmentSize) {
				result.segmentPairResults.push_back(
					SegmentPairResult<Graph>(contig1, contig2, g, i, j, 
					segmentSize, segmentDist, 0.0f, 0.0f,
					SPRC_SEGMENTS_TOO_CLOSE));
				continue;
			}

			/*
			 * if there are insufficient Jaccard samples
			 * for given distance
			 */

			if (distModel.samples(roundedSegmentDist) < MIN_SAMPLES) {
				result.segmentPairResults.push_back(
					SegmentPairResult<Graph>(contig1, contig2, g, i, j,
						segmentSize, segmentDist, 0.0f, 0.0f,
						SPRC_NOT_ENOUGH_SAMPLES));
				continue;
			}

			/* calculate p-value for segment pair at given distance */

			Segment segment1(contig1.id(), i);
			Segment segment2(contig2.id(), j);
			double jacc = jaccard(segment1, segment2, segmentToBarcode);
			double p = distModel.p(roundedSegmentDist, jacc);

			SegmentPairResult<Graph> pairResult(contig1, contig2, g, i, j,
				segmentSize, segmentDist, jacc, p, SPRC_COMPUTED_P_VALUE);

			/* add to sum for Fisher combined probability test */

			if (p == 0) {
				pairResult.code = SPRC_OVERFLOW;
				result.code = FRC_OVERFLOW;
			}

			logsum += log(p);

			result.segmentPairResults.push_back(pairResult);
			segmentPairs++;
		}
	}

	if (segmentPairs == 0 || result.code == FRC_OVERFLOW)
		return result;

	/* combine p-values using *Fisher's method* */
	double fisherp = 1.0 -
		boost::math::gamma_p(segmentPairs, -logsum / 2.0);

	/* multiple testing correction */
	double corrected = double(segmentPairs + 1) / (2 * segmentPairs)
		* fisherp;

	result.p = corrected;
	result.code = FRC_COMPUTED_P_VALUE;

	return result;
}

#endif
