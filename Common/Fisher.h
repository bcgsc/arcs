#ifndef _ARCS_FISHER_H_
#define _ARCS_FISHER_H_ 1

#include "DataStructures/Contig.h"
#include "DataStructures/Segment.h"
#include "DataStructures/SegmentToBarcode.h"
#include "DistanceEst/DistanceModel.h"

#include <boost/math/special_functions/gamma.hpp>
#include <cmath>
#include <iostream>
#include <limits>

/**
 * Compute a p-value that two contigs are located next
 * to each other with the given orientations and distance.
 * The p-value is calculated using a Fisher combined
 * probability test for all segment pairs across
 * contig1 and contig2.
 */
static inline double fisher(
	const ContigSense& contig1, unsigned length1,
	const ContigSense& contig2, unsigned length2,
	const SegmentToBarcode& segmentToBarcode,
	const DistanceModel& distModel, int dist,
	std::ostream& segmentStatsOut)
{
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

			unsigned start1 = calc.start(length1, i,
				contig1.second == SENSE_REVERSE);
			unsigned start2 = calc.start(length2, j,
				contig2.second == SENSE_REVERSE);
			unsigned end1 = start1 + segmentSize - 1;
			unsigned end2 = start2 + segmentSize - 1;

			assert(dist > -int(length1));
			unsigned segmentDist = length1 + dist + start2 - start1;

			/* round distance to a multiple of segment size */

			unsigned roundedSegmentDist =
				(unsigned)std::round((double)segmentDist / segmentSize)
				* segmentSize;

			if (roundedSegmentDist < segmentSize)
				continue;

			/*
			 * if there are insufficient Jaccard samples
			 * for given distance
			 */

			if (distModel.samples(roundedSegmentDist) < MIN_SAMPLES)
				continue;

			/* calculate p-value for segment pair at given distance */

			Segment segment1(contig1.first, i);
			Segment segment2(contig2.first, j);
			double jacc = jaccard(segment1, segment2, segmentToBarcode);
			double p = distModel.p(roundedSegmentDist, jacc);

			segmentStatsOut
				<< contig1.first << '\t'  /* id1 */
				<< length1 << '\t'        /* length1 */
				<< contig2.first << '\t'  /* id2 */
				<< length2 << '\t'     /* l2 */
				<< dist << '\t'        /* path_length */
				<< i << '\t'           /* segment1 */
				<< start1 << '\t'      /* start1 */
				<< end1 << '\t'        /* end1 */
				<< j << '\t'           /* segment2 */
				<< start2 << '\t'      /* start2 */
				<< end2 << '\t'        /* end2 */
				<< segmentDist << '\t' /* dist */
				<< roundedSegmentDist << '\t' /* roundedDist */
				<< jacc << '\t'        /* jaccard */
				<< p << '\n';          /* p */

			/* add to sum for Fisher combined probability test */

			if (p == 0)
				overflow = true;

			logsum += log(p);
			segmentPairs++;

		}
	}

	if (segmentPairs == 0 || overflow)
		return 0.0;

	/* combine p-values using *Fisher's mthod* */
	double fisherp = 1.0 -
		boost::math::gamma_p(segmentPairs, -logsum / 2.0);

	/* multiple testing correction */
	double corrected = double(segmentPairs + 1) / (2 * segmentPairs)
		* fisherp;

	return corrected;
}

#endif
