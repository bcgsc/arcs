#ifndef _EDGE_LIST_H_
#define _EDGE_LIST_H_ 1

#include "Common/Fisher.h"
#include "Common/PairHash.h"
#include "DataStructures/BarcodeToSegment.h"
#include "DataStructures/SharedBarcodeMap.h"
#include "DataStructures/Segment.h"
#include "DataStructures/SegmentToBarcode.h"
#include "DistanceEst/DistanceModel.h"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <cassert>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

enum BuildEdgeResultCode {
	ER_SUCCESS = 0,
	ER_NO_SOLUTION,
	ER_AMBIGUOUS
};

static inline std::string BuildEdgeResultToString(BuildEdgeResultCode code)
{
	switch(code) {
	case ER_SUCCESS: return "SUCCESS";
	case ER_NO_SOLUTION: return "NO_SOLUTION";
	case ER_AMBIGUOUS: return "AMBIGUOUS";
	}
	assert(false);
}

struct BuildEdgeParams
{
	unsigned k;
	unsigned segmentSize;
	unsigned maxDist;
	double minp;

	BuildEdgeParams() : k(0), segmentSize(0), maxDist(0), minp(0) {}
};

struct BuildEdgeResult
{
	ContigSense contig1;
	ContigSense contig2;
	unsigned dist;
	double p;

	BuildEdgeResultCode code;

	BuildEdgeResult() : dist(0), p(0), code(ER_NO_SOLUTION) {}
};

typedef std::vector<BuildEdgeResult> EdgeList;
typedef typename EdgeList::const_iterator EdgeListConstIt;

struct BuildEdgeCounters
{
	size_t processed;
	size_t success;
	size_t noSolution;
	size_t ambiguous;

	BuildEdgeCounters()
		: processed(0), success(0), noSolution(0), ambiguous(0) {}

	friend std::ostream& operator <<(std::ostream& out,
		const BuildEdgeCounters& o)
	{
		out << "processed: " << o.processed << ", "
			<< "success: " << o.success << ", "
			<< "no solution: " << o.noSolution << ", "
			<< "ambiguous: " << o.ambiguous << "\n";
		return out;
	}
};

struct BuildEdgeStreams
{
	std::ostream& segmentStatsOut;
	std::ostream& fisherStatsOut;
	std::ostream& pairStatsOut;

	BuildEdgeStreams(
		std::ostream& segmentStatsOut,
		std::ostream& fisherStatsOut,
		std::ostream& pairStatsOut)
		: segmentStatsOut(segmentStatsOut),
		fisherStatsOut(fisherStatsOut),
		pairStatsOut(pairStatsOut)
		{}
};

static inline BuildEdgeResult buildEdge(
	const ContigName& name1, unsigned length1,
	const ContigName& name2, unsigned length2,
	const SegmentToBarcode& segmentToBarcode,
	const DistanceModel& distModel,
	const BuildEdgeParams& params,
	BuildEdgeStreams& streams)
{
	assert(params.segmentSize >= params.k);
	assert(params.maxDist >= params.segmentSize);

	BuildEdgeResult result;
	result.p = 0.0;
	result.code = ER_NO_SOLUTION;
	bool ambiguous = false;

	for (unsigned dist = params.segmentSize; dist <= params.maxDist;
		dist += params.segmentSize) {
		for (int sense1 = SENSE_FORWARD; sense1 < NUM_SENSE; ++sense1) {
			for (int sense2 = SENSE_FORWARD; sense2 < NUM_SENSE; ++sense2) {
				ContigSense contig1(name1, (Sense)sense1);
				ContigSense contig2(name2, (Sense)sense2);
				double p = fisher(contig1, length1, contig2, length2,
					segmentToBarcode, distModel, dist, streams.segmentStatsOut);
				streams.fisherStatsOut
					<< ContigSenseToString(contig1) << '\t'
					<< length1 << '\t'
					<< ContigSenseToString(contig2) << '\t'
					<< length2 << '\t'
					<< dist << '\t'
					<< p << '\n';
				if (p > params.minp) {
					if (result.p > params.minp)
						ambiguous = true;
					result.code = ER_SUCCESS;
					result.contig1 = contig1;
					result.contig2 = contig2;
					result.dist = dist - params.segmentSize + params.k - 1;
					result.p = p;
				}
			}
		}
	}

	if (ambiguous) {
		result.code = ER_AMBIGUOUS;
		result.p = 0.0;
	}

	return result;
}

static inline void buildEdgeList(
	const SegmentToBarcode& segmentToBarcode,
	const SharedBarcodeMap& sharedBarcodeMap,
	const DistanceModel& distModel,
	const ScaffSizeMap& scaffSizeMap,
	const BuildEdgeParams& params, EdgeList& edgeList,
	BuildEdgeStreams& streams)
{
	/* log file headers */

	streams.segmentStatsOut
		<< "id1" << '\t'
		<< "l1" << '\t'
		<< "id2" << '\t'
		<< "l2" << '\t'
		<< "gap" << '\t'
		<< "segment1" << '\t'
		<< "start1" << '\t'
		<< "end1" << '\t'
		<< "segment2" << '\t'
		<< "start2" << '\t'
		<< "end2" << '\t'
		<< "dist" << '\t'
		<< "rounded_dist" << '\t'
		<< "jaccard" << '\t'
		<< "p" << '\n';

	streams.fisherStatsOut
		<< "v1" << '\t'
		<< "l1" << '\t'
		<< "v2" << '\t'
		<< "l2" << '\t'
		<< "dist" << '\t'
		<< "p" << '\n';

	streams.pairStatsOut
		<< "v1" << '\t'
		<< "v2" << '\t'
		<< "result" << '\n';

	/* contig ends already used in an edge */
	std::unordered_set<ContigSense, ContigSenseHash> paired;

	BuildEdgeCounters edgeCounts;
	const size_t edgeProgressStep = 100;
	BuildEdgeCounters vertexCounts;
	const size_t vertexProgressStep = 100;

	for (SharedBarcodeMapConstIt it1 = sharedBarcodeMap.begin();
		it1 != sharedBarcodeMap.end(); ++it1, vertexCounts.processed++)
	{
		if (vertexCounts.processed % vertexProgressStep == 0)
			std::cerr << "vertex stats: " << edgeCounts;

		const ContigName& name1 = it1->first;
		unsigned length1 = scaffSizeMap.at(name1);

		/*
		 * find closest left/right neighbour contigs,
		 * along with their distance and orientation
		 */

		BuildEdgeResult closest[NUM_SENSE];
		/*
		 * if there are multiple closest neighbours on the left/right
		 * of contig1 with equal distances
		 */
		bool ambiguous[NUM_SENSE] = { false, false };

		const ContigToCount& contigToCount = it1->second;
		for (ContigToCountConstIt it2 = contigToCount.begin();
			it2 != contigToCount.end(); ++it2, edgeCounts.processed++)
		{
			if (edgeCounts.processed % edgeProgressStep == 0)
				std::cerr << "edge stats: " << edgeCounts;

			const ContigName& name2 = it2->first;
			unsigned length2 = scaffSizeMap.at(name2);

			BuildEdgeResult result = buildEdge(
				name1, length1, name2, length2, segmentToBarcode,
				distModel, params, streams);

			switch(result.code) {
			case ER_AMBIGUOUS: edgeCounts.ambiguous++; continue;
			case ER_NO_SOLUTION: edgeCounts.noSolution++; continue;
			case ER_SUCCESS: edgeCounts.success++; break;
			}

			assert(result.code == ER_SUCCESS);
			assert(result.p >= params.minp);

			/* if either contig end already has a mate */
			if (paired.find(result.contig1) != paired.end()
				|| paired.find(result.contig2) != paired.end())
				continue;

			Sense sense1 = result.contig1.second;
			if (closest[sense1].p < params.minp
				|| result.dist < closest[sense1].dist) {
				closest[sense1] = result;
			} else if (result.dist == closest[sense1].dist) {
				ambiguous[sense1] = true;
			}

		}

		/* mark contig1 as visited */

		ContigSense contig1Fwd(name1, SENSE_FORWARD);
		ContigSense contig1Rev(name1, SENSE_REVERSE);
		paired.insert(contig1Fwd);
		paired.insert(contig1Rev);

		/* add edges to list and mark contig 2 as paired */

		for (int sense = SENSE_FORWARD; sense < NUM_SENSE; ++sense)
		{
			Sense sense1 = (Sense)sense;
			BuildEdgeResult& _closest = closest[sense1];

			if (ambiguous[sense1]) {
				vertexCounts.ambiguous++;
				continue;
			}

			switch (_closest.code) {
			case ER_AMBIGUOUS:
				vertexCounts.noSolution++;
				streams.pairStatsOut
					<< ContigSenseToString(ContigSense(name1, sense1)) << '\t'
					<< "NA" << '\t'
					<< BuildEdgeResultToString(ER_AMBIGUOUS) << '\n';
				continue;
			case ER_NO_SOLUTION:
				vertexCounts.noSolution++;
				streams.pairStatsOut
					<< ContigSenseToString(ContigSense(name1, sense1)) << '\t'
					<< "NA" << '\t'
					<< BuildEdgeResultToString(ER_NO_SOLUTION) << '\n';
				continue;
			case ER_SUCCESS: vertexCounts.success++;
				streams.pairStatsOut
					<< ContigSenseToString(_closest.contig1) << '\t'
					<< ContigSenseToString(_closest.contig2) << '\t'
					<< BuildEdgeResultToString(ER_SUCCESS) << '\n';
				break;
			}

			/* add edge to list */

			assert(_closest.contig1.second == sense1);
			assert(_closest.code == ER_SUCCESS && !ambiguous[sense1]);
			edgeList.push_back(_closest);

			/* mark contig2 as paired */

			ContigSense contig2 = _closest.contig2;
			contig2.second = sense1 == SENSE_FORWARD ?
				contig2.second == SENSE_FORWARD ? SENSE_REVERSE : SENSE_FORWARD :
				contig2.second == SENSE_FORWARD ? SENSE_FORWARD : SENSE_REVERSE;
			paired.insert(contig2);
		}
	}

	std::cerr << "edge stats: " << edgeCounts;
	std::cerr << "vertex stats: " << vertexCounts;
}

#endif
