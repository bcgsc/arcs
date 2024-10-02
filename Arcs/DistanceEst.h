#ifndef _DISTANCE_EST_H_
#define _DISTANCE_EST_H_ 1

#include "Arcs/Arcs.h"
#include "Common/MapUtil.h"
#include "Common/PairHash.h"
#include "Common/StatUtil.h"
#include <array>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>

/** min/max distance estimate for a pair contigs */
struct DistanceEstimate
{
	int minDist;
	int dist;
	int maxDist;
	double jaccard;

	DistanceEstimate()
	  : minDist(0)
	  , dist(0)
	  , maxDist(0)
	  , jaccard(0.0)
	{
	}
};

/**
 * Records the distance between the head/tail regions of the same
 * contig vs. barcode union size, barcode intersection size,
 * and number of distinct barcodes mapped to each end.
 */
struct DistSample
{
	unsigned distance;
	unsigned barcodesHead;
	unsigned barcodesTail;
	unsigned barcodesUnion;
	unsigned barcodesIntersect;

	DistSample()
	  : distance(std::numeric_limits<unsigned>::max())
	  , barcodesHead(0)
	  , barcodesTail(0)
	  , barcodesUnion(0)
	  , barcodesIntersect(0)
	{
	}
};

/** maps contig ID => intra-contig distance/barcode sample */
typedef std::unordered_map<std::string, DistSample> DistSampleMap;
typedef typename DistSampleMap::const_iterator DistSampleConstIt;

/** maps barcode Jaccard index => intra-contig distance sample */
typedef std::map<double, DistSample> JaccardToDist;
typedef typename JaccardToDist::const_iterator JaccardToDistConstIt;

/** Barcode stats for a candidate pair of contig ends */
struct BarcodeStats
{
	unsigned barcodes1;
	unsigned barcodes2;
	unsigned barcodesUnion;
	unsigned barcodesIntersect;

	BarcodeStats()
	  : barcodes1(0)
	  , barcodes2(0)
	  , barcodesUnion(0)
	  , barcodesIntersect(0)
	{
	}
};

/** possible contig pair orientations (e.g. HH = head-to-head) */
enum PairOrientation
{
	HH = 0,
	HT,
	TH,
	TT,
	NUM_ORIENTATIONS
};

/** barcode stats for each possible contig pair orientation (HH,HT,TH,HH) */
typedef std::array<BarcodeStats, NUM_ORIENTATIONS> BarcodeStatsArray;

/** barcode stats each possible orientation of a contig pair */
typedef std::map<ARCS::ContigPair, BarcodeStatsArray> PairToBarcodeStats;
typedef typename PairToBarcodeStats::iterator PairToBarcodeStatsIt;

/**
 * Measure distance between contig ends vs.
 * barcode intersection size and barcode union size.
 */
void
calcDistSamples(
    const ARCS::IndexMap& imap,
    const ARCS::ContigToLength& contigToLength,
    const std::unordered_map<std::string, int>& indexMultMap,
    const ARCS::ArcsParams& params,
    DistSampleMap& distSamples)
{
	/* for each chromium barcode */
	for (auto barcodeIt = imap.begin(); barcodeIt != imap.end(); ++barcodeIt) {
		/* skip barcodes outside of min/max multiplicity range */
		std::string index = barcodeIt->first;
		int indexMult = indexMultMap.at(index);
		if (indexMult < params.min_mult || indexMult > params.max_mult)
			continue;

		/* contig head/tail => number of mapped read pairs */
		const ARCS::ScafMap& contigToCount = barcodeIt->second;

		for (auto contigIt = contigToCount.begin(); contigIt != contigToCount.end(); ++contigIt) {
			std::string contigID;
			bool isHead;
			std::tie(contigID, isHead) = contigIt->first;
			int readPairs = contigIt->second;

			/*
			 * skip contigs with less than required number of
			 * mapped read pairs (-c option)
			 */
			if (readPairs < params.min_reads)
				continue;

			/*
			 * skip contigs shorter than 2 times the contig
			 * end length, because we want our distance samples
			 * to be based on a uniform head/tail length
			 */

			unsigned l = contigToLength.at(contigID);
			if (l < (unsigned)2 * params.end_length)
				continue;

			DistSample& distSample = distSamples[contigID];
			distSample.distance = l - 2 * params.end_length;

			if (isHead)
				distSample.barcodesHead++;
			else
				distSample.barcodesTail++;

			/*
			 * Check if barcode also maps to other end of contig
			 * with sufficient number of read pairs.
			 *
			 * The `isHead` part of the `if` condition prevents
			 * double-counting when a barcode maps to both
			 * ends of a contig.
			 */

			ARCS::CI otherEnd(contigID, !isHead);
			ARCS::ScafMapConstIt otherIt = contigToCount.find(otherEnd);
			bool foundOther = otherIt != contigToCount.end() && otherIt->second >= params.min_reads;

			if (foundOther && isHead) {
				distSample.barcodesIntersect++;
				distSample.barcodesUnion++;
			} else if (!foundOther) {
				distSample.barcodesUnion++;
			}
		}
	}
}

/**
 * Build a ordered map from barcode Jaccard index to
 * distance sample. Each distance sample comes from
 * measuring the distance between the head/tail of the
 * same contig, along with associated head/tail barcode
 * counts.
 */
static inline void
buildJaccardToDist(const DistSampleMap& distSamples, JaccardToDist& jaccardToDist)
{
	for (DistSampleConstIt it = distSamples.begin(); it != distSamples.end(); ++it) {
		const DistSample& sample = it->second;
		double jaccard = double(sample.barcodesIntersect) / sample.barcodesUnion;
		jaccardToDist.insert(JaccardToDist::value_type(jaccard, sample));
	}
}

/**
 * Check requirements for using the given barcode-to-contig-end mapping
 * in distance estimates. Return true if we should the given
 * mapping in our calculations.
 */
static inline bool
validBarcodeMapping(unsigned contigLength, int pairs, const ARCS::ArcsParams& params)
{
	/*
	 * skip contigs with less than required number of
	 * mapped read pairs (-c option)
	 */

	if (pairs < params.min_reads)
		return false;

	/*
	 * skip contigs shorter than 2 times the contig
	 * end length, our distance samples are based
	 * based on a uniform head/tail length
	 */

	if (contigLength < unsigned(2 * params.end_length))
		return false;

	return true;
}

/** calculate shared barcode stats for candidate contig pairs */
static inline void
buildPairToBarcodeStats(
    const ARCS::IndexMap& imap,
    const std::unordered_map<std::string, int>& indexMultMap,
    const ARCS::ContigToLength& contigToLength,
    const ARCS::ArcsParams& params,
    PairToBarcodeStats& pairToStats)
{
	typedef std::unordered_map<ARCS::CI, size_t, PairHash> ContigEndToBarcodeCount;
	typedef typename ContigEndToBarcodeCount::const_iterator BarcodeCountConstIt;
	ContigEndToBarcodeCount contigEndToBarcodeCount;

	/* calculate number of shared barcodes for candidate contig end pairs */

	for (auto barcodeIt = imap.begin(); barcodeIt != imap.end(); ++barcodeIt) {
		/* skip barcodes outside of min/max multiplicity range */
		std::string index = barcodeIt->first;
		int indexMult = indexMultMap.at(index);
		if (indexMult < params.min_mult || indexMult > params.max_mult)
			continue;

		/* contig head/tail => number of mapped read pairs */
		const ARCS::ScafMap& contigEndToPairCount = barcodeIt->second;

		for (auto endIt1 = contigEndToPairCount.begin(); endIt1 != contigEndToPairCount.end();
		     ++endIt1) {
			/* get contig ID and head/tail flag */
			std::string id1;
			bool head1;
			std::tie(id1, head1) = endIt1->first;

			/* check requirements for calculating distance estimates */
			unsigned length1 = contigToLength.at(id1);
			int pairs1 = endIt1->second;
			if (!validBarcodeMapping(length1, pairs1, params))
				continue;

			/* count distinct barcodes mapped to head/tail of each contig */
			contigEndToBarcodeCount[endIt1->first]++;

			for (auto endIt2 = contigEndToPairCount.begin(); endIt2 != contigEndToPairCount.end();
			     ++endIt2) {
				/* get contig ID and head/tail flag */
				std::string id2;
				bool head2;
				std::tie(id2, head2) = endIt2->first;

				/* check requirements for calculating distance estimates */
				unsigned length2 = contigToLength.at(endIt2->first.first);
				int pairs2 = endIt2->second;
				if (!validBarcodeMapping(length2, pairs2, params))
					continue;

				/* avoid double-counting contig end pairs */
				if (id1 > id2)
					continue;

				/* initialize barcode/weight data for contig end pair */
				ARCS::ContigPair pair(id1, id2);
				if (pairToStats.count(pair) == 0)
					pairToStats[pair].fill(BarcodeStats());

				// Head - Head
				if (head1 && head2) {
					pairToStats[pair][0].barcodesIntersect++;
					// Head - Tail
				} else if (head1 && !head2) {
					pairToStats[pair][1].barcodesIntersect++;
					// Tail - Head
				} else if (!head1 && head2) {
					pairToStats[pair][2].barcodesIntersect++;
					// Tail - Tail
				} else if (!head1 && !head2) {
					pairToStats[pair][3].barcodesIntersect++;
				}
			}
		}
	}

	/*
	 * Compute/store further barcode stats for each candidate
	 * contig pair:
	 *
	 * (1) number of distinct barcodes mapping to contig A (|A|)
	 * (2) number of distinct barcodes mapping to contig B (|B|)
	 * (3) barcode union size for contigs A and B (|A union B|)
	 */

	for (PairToBarcodeStatsIt it = pairToStats.begin(); it != pairToStats.end(); ++it) {
		for (PairOrientation i = HH; i < NUM_ORIENTATIONS; i = PairOrientation(i + 1)) {
			BarcodeStats& stats = it->second.at(i);

			const std::string& id1 = it->first.first;
			const std::string& id2 = it->first.second;

			ARCS::CI tail1(id1, i == HH || i == HT);
			ARCS::CI tail2(id2, i == HH || i == TH);

			BarcodeCountConstIt countIt1 = contigEndToBarcodeCount.find(tail1);
			if (countIt1 == contigEndToBarcodeCount.end())
				continue;
			stats.barcodes1 = countIt1->second;
			assert(stats.barcodes1 > 0);

			BarcodeCountConstIt countIt2 = contigEndToBarcodeCount.find(tail2);
			if (countIt2 == contigEndToBarcodeCount.end())
				continue;
			stats.barcodes2 = countIt2->second;
			assert(stats.barcodes2 > 0);

			assert(stats.barcodes1 + stats.barcodes2 >= stats.barcodesIntersect);
			stats.barcodesUnion = stats.barcodes1 + stats.barcodes2 - stats.barcodesIntersect;
		}
	}
}

/** estimate min/max distance between a pair of contigs */
std::pair<DistanceEstimate, bool>
estimateDistance(
    const BarcodeStats& stats,
    const JaccardToDist& jaccardToDist,
    const ARCS::ArcsParams& params)
{
	DistanceEstimate result;

	/*
	 * if distance estimation was not enabled (`-D`) or input contigs
	 * were too short to provide any training data
	 */

	if (jaccardToDist.empty())
		return std::make_pair(result, false);

	/*
	 * barcodesUnion == 0 when a pair doesn't
	 * meet the requirements for distance estimation,
	 * (e.g. contig length < 2 * params.end_length)
	 */

	if (stats.barcodesUnion == 0)
		return std::make_pair(result, false);

	/* calc jaccard score for current contig pair */

	result.jaccard = double(stats.barcodesIntersect) / stats.barcodesUnion;
	assert(result.jaccard >= 0.0 && result.jaccard <= 1.0);

	/*
	 * get intra-contig distance samples with
	 * with closest Jaccard scores
	 */

	JaccardToDistConstIt lowerIt, upperIt;
	std::tie(lowerIt, upperIt) = closestKeys(jaccardToDist, result.jaccard, params.dist_bin_size);

	std::vector<unsigned> distances;
	for (JaccardToDistConstIt sampleIt = lowerIt; sampleIt != upperIt; ++sampleIt) {
		distances.push_back(sampleIt->second.distance);
	}

	std::sort(distances.begin(), distances.end());

	/* use 1st percentile, median, and 99th percentile */

	result.minDist = (int)floor(quantile(distances.begin(), distances.end(), 0.01));
	result.dist = (int)round(quantile(distances.begin(), distances.end(), 0.5));
	result.maxDist = (int)ceil(quantile(distances.begin(), distances.end(), 0.99));

	return std::make_pair(result, true);
}

/** add distance estimates to output graph edges */
static inline void
addEdgeDistances(
    const PairToBarcodeStats& pairToStats,
    const JaccardToDist& jaccardToDist,
    const ARCS::ArcsParams& params,
    ARCS::Graph& g)
{
	if (jaccardToDist.empty())
		return;

	for (const auto e : boost::make_iterator_range(boost::edges(g))) {

		auto v1 = source(e, g);
		auto v2 = target(e, g);

		auto id1 = g[v1].id;
		auto id2 = g[v2].id;

		auto orientation = g[e].orientation;

		auto pair = std::make_pair(id1, id2);
		auto statsIt = pairToStats.find(pair);
		if (statsIt == pairToStats.end())
			continue;
		const BarcodeStats& stats = statsIt->second.at(orientation);

		DistanceEstimate est;
		bool success;

		std::tie(est, success) = estimateDistance(stats, jaccardToDist, params);
		if (!success)
			continue;

		g[e].minDist = est.minDist;
		g[e].dist = est.dist;
		g[e].maxDist = est.maxDist;
		g[e].jaccard = est.jaccard;
	}
}

/** dump distance estimates and barcode data to TSV */
static inline void
writeDistTSV(const std::string& path, const PairToBarcodeStats& pairToStats, const ARCS::Graph& g)
{
	assert(!path.empty());

	/* open output TSV file */

	ofstream tsvOut;
	tsvOut.open(path.c_str());

	/* write TSV headers */

	tsvOut << "contig1" << '\t' << "contig2" << '\t' << "min_dist" << '\t' << "dist" << '\t'
	       << "max_dist" << '\t' << "barcodes1" << '\t' << "barcodes2" << '\t' << "barcodes_union"
	       << '\t' << "barcodes_intersect" << '\n';
	assert(tsvOut);

	for (const auto e : boost::make_iterator_range(boost::edges(g))) {

		auto v1 = source(e, g);
		auto v2 = target(e, g);

		auto id1 = g[v1].id;
		auto id2 = g[v2].id;

		auto orientation = g[e].orientation;

		auto pair = std::make_pair(id1, id2);
		auto statsIt = pairToStats.find(pair);
		if (statsIt == pairToStats.end())
			continue;
		const BarcodeStats& stats = statsIt->second.at(orientation);

		bool sense1 = orientation < 2;
		bool sense2 = orientation % 2;

		tsvOut << pair.first << (sense1 ? '-' : '+') << '\t' << pair.second << (sense2 ? '-' : '+')
		       << '\t';
		if (g[e].jaccard >= 0) {
			tsvOut << g[e].minDist << '\t' << g[e].dist << '\t' << g[e].maxDist << '\t';
		} else {
			tsvOut << "NA" << '\t' << "NA" << '\t' << "NA" << '\t';
		}
		tsvOut << stats.barcodes1 << '\t' << stats.barcodes2 << '\t' << stats.barcodesUnion << '\t'
		       << stats.barcodesIntersect << '\n';

		tsvOut << pair.second << (sense2 ? '+' : '-') << '\t' << pair.first << (sense1 ? '+' : '-')
		       << '\t';
		if (g[e].jaccard >= 0) {
			tsvOut << g[e].minDist << '\t' << g[e].dist << '\t' << g[e].maxDist << '\t';
		} else {
			tsvOut << "NA" << '\t' << "NA" << '\t' << "NA" << '\t';
		}
		tsvOut << stats.barcodes2 << '\t' << stats.barcodes1 << '\t' << stats.barcodesUnion << '\t'
		       << stats.barcodesIntersect << '\n';

		assert(tsvOut);
	}

	tsvOut.close();
}

/**
 * Write distance samples to an output stream.  The distance
 * samples record the distance between the head and tail regions
 * of the same contig with associated barcode stats (e.g.
 * barcode intersection size).
 */
static inline std::ostream&
writeDistSamplesTSV(std::ostream& out, const DistSampleMap& distSamples)
{
	out << "contig_id" << '\t' << "distance" << '\t' << "barcodes_head" << '\t' << "barcodes_tail"
	    << '\t' << "barcodes_union" << '\t' << "barcodes_intersect" << '\n';

	for (DistSampleConstIt it = distSamples.begin(); it != distSamples.end(); ++it) {
		const std::string& contigID = it->first;
		const DistSample& sample = it->second;

		out << contigID << '\t' << sample.distance << '\t' << sample.barcodesHead << '\t'
		    << sample.barcodesTail << '\t' << sample.barcodesUnion << '\t'
		    << sample.barcodesIntersect << '\n';
	}

	return out;
}

/**
 * Write intra-contig distance samples and barcode counts to a
 * TSV file.
 */
static inline void
writeDistSamplesTSV(const std::string& path, const DistSampleMap& distSamples)
{
	if (path.empty())
		return;

	ofstream samplesOut;
	samplesOut.open(path.c_str());
	assert(samplesOut);
	writeDistSamplesTSV(samplesOut, distSamples);
	assert(samplesOut);
	samplesOut.close();
}

#endif
