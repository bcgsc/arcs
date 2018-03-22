#ifndef _LOAD_H_
#define _LOAD_H_ 1

#include "DataLayer/SAM.h"
#include "DataLayer/SAMUtil.h"
#include "DataStructures/Barcode.h"
#include "DataStructures/ScaffSize.h"
#include "DataStructures/Segment.h"
#include "DataStructures/SegmentToBarcode.h"

#include <boost/tuple/tuple.hpp>
#include <fstream>
#include <iostream>

struct LoadParams
{
	/**
	 * If true, extract Chromium barcode from BX tag.
	 * If false, extract Chromium barcode from string
	 * following last underscore ("_") in read name
	 * (legacy ARCS behaviour).
	 */
	bool bx;
	/** min sequence identity for read alignments */
	float identity;
	/** number of alignments between progress messages */
	size_t progressStep;
	/** contig segment length */
	unsigned segmentSize;
	/** verbose logging */
	bool verbose;

	LoadParams()
		: bx(true), identity(98.0f), progressStep(100000),
		segmentSize(1000), verbose(false) {}
};

struct LoadCounters
{
	size_t processed;
	size_t unmapped;
	size_t noBarcode;
	size_t belowSeqIdentity;

	LoadCounters() : processed(0), unmapped(0), noBarcode(0),
		belowSeqIdentity(0) {}

	friend std::ostream& operator <<(std::ostream& out,
		LoadCounters& o)
	{
		out << "Processed " << o.processed
			<< " alignments ("
			<< "unmapped: " << o.unmapped << ", "
			<< "no barcode: " << o.noBarcode << ", "
			<< "below seq identity: " << o.belowSeqIdentity << ")\n";
		return out;
	}
};

static inline std::istream& readSAM(std::istream& in,
	const ScaffSizeMap& scaffSizes,
	const LoadParams& params,
	BarcodeToIndexMap& barcodeToIndex,
	BarcodeMultMap& barcodeMult,
	SegmentToBarcode& segmentToBarcode)
{
	if (params.verbose)
		std::cerr << "Loading alignments..." << std::endl;

	/* skip SAM header lines */

	while (in.peek() == '@') {
		std::string line;
		getline(in, line);
		assert(in);
	}

	/* read SAM alignment lines */

	LoadCounters counters;

	SAMRecord sam;
	for (; in >> sam; counters.processed++) {

		/* progress message */

		if (params.verbose && counters.processed > 0
			&& counters.processed % params.progressStep == 0) {
			std::cerr << counters;
		}

		/* skip unmapped reads */

		if (sam.isUnmapped()) {
			counters.unmapped++;
			continue;
		}

		/* get Chromium barcode for read */

		std::string barcode;
		if (params.bx) {
			barcode = getBarcodeSeq(sam.tags);
		} else {
			std::size_t found = sam.qname.rfind("_");
			if (found!=std::string::npos)
				barcode = sam.qname.substr(found+1);
		}
		if (barcode.empty() || !checkIndex(barcode)) {
			counters.noBarcode++;
			continue;
		}

		/* get integer Index for barcode sequence */

		BarcodeIndex barcodeIndex = barcodeToIndex.getIndex(barcode);

		/* increment barcode multiplicity (read count) */

		if (barcodeIndex >= barcodeMult.size())
			barcodeMult.resize(barcodeIndex + 1);
		barcodeMult.at(barcodeIndex)++;

		/* check for sufficient sequence identity */

		int si = calcSequenceIdentity(sam.tags, sam.cigar, sam.seq);
		if ((unsigned)si < params.identity)
			continue;

		/* calculate contig segments overlapping alignment */

		assert(sam.cigar != "*");
		SAMAlignment::CigarCoord cigar(sam.cigar);
		SegmentCalc calc(params.segmentSize);
		unsigned length = scaffSizes.at(sam.rname);
		std::pair<SegmentIndex, SegmentIndex> range;
		bool valid;

		/*
		 * `SAMRecord` provides 0-based pos but `SegmentCalc`
		 * requires 1-based pos
		 */
		sam.pos++;

		boost::tie(range, valid) = calc.indexRange(
			sam.pos, sam.pos + cigar.tspan - 1, length);
		if (!valid)
			continue;

		/* discard alignments that span multiple segments */
		if (range.first != range.second)
			continue;

		Segment segment(sam.rname, range.first);
		segmentToBarcode[segment][barcodeIndex]++;

	}
	assert(in.eof());

	if (params.verbose)
		std::cerr << counters;

	return in;
}

static inline void readSAM(const std::string& path,
	const ScaffSizeMap& scaffSizes,
	const LoadParams& params,
	BarcodeToIndexMap& barcodeToIndex,
	BarcodeMultMap& barcodeMult,
	SegmentToBarcode& segmentToBarcode)
{
	if (path == "-") {
		readSAM(std::cin, scaffSizes, params, barcodeToIndex,
			barcodeMult, segmentToBarcode);
	} else {
		std::ifstream in(path);
		assert_good(in, path);
		readSAM(in, scaffSizes, params, barcodeToIndex,
			barcodeMult, segmentToBarcode);
		in.close();
	}
}

#endif
