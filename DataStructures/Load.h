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
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iostream>
#include <cstdlib>

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
	size_t nmTag;
	size_t belowSeqIdentity;

	LoadCounters() : processed(0), unmapped(0), noBarcode(0),
		belowSeqIdentity(0) {}

	friend std::ostream& operator <<(std::ostream& out,
		LoadCounters& o)
	{
		out << "processed " << o.processed
			<< " alignments ("
			<< "unmapped: " << o.unmapped << ", "
			<< "no barcode: " << o.noBarcode << ", "
			<< "below seq identity: " << o.belowSeqIdentity << ")\n";
		return out;
	}
};

/** calculate percent sequence identity of a SAM alignment record */
static inline std::pair<double, bool> identity(const bam1_t* rec)
{
	int qalen = bam_cigar2qlen(rec->core.n_cigar, bam_get_cigar(rec));
	assert(qalen > 0);

	uint8_t* tagData = bam_aux_get(rec, "NM");
	if (tagData == NULL)
		return std::make_pair(100.0, false);

	int64_t editDist = bam_aux2i(tagData);
	if (errno == EINVAL)
		return std::make_pair(100.0, false);

	assert(editDist >= 0);
	assert(rec->core.l_qseq > 0);
	double _identity = double(qalen - editDist) / rec->core.l_qseq * 100.0;

	return std::make_pair(_identity, true);
}

/** load SAM data into ARCS data structures */
static inline void readSAM(const std::string& path,
	const LoadParams& params,
	BarcodeToIndexMap& barcodeToIndex,
	BarcodeMultMap& barcodeMult,
	SegmentToBarcode& segmentToBarcode)
{
	/* open SAM file */
	samFile *in = sam_open(path.c_str(), "r");
	if (in == NULL) {
		std::cerr << "error: failed to open SAM file `" << path << "`\n";
		exit(EXIT_FAILURE);
	}

	/* read SAM headers */
	bam_hdr_t* header = sam_hdr_read(in);

	/* read SAM records */
	LoadCounters counters;
	bam1_t* rec = bam_init1();
	int result;
	while ((result = sam_read1(in, header, rec)) >= 0) {

		/* progress message */
		if (params.verbose && counters.processed > 0
			&& counters.processed % params.progressStep == 0) {
			std::cerr << counters;
		}
		counters.processed++;

		/* skip unmapped reads */
		if (rec->core.flag & BAM_FUNMAP) {
			counters.unmapped++;
			continue;
		}

		/* get Chromium barcode */
		uint8_t* barcodePtr = bam_aux_get(rec, "BX");
		if (barcodePtr == NULL) {
			counters.noBarcode++;
			continue;
		}

		/* get integer index for barcode sequence */
		std::string barcode((char*)barcodePtr);
		BarcodeIndex barcodeIndex = barcodeToIndex.getIndex(barcode);

		/* increment barcode multiplicity (read count) */
		if (barcodeIndex >= barcodeMult.size())
			barcodeMult.resize(barcodeIndex + 1);
		barcodeMult.at(barcodeIndex)++;

		/* check for sufficient sequence identity (if 'NM:i' tag is available) */
		double _identity;
		bool nmTag;
		boost::tie(_identity, nmTag) = identity(rec);
		if (nmTag) {
			counters.nmTag++;
			if (_identity < params.identity) {
				counters.belowSeqIdentity++;
				continue;
			}
		}

		/* identify contig segments overlapping alignment */

		SegmentCalc calc(params.segmentSize);

		char* rname = header->target_name[rec->core.tid];
		uint32_t rlen = header->target_len[rec->core.tid];
		assert(rlen > 0);
		int32_t pos = rec->core.pos + 1;
		assert(pos > 0);
		int ralen = bam_cigar2rlen(rec->core.n_cigar, bam_get_cigar(rec));
		assert(ralen > 0);

		std::pair<SegmentIndex, SegmentIndex> range;
		bool valid;
		boost::tie(range, valid) = calc.indexRange(pos, pos + ralen - 1, rlen);

		/* discard alignment to middle remainder segment */
		if (!valid)
			continue;

		/* discard alignments that span multiple segments */
		if (range.first != range.second)
			continue;

		/* record segment-to-barcode mappings */
		Segment segment(rname, range.first);
		segmentToBarcode[segment][barcodeIndex]++;

	}

	/* check if we successfully read to the end of the SAM stream/file */
	if (result < -1) {
		std::cerr << "error: truncated SAM file\n";
		exit(EXIT_FAILURE);
	}

	/* print final progress message */
	if (params.verbose)
		std::cerr << counters;
}

#endif
