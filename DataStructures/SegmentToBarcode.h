#ifndef _BARCODE_INDEX_H_
#define _BARCODE_INDEX_H_ 1

#include "Barcode.h"
#include "Segment.h"

#include <map>
#include <unordered_map>

/** Maps barcode => number of read mappings */
typedef std::map<BarcodeIndex, unsigned> BarcodeToCount;
typedef typename BarcodeToCount::const_iterator BarcodeToCountConstIt;

/** Maps contig segment => barcodes / mapping counts */
typedef std::unordered_map<Segment, BarcodeToCount, HashSegment>
	SegmentToBarcode;
typedef typename SegmentToBarcode::const_iterator SegmentToBarcodeConstIt;
typedef typename SegmentToBarcode::iterator SegmentToBarcodeIt;

/** get barcodes for the given segment and append to `out` */
static inline void addBarcodes(const Segment& segment,
	const SegmentToBarcode& segmentToBarcode,
	std::vector<BarcodeIndex>& out)
{
	SegmentToBarcodeConstIt segmentIt = segmentToBarcode.find(segment);
	if (segmentIt == segmentToBarcode.end())
		return;
	for (const auto& rec : segmentIt->second)
		out.push_back(rec.first);
}

#endif
