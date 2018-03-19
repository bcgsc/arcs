#ifndef _BARCODE_TO_SEGMENT_H_
#define _BARCODE_TO_SEGMENT_H_ 1

#include "DataStructures/Barcode.h"
#include "DataStructures/Segment.h"
#include "DataStructures/SegmentToBarcode.h"

/** maps barcode index => contig segments */
typedef std::vector<SegmentList> BarcodeToSegment;
typedef typename BarcodeToSegment::const_iterator BarcodeToSegmentConstIt;

static inline void buildBarcodeToSegment(
	const SegmentToBarcode& segmentToBarcode,
	BarcodeToSegment& barcodeToSegment)
{
	for (SegmentToBarcodeConstIt segmentIt = segmentToBarcode.begin();
		segmentIt != segmentToBarcode.end(); ++segmentIt)
	{
		const Segment& segment = segmentIt->first;
		for (BarcodeToCountConstIt barcodeIt = segmentIt->second.begin();
			barcodeIt != segmentIt->second.end(); ++barcodeIt)
		{
			BarcodeIndex barcode = barcodeIt->first;
			if (barcode >= barcodeToSegment.size())
				barcodeToSegment.resize(barcode + 1);
			barcodeToSegment.at(barcode).push_back(segment);
		}
	}
}

#endif
