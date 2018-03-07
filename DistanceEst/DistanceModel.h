#ifndef _DISTANCE_MODEL_H_
#define _DISTANCE_MODEL_H_ 1

#include "Barcode.h"
#include "ScaffSize.h"
#include "Segment.h"
#include "SetUtil.h"

#include <cassert>

class DistanceModel
{
public:

	/** stores a histogram of inter-segment barcode Jaccard scores */
	typedef std::vector<unsigned> Histogram;

	/**
	* stores Jaccard score histograms for integer multiples
	* of the segment length
	*/
	typedef std::vector<Histogram> Model;

	/** a contig name and length */
	typedef std::pair<std::string, int> ContigSizeRec;
	typedef std::vector<ContigSizeRec> ContigSizeList;
	typedef typename ContigSizeList::const_iterator ContigSizeIt;

	/** constructor */
	DistanceModel(unsigned segmentSize, unsigned histBins)
		: m_segmentSize(segmentSize), m_histBins(histBins) {}

	void loadBarcodeMappings(
		const ScaffSizeList& scaffSizes,
		const SegmentToBarcode& segmentToBarcode)
	{
		float binWidth = 1.0f / m_histBins;

		SegmentCalc calc(m_segmentSize);

		for (ScaffSizeConstIt it = scaffSizes.begin();
			 it != scaffSizes.end(); ++it)
		{
			const ScaffSizeRec& rec = *it;
			const std::string& id = rec.first;
			unsigned length = rec.second;

			SegmentPairIterator pairIt(id, length, m_segmentSize);
			SegmentPairIterator pairEnd;
			for (; pairIt != pairEnd; ++pairIt) {

				/* extract barcode sets for segments */

				BarcodeList barcodes1, barcodes2;

				barcodes1.reserve(1000);
				barcodes2.reserve(1000);

				const Segment& segment1 = pairIt->first;
				const Segment& segment2 = pairIt->second;

				addBarcodes(segment1, segmentToBarcode, barcodes1);
				addBarcodes(segment2, segmentToBarcode, barcodes2);

				if (barcodes1.size() == 0 || barcodes2.size() == 0)
					continue;

				/* calculate distance and Jaccard score */

				unsigned dist = calc.start(length, segment2.second)
					- calc.start(length, segment1.second);

				float jaccard = (float)intersectionSize(barcodes1, barcodes2)
					/ unionSize(barcodes1, barcodes2);

				/* increment count in histogram bin */

				assert(dist % m_segmentSize == 0);
				unsigned histIndex = dist / m_segmentSize - 1;
				if (histIndex >= m_model.size())
					m_model.resize(histIndex + 1);

				assert(jaccard >= 0.0f && jaccard <= 1.0f);
				unsigned binIndex = std::floor(jaccard / binWidth);

				/* special case: jaccard scores of 1.0f go in last bin */
				binIndex = std::min(binIndex, m_histBins - 1);

				Histogram& hist = m_model.at(histIndex);
				if (binIndex >= hist.size())
					hist.resize(binIndex + 1);

				hist.at(binIndex)++;

			}

		}
	}

	friend std::ostream& operator <<(std::ostream& out,
		const DistanceModel& o)
	{
		out << std::fixed << std::setprecision(2);

		out << "dist" << '\t'
			<< "start" << '\t'
			<< "end" << '\t'
			<< "count" << '\n';

		float binWidth = 1.0f / o.m_histBins;
		unsigned dist = 0;
		for (const auto& hist : o.m_model) {
			dist += o.m_segmentSize;
			float start = 0.0f;
			for (const auto& count : hist) {
				float end = start + binWidth;
				out << dist << '\t'
					<< start << '\t'
					<< end << '\t'
					<< count << '\n';
				start += binWidth;
			}
		}

		return out;
	}

protected:

	/** contig segment size */
	unsigned m_segmentSize;
	/** number of bins for Jaccard score histograms */
	unsigned m_histBins;

	/**
	 * Jaccard score histograms for distances at exact multiples
	 * of the segment length
	 */
	Model m_model;

};

#endif
