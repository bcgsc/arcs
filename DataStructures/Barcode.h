#ifndef _BARCODE_H_
#define _BARCODE_H_ 1

#include <boost/tuple/tuple.hpp>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/** Chromium barcode sequence */
typedef std::string BarcodeSeq;

/** Unique index assigned to each Chromium barcode sequence */
typedef uint32_t BarcodeIndex;

/** Number of barcode mappings per segment */
typedef uint16_t BarcodeCount;

/** A set of Chromium barcode indices */
typedef std::unordered_set<BarcodeIndex> BarcodeSet;
typedef typename BarcodeSet::const_iterator BarcodeSetConstIt;

/** A list of Chromium barcodes */
typedef std::vector<BarcodeIndex> BarcodeList;
typedef typename BarcodeList::iterator BarcodeListIt;

/** Maps barcode => number of reads */
typedef std::vector<BarcodeIndex> BarcodeMultMap;

/** Assign/query unique integer indices for Chromium barcode sequences */
class BarcodeToIndexMap
{
public:

    BarcodeToIndexMap() : m_nextBarcodeIndex(0) {}

    BarcodeIndex getIndex(const BarcodeSeq& barcode)
    {
        BarcodeToIndexIt barcodeIt = m_barcodeToIndex.find(barcode);
        if (barcodeIt == m_barcodeToIndex.end()) {
            bool inserted;
            boost::tie(barcodeIt, inserted) =
                m_barcodeToIndex.insert(
                std::make_pair(barcode, m_nextBarcodeIndex++));
            assert(inserted);
        }
        return barcodeIt->second;
    }

protected:

    /** Maps Chromium barcode sequence => unique index */
    typedef std::unordered_map<BarcodeSeq, BarcodeIndex> BarcodeToIndex;
    typedef typename BarcodeToIndex::iterator BarcodeToIndexIt;

    BarcodeToIndex m_barcodeToIndex;
    unsigned m_nextBarcodeIndex;
};


#endif
