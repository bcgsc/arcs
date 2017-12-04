#ifndef _STAT_UTIL_H_
#define _STAT_UTIL_H_ 1

#include <algorithm>
#include <cassert>
#include <cmath>

/** compute the qth quantile, where `q` is in the range [0,1] */
template <class IteratorT>
static inline double quantile(IteratorT it1, IteratorT it2, double q)
{
	assert(q >= 0.0);
	assert(q <= 1.0);

	size_t length = it2 - it1;
	assert(length >= 1);

	size_t lastPos = length - 1;

	/* get elements bordering quantile boundary */

	size_t beforePos = (size_t)floor(q * lastPos);
	size_t before = *(it1 + beforePos);

	size_t afterPos = (size_t)ceil(q * lastPos);
	size_t after = *(it1 + afterPos);

	/* weight elements by distance to quantile boundary */

	double weight = (q * lastPos - beforePos) / 1.0;
	return weight * before + (1.0 - weight) * after;
}

#endif
