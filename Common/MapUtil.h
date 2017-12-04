#ifndef _MAP_UTIL_H_
#define _MAP_UTIL_H_ 1

#include <map>
#include <cmath>
#include <utility>

/** Return iterator to the element with the closest key */
template <typename MapT>
typename MapT::const_iterator
closestKey(const MapT& map, double key)
{
	typedef typename MapT::const_iterator It;

	if (map.size() == 0)
		return map.end();

	// find first key greater than `key`
	It it = map.lower_bound(key);

	// first element is greater than `key`
	if (it == map.begin())
		return it;

	// all elements are less than `key`, so return last
	if (it == map.end()) {
		--it;
		return it;
	}

	// compare first element greater than `key` with
	// preceding element

    It prevIt = it;
	--prevIt;

	double diff1 = fabs(key - prevIt->first);
	double diff2 = fabs(key - it->first);

	if (diff1 > diff2)
		return it;
	else
		return prevIt;
}

/** Return iterator range for elements with closest n keys */
template <typename MapT>
std::pair<typename MapT::const_iterator, typename MapT::const_iterator>
closestKeys(const MapT& map, double key, size_t n)
{
	typedef typename MapT::const_iterator It;

	if (map.size() == 0)
		return std::make_pair(map.end(), map.end());

	// initialize range to closest key

	It first, last;
	first = closestKey(map, key);
	last = first;
	++last;

	// iteratively expand range left/right by next closest key

	size_t count = 1;
	while (count < n)
	{
		if (first == map.begin() && last == map.end())
			break;

		if (first == map.begin()) {
			++last;
		} else if (last == map.end()) {
			--first;
		} else {
			assert(first != map.begin());
			assert(last != map.end());

			It prevIt = first;
			--prevIt;

			double diff1 = fabs(key - prevIt->first);
			double diff2 = fabs(key - last->first);
			if (diff1 < diff2)
				first = prevIt;
			else
				++last;
		}

		++count;
	}

	return std::make_pair(first, last);
}

#endif
