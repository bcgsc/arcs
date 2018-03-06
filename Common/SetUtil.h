#ifndef _SET_UTIL_H_
#define _SET_UTIL_H_ 1

template <class List>
static inline unsigned unionSize(List& list1, List& list2)
{
	List _union(list1.size() + list2.size());
	typename List::iterator end;

	std::sort(list1.begin(), list1.end());
	std::sort(list2.begin(), list2.end());

	end = std::set_union(list1.begin(), list1.end(),
		list2.begin(), list2.end(), _union.begin());

	return end - _union.begin();
}

template <class List>
static inline unsigned intersectionSize(List& list1, List& list2)
{
	List intersection(list1.size() + list2.size());
	typename List::iterator end;

	std::sort(list1.begin(), list1.end());
	std::sort(list2.begin(), list2.end());

	end = std::set_intersection(list1.begin(), list1.end(),
		list2.begin(), list2.end(), intersection.begin());

	return end - intersection.begin();
}

#endif
