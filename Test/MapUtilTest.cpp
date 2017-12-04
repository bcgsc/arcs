#define CATCH_CONFIG_MAIN
#include "ThirdParty/Catch/catch.hpp"

#include "Common/MapUtil.h"
#include <boost/tuple/tuple.hpp>

using namespace std;

TEST_CASE("closestKey", "[MapUtil]")
{
	map<int, int> data = {{1, 1}, {2, 2}, {3, 3}, {4, 4}};
	map<int, int> empty;
	map<int, int>::const_iterator it;

	// empty map => return end()

	it = closestKey(empty, 1);
	REQUIRE(it == empty.end());

	// exact match

	it = closestKey(data, 3);
	REQUIRE(it->first == 3);

	// tie => return element that occurs first in the map

	it = closestKey(data, 2.5);
	REQUIRE(it->first == 2);

	// smaller than all keys => return first element

	it = closestKey(data, 0);
	REQUIRE(it == data.begin());

	// greater than all keys => return last element

	it = closestKey(data, 5);
	REQUIRE(it->first == 4);
}

TEST_CASE("closestKeys", "[MapUtil]")
{
	map<int, int> data = {{1, 1}, {2, 2}, {3, 3}, {4, 4}};
	map<int, int> empty;
	map<int, int>::const_iterator first, last;

	// empty map => return (end(), end())

	boost::tie(first, last) = closestKeys(empty, 3, 3);
	REQUIRE(first == empty.end());
	REQUIRE(last == empty.end());

	// exact match, n = 1

	boost::tie(first, last) = closestKeys(data, 2, 1);
	REQUIRE(std::distance(first, last) == 1);
	REQUIRE(first->first == 2);

	// inexact match, n = 2

	boost::tie(first, last) = closestKeys(data, 2.5, 2);
	REQUIRE(std::distance(first, last) == 2);
	REQUIRE(first->first == 2);
	++first;
	REQUIRE(first->first == 3);

	// less than all keys => first n elements

	boost::tie(first, last) = closestKeys(data, 0, 2);
	REQUIRE(std::distance(first, last) == 2);
	REQUIRE(first->first == 1);
	++first;
	REQUIRE(first->first == 2);

	// greater than all keys => last n elements

	boost::tie(first, last) = closestKeys(data, 5, 2);
	REQUIRE(std::distance(first, last) == 2);
	REQUIRE(first->first == 3);
	++first;
	REQUIRE(first->first == 4);
}
