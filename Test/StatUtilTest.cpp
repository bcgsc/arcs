#define CATCH_CONFIG_MAIN
#include "ThirdParty/Catch/catch.hpp"

#include "Common/StatUtil.h"
#include <ctgmath>
#include <iostream>

using namespace std;

template <class T>
bool approxEqual(const T& a, const T& b, const T& epsilon) {
	return fabs(a - b) < epsilon;
}

TEST_CASE("quantile calculation", "[StatUtil]")
{
	// allowable error for floating point comparison

	const double epsilon = 0.001;

	// array with even number of elements

	array<int, 4> data1 = { 1, 2, 3, 4 };
	REQUIRE(approxEqual(2.5, quantile(data1.begin(), data1.end(), 0.5), epsilon));

	// array with odd number of elements

	array<int, 5> data2 = { 1, 2, 3, 4, 5 };
	REQUIRE(approxEqual(3.0, quantile(data2.begin(), data2.end(), 0.5), epsilon));

}
