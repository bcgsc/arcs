#define CATCH_CONFIG_MAIN
#include "ThirdParty/Catch/catch.hpp"

#include "DataStructures/Segment.h"

#include <limits>

TEST_CASE("SegmentCalc.index() (no remainder)", "[Segment]")
{
    const unsigned seqSize = 9;
    const unsigned segmentSize = 3;

    SegmentCalc calc(segmentSize);

    /* position lands in left half of seq */
    REQUIRE(calc.index(3, seqSize) == 0);

    /* position lands in middle remainder seg */
    REQUIRE(calc.index(4, seqSize) == 1);

    /* position lands in right half of seq */
    REQUIRE(calc.index(7, seqSize) == 2);
}

TEST_CASE("SegmentCalc.index() (with remainder)", "[Segment]")
{
    const unsigned NO_INDEX = std::numeric_limits<unsigned>::max();

    const unsigned seqSize = 11;
    const unsigned segmentSize = 3;

    SegmentCalc calc(segmentSize);

    /* position lands in left half of seq */
    REQUIRE(calc.index(2, seqSize) == 0);

    /* position lands in middle remainder seg */
    REQUIRE(calc.index(4, seqSize) == NO_INDEX);

    /* position lands in right half of seq */
    REQUIRE(calc.index(10, seqSize) == 1);
}

TEST_CASE("SegmentPairIterator (no remainder)", "[Segment]")
{
    const unsigned NO_INDEX = std::numeric_limits<unsigned>::max();

	const std::string id("id");
    const unsigned seqSize = 9;
    const unsigned segmentSize = 3;

    SegmentPairIterator it(id, seqSize, segmentSize);
	SegmentPairIterator end;

    REQUIRE(it != end);
	REQUIRE(it->first.second == 0);
	REQUIRE(it->second.second == 1);
	++it;
	REQUIRE(it != end);
	REQUIRE(it->first.second == 0);
	REQUIRE(it->second.second == 2);
	++it;
	REQUIRE(it != end);
	REQUIRE(it->first.second == 1);
	REQUIRE(it->second.second == 2);
	++it;
	REQUIRE(it == end);
}

TEST_CASE("SegmentPairIterator (with remainder)", "[Segment]")
{
    const unsigned NO_INDEX = std::numeric_limits<unsigned>::max();

	const std::string id("id");
    const unsigned seqSize = 13;
    const unsigned segmentSize = 2;

    SegmentPairIterator it(id, seqSize, segmentSize);
	SegmentPairIterator end;

    REQUIRE(it != end);
	REQUIRE(it->first.second == 0);
	REQUIRE(it->second.second == 1);
	++it;
	REQUIRE(it != end);
	REQUIRE(it->first.second == 0);
	REQUIRE(it->second.second == 2);
	++it;
	REQUIRE(it != end);
	REQUIRE(it->first.second == 1);
	REQUIRE(it->second.second == 2);
	++it;
	REQUIRE(it != end);
	REQUIRE(it->first.second == 3);
	REQUIRE(it->second.second == 4);
	++it;
	REQUIRE(it != end);
	REQUIRE(it->first.second == 3);
	REQUIRE(it->second.second == 5);
	++it;
	REQUIRE(it != end);
	REQUIRE(it->first.second == 4);
	REQUIRE(it->second.second == 5);
	++it;
	REQUIRE(it == end);
}
