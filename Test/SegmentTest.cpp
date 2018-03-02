#define CATCH_CONFIG_MAIN
#include "ThirdParty/Catch/catch.hpp"

#include "Common/Segment.h"

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

