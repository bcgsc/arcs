#define CATCH_CONFIG_MAIN
#include "ThirdParty/Catch/catch.hpp"

#include "Common/Segment.h"

TEST_CASE("posToSegmentIndex (no remainder)", "[SegmentUtil]")
{
    const unsigned seqLen = 9;
    const unsigned segLen = 3;

    /* position lands in left half of seq */
    REQUIRE(posToSegmentIndex(3, seqLen, segLen) == 0);

    /* position lands in middle remainder seg */
    REQUIRE(posToSegmentIndex(4, seqLen, segLen) == 1);

    /* position lands in right half of seq */
    REQUIRE(posToSegmentIndex(7, seqLen, segLen) == 2);
}

TEST_CASE("posToSegmentIndex (with remainder)", "[SegmentUtil]")
{
    const unsigned seqLen = 11;
    const unsigned segLen = 3;

    /* position lands in left half of seq */
    REQUIRE(posToSegmentIndex(2, seqLen, segLen) == 0);

    /* position lands in middle remainder seg */
    REQUIRE(posToSegmentIndex(4, seqLen, segLen)
        == std::numeric_limits<unsigned>::max());

    /* position lands in right half of seq */
    REQUIRE(posToSegmentIndex(10, seqLen, segLen) == 1);
}

