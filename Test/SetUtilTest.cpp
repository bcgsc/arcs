#define CATCH_CONFIG_MAIN
#include "ThirdParty/Catch/catch.hpp"

#include "Common/SetUtil.h"
#include <array>

TEST_CASE("unionSize / intersectionSize", "[SetUtil]")
{
	std::vector<int> forward = {1, 2, 3, 4, 5};
	std::vector<int> reverse = {5, 4, 3, 2, 1};

	unsigned _union = unionSize(forward, reverse);
	unsigned intersect = intersectionSize(forward, reverse);

	REQUIRE(_union == 5);
	REQUIRE(intersect == 5);
}
