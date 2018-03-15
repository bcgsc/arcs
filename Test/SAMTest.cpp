#define CATCH_CONFIG_MAIN
#include "ThirdParty/Catch/catch.hpp"

#include "Common/SAM.h"

using namespace std;

TEST_CASE("parseBXTag", "[SAM]")
{
    const string tags1("QT:Z:AA<FFKKK BX:Z:CGTCAGGTCAGAGGTG-1 XT:i:0");
    const string tags2("QT:Z:AA<FFKKK XT:i:0");

    REQUIRE(parseBXTag(tags1) == "CGTCAGGTCAGAGGTG-1");
    REQUIRE(parseBXTag(tags2).empty());
}
