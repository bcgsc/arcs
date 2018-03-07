#define CATCH_CONFIG_MAIN
#include "ThirdParty/Catch/catch.hpp"

#include "DataLayer/SAMUtil.h"

using namespace std;

TEST_CASE("getBarcodeSeq", "[SAMUtil]")
{
    const string tags1("QT:Z:AA<FFKKK BX:Z:CGTCAGGTCAGAGGTG-1 XT:i:0");
    const string tags2("QT:Z:AA<FFKKK XT:i:0");

    REQUIRE(getBarcodeSeq(tags1) == "CGTCAGGTCAGAGGTG");
    REQUIRE(getBarcodeSeq(tags2).empty());
}
