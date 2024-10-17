#include <iostream>
#include <array>
#include <cstdio>
#include <unordered_map>
#include <filesystem>
#include <string>

#include "test/Util.h"
#include "src/extract/CNFBaseFeatures.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

namespace fs = std::filesystem;

// TEST_CASE("GBDLib")
// {
//     SUBCASE("BaseFeature_Names")
//     {
//         auto ex = CNF::BaseFeatures("");
//         auto super = feature_names<CNF::BaseFeatures<>>();
//         auto sub = ex.getNames();
//         check_subset(sub, super);
//         CHECK_EQ(sub.size(), super.size() - 1);
//         CHECK_EQ(super.back(), ex.getRuntimeDesc());
//     }

