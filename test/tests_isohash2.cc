#include <string>
#include <vector>
#include <random>
#include <filesystem>
#include <algorithm>
#include <exception>
#include <iostream>

#include "src/identify/ISOHash2.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

namespace fs = std::filesystem;

static fs::path find_scrambled_root() {
    const fs::path candidates[] = {
        fs::path("test/resources/scrambled/indepth"),
        fs::path("../test/resources/scrambled/indepth"),
        fs::path("resources/scrambled/indepth"),
        fs::path("../resources/scrambled/indepth"),
        fs::path("../../test/resources/scrambled/indepth"),
        fs::path("../../resources/scrambled/indepth"),
    };
    for (const auto& p : candidates) {
        if (fs::exists(p) && fs::is_directory(p)) return p;
    }
    return {};
}

static std::vector<fs::path> list_sorted_files(const fs::path& dir) {
    std::vector<fs::path> files;
    for (const auto& e : fs::directory_iterator(dir)) {
        if (e.is_regular_file()) files.push_back(e.path());
    }
    std::sort(files.begin(), files.end(),
              [](const fs::path& a, const fs::path& b) { return a.string() < b.string(); });
    return files;
}

static constexpr std::size_t SAMPLE_SIZE = 10;

static std::vector<fs::path> sample_files(std::vector<fs::path> files, std::mt19937& rng) {
    files.erase(std::remove_if(files.begin(), files.end(),
        [](const fs::path& p) { const auto s = p.filename().string(); return !s.empty() && s[0] == '.'; }),
        files.end());
    if (files.size() > SAMPLE_SIZE) {
        std::shuffle(files.begin(), files.end(), rng);
        files.resize(SAMPLE_SIZE);
    }
    return files;
}

TEST_CASE("IsoHash2 Robustness") {
    const fs::path scrambled_root = find_scrambled_root();
    REQUIRE_MESSAGE(!scrambled_root.empty(),
        "Cannot find scrambled test resources, tried several relative paths");

    CNF::IsoHash2Settings config;
    config.max_iterations = 6;

    bool saw_any_family = false;

    std::mt19937 rng(42);

    std::vector<fs::path> families;
    for (const auto& e : fs::directory_iterator(scrambled_root)) {
        if (e.is_directory()) families.push_back(e.path());
    }
    std::sort(families.begin(), families.end(),
              [](const fs::path& a, const fs::path& b) { return a.string() < b.string(); });

    for (const auto& fam_dir : families) {
        saw_any_family = true;

        const std::string instance_name = fam_dir.filename().string();
        const std::string subcase_name = "Instance: " + instance_name;

        SUBCASE(subcase_name.c_str()) {
            const auto files = sample_files(list_sorted_files(fam_dir), rng);
            REQUIRE_MESSAGE(!files.empty(), ("No files found in " + fam_dir.string()).c_str());

            std::string expected_hash, reference_file;
            std::size_t tested_ok = 0;

            for (const auto& path : files) {
                const std::string filepath = path.string();
                std::string current_hash;
                try {
                    current_hash = CNF::isohash2(filepath.c_str(), config);
                } catch (const std::exception& e) {
                    std::cerr << "[IsoHash2] EXCEPTION instance=" << instance_name
                              << " file=" << filepath
                              << " what=" << e.what() << std::endl;
                    FAIL_CHECK(("Exception during hashing: " + filepath + " : " + e.what()).c_str());
                    continue;
                }
                if (expected_hash.empty()) {
                    expected_hash = current_hash;
                    reference_file = filepath;
                } else if (current_hash != expected_hash) {
                    std::cerr << "[IsoHash2] MISMATCH instance=" << instance_name
                              << "\n  ref_file=" << reference_file
                              << "\n  ref_hash=" << expected_hash
                              << "\n  cur_file=" << filepath
                              << "\n  cur_hash=" << current_hash
                              << std::endl;
                    CHECK_MESSAGE(false,
                        ("\nHash mismatch!"
                         "\nReference: " + reference_file + " -> " + expected_hash +
                         "\nCurrent:   " + filepath + " -> " + current_hash + "\n").c_str());
                }
                ++tested_ok;
            }

            REQUIRE_MESSAGE(tested_ok > 0, ("No hashable files in " + fam_dir.string()).c_str());

            std::cerr << "[IsoHash2] SUMMARY instance=" << instance_name
                      << " files=" << tested_ok
                      << " hash=" << (expected_hash.empty() ? std::string("<none>") : expected_hash)
                      << std::endl;
        }
    }

    REQUIRE_MESSAGE(saw_any_family, ("No family directories under " + scrambled_root.string()).c_str());
}
