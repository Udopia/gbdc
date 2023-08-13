/**
 * Some tests for gbdc
 * 
 * @author Markus Iser 
 */

#include <stdio.h>
#include <filesystem>

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "src/util/StreamBuffer.h"

bool tempfile(FILE** file, char** name) {
    *name = tempnam("/tmp", "gbdc.test");
    *file = fopen(*name, "w+");
    return *file != nullptr;
}

TEST_CASE("StreamBuffer") {
    std::FILE* file;
    char* name;

    SUBCASE("read hello world") {
        CHECK(tempfile(&file, &name));
        std::fputs("Hello World!", file);
        std::fclose(file);
        StreamBuffer reader(name);
        CHECK(reader.skipString("Hello"));
        CHECK(reader.skipWhitespace());
        CHECK(!reader.skipString("World!"));
    }

    SUBCASE("read mixed: integers, strings, whitespace, linebreaks") {
        CHECK(tempfile(&file, &name));
        std::fputs("123 137 no   -7 mer\n ci\n", file);
        std::fclose(file);
        StreamBuffer reader(name);
        int num;
        CHECK(reader.readInteger(&num));
        CHECK(num == 123);
        CHECK(reader.skipWhitespace());
        CHECK(reader.readInteger(&num));
        CHECK(num == 137);
        CHECK(reader.skipWhitespace());
        CHECK(reader.skipString("no"));
        CHECK(reader.skipWhitespace());
        CHECK(reader.readInteger(&num));
        CHECK(num == -7);
        CHECK(reader.skipLine());
        CHECK(reader.skipString("ci"));
        CHECK(!reader.eof());
        CHECK(!reader.skipWhitespace());
        CHECK(reader.eof());
    }
}

// int main() {
//     return 0;
// }