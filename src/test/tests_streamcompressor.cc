#include <stdio.h>
#include <filesystem>
#include <iostream>
#include <fstream>
#include "src/util/StreamBuffer.h"
#include "src/util/StreamCompressor.h"
#include "src/test/Util.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

TEST_CASE("StreamCompressor")
{
    SUBCASE("Write to archive")
    {
        std::cerr << std::filesystem::current_path() << "\n";
        auto tmp_file = tmp_filename("src/test/resources", ".cnf.xz");
        StreamCompressor c(tmp_file.c_str(), 30);
        const char *data = "p cnf 1 2\n1 2 0\n1 0\n-2 3 0";
        c.write(data, strlen(data));
        c.close();
        Cl cl1 = {Lit(1, false), Lit(2, false)};
        Cl cl2 = {Lit(1, false)};
        Cl cl3 = {Lit(2, true), Lit(3, false)};
        std::vector<Cl> clauses{cl1, cl2, cl3};
        StreamBuffer b(tmp_file.c_str());
        Cl clause;
        for (int i = 0; b.readClause(clause); ++i)
        {
            for (Lit l : clause)
            {
                CHECK(clause == clauses[i]);
            }
        }
        remove(tmp_file.c_str());
    }

    SUBCASE("Write from istream to archive")
    {
        auto tmp_file = tmp_filename("src/test/resources", ".cnf.xz");
        const char *cnf_file = "src/test/resources/test_files/ibm-2004-03-k70.cnf";
        StreamBuffer cnf_buf(cnf_file);

        StreamCompressor cmpr(tmp_file.c_str());
        std::ifstream in(cnf_file);
        in >> cmpr;
        cmpr.close();
        in.close();
        StreamBuffer tmp_buf(tmp_file.c_str());

        Cl tmp_cl, cnf_cl;
        bool cnf_cl_read = true, tmp_cl_read = true;
        while (cnf_cl_read && tmp_cl_read)
        {
            cnf_cl_read = cnf_buf.readClause(cnf_cl);
            tmp_cl_read = tmp_buf.readClause(tmp_cl);
            CHECK(cnf_cl == tmp_cl);
        }
        CHECK(cnf_cl_read == tmp_cl_read);
        remove(tmp_file.c_str());
    }
}

// int main()
// {
// }