#define CATCH_CONFIG_MAIN GenepopTest
#include <catch2/catch.hpp>

#include "input.hpp"

TEST_CASE("haplo_gen_pop_input_test")
{
    SECTION("read_test")
    {
        genepop_input_c<1> input("genotype_genepop_format.txt");
        //Verify ADH-4 , ADH-5 \n mtDNA have been well separate
        REQUIRE(input.Locus_name.size() == 11);
        REQUIRE(input.Locus_name[5] == "4");
        REQUIRE(input.Locus_name[6] == "5");
        //pop name
        REQUIRE(input.Pop_name.size() == 5);
        REQUIRE(input.Indiv_name.size() == 5);
        REQUIRE(input.Pop_name[0] == "00");
        REQUIRE(input.Pop_name[3] == "03");
        //Name of indiv
        REQUIRE(input.Indiv_name.size() == 5);
        REQUIRE(input.Indiv_name[0][0] == "00");
        REQUIRE(input.Indiv_name[2][1] == "02");
        REQUIRE(input.Indiv_name[3][2] == "03");
        //Number of pop
        REQUIRE(input.Genotype.size() == 5);
        //Pop size
        REQUIRE(input.Genotype[0].size() == 10);
        //First indiv
        REQUIRE(input.Genotype[0][0].at(0) == std::array{37});
        REQUIRE(input.Genotype[0][0].at(5) == std::array{39});
        //Last indiv -1
        REQUIRE(input.Genotype[3][2].at(0) == std::array{1});
        REQUIRE(input.Genotype[3][2].at(4) == std::array{51});
    }
}