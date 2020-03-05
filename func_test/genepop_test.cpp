#define CATCH_CONFIG_MAIN GenepopTest
#include <catch2/catch.hpp>
#include <iostream>

#include "input.hpp"
#include "common_tools.hpp"
#include "calc_stat.hpp"

TEST_CASE("haplo_gen_pop_input_test")
{
    SECTION("read_test")
    {
        genepop_input_c<1> input("genotype_genepop_format.txt");
        //Verify ADH-4 , ADH-5 \n mtDNA have been well separate
        REQUIRE(input.Locus_name.size() == 10);
        REQUIRE(input.Locus_name[4] == "4");
        REQUIRE(input.Locus_name[5] == "5");
        //pop name
        REQUIRE(input.Pop_name.size() == 5);
        REQUIRE(input.Indiv_name.size() == 5);
        REQUIRE(input.Pop_name[0][0] == "-0.1");
        REQUIRE(input.Pop_name[0][1] == "0.1");
        REQUIRE(input.Pop_name[3][0] == "-0.1");
        REQUIRE(input.Pop_name[3][1] == "3.1");
        //Name of indiv
        REQUIRE(input.Indiv_name.size() == 5);
        REQUIRE(input.Indiv_name[0][0] == "-0.1 0.1");
        REQUIRE(input.Indiv_name[2][1] == "-0.1 2.1");
        REQUIRE(input.Indiv_name[3][2] == "-0.1 3.1");
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

    SECTION("Dist_btw_pop")
    {
        genepop_input_c<1> input("genotype_genepop_format.txt");
        REQUIRE(input.Dist_btw_pop[0][0] == 0);
        REQUIRE(input.Dist_btw_pop[0][1] == 1);
        REQUIRE(input.Dist_btw_pop[0][2] == 2);
        REQUIRE(input.Dist_btw_pop[0][3] == 3);
        REQUIRE(input.Dist_btw_pop[0][4] == 4);
    }

    SECTION("regres_ar_without_missing_value")
    {
        genepop_input_c<2> input("test_ar.txt");
        data_plane_vec_c data_plane_vec(input);

        auto ar = ar_by_pair(data_plane_vec);
        auto regr_ar_loc = linear_regres_X_Y(ar);

        REQUIRE(Approx(regr_ar_loc.at(0)).margin(0.000001) == -0.0184381);
        REQUIRE(Approx(regr_ar_loc.at(1)).margin(0.000001) == 0.238358);
    }

    SECTION("regres_ar_with_missing_value")
    {
        genepop_input_c<2> input("test_ar_mv.txt");
        data_plane_vec_c data_plane_vec(input);

        auto ar = ar_by_pair(data_plane_vec);
        auto regr_ar_loc = linear_regres_X_Y(ar);

        REQUIRE(Approx(regr_ar_loc.at(0)).margin(0.000001) == -0.0130383);
        REQUIRE(Approx(regr_ar_loc.at(1)).margin(0.00001) == 0.23328);
    }
}

TEST_CASE("reptile_test")
{
    SECTION("read_test")
    {
        genepop_input_c<2> input("Test_reptile_genepop.txt");
        //Verify ADH-4 , ADH-5 \n mtDNA have been well separate
        REQUIRE(input.Locus_name.size() == 6043);
        //pop name
        REQUIRE(input.Pop_name.size() == 186);
        REQUIRE(input.Pop_name[0][0] == "2.65613");
        REQUIRE(input.Pop_name[0][1] == "14.7019");
        //Name of indiv
        REQUIRE(input.Indiv_name.size() == 186);
        //Pop size
        REQUIRE(input.Genotype.size() == 186);
        REQUIRE(input.Genotype[0].size() == 1);
        //First indiv
        REQUIRE(input.Genotype[0][0][0] == std::array{12, 12});
        REQUIRE(input.Genotype[0][0][1] == std::array{10, 10});
        REQUIRE(input.Genotype[0][0].at(6042) == std::array{0, 0});
        //Last indiv -1
        REQUIRE(input.Genotype[185][0][0] == std::array{10, 12});
        REQUIRE(input.Genotype.at(185)[0].at(6042) == std::array{12, 12});
    }

    SECTION("Fst_by_dist")
    {
        genepop_input_c<2> input("Test_reptile_genepop.txt");

        data_plane_vec_c data_plane_vec(input);
        auto result = Fstat_by_loc_with_probid(data_plane_vec, 0);
        std::cout<<"RESULTAT : "<<result.at(1)<<std::endl;
        REQUIRE(0);
    }
}