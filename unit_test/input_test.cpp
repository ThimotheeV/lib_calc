#define CATCH_CONFIG_MAIN InputTest
#include <catch2/catch.hpp>

#include "input.hpp"

TEST_CASE("input_handler_test")
{
    SECTION("remove_spaces_in_range")
    {
        std::string str = "Grange des Peres  ";
        auto result = remove_spaces_in_range(str, 15, str.size());
        result = lower_case(result);

        REQUIRE(result == "grange des peres");
    }

    SECTION("remove_spaces_underscores_andlower")
    {
        auto result = remove_spaces_underscores("Grange_des Peres  ");
        result = lower_case(result);

        REQUIRE(result == "grangedesperes");
    }

    SECTION("trim_by_space")
    {
        auto result = trim_by_char(" 0201  003003 0102   0302 1011 01", ' ');

        REQUIRE(result.size() == 6);
        REQUIRE(result[0] == "0201");
        REQUIRE(result[1] == "003003");
        REQUIRE(result[2] == "0102");
        REQUIRE(result[3] == "0302");
        REQUIRE(result[4] == "1011");
        REQUIRE(result[5] == "01");
    }

    SECTION("trim_by_comma")
    {
        auto result = trim_by_char("Grange des Peres  ,  0201 003003 0102 0302 1011 01", ',');

        REQUIRE(result[0] == "Grange des Peres  ");
        REQUIRE(result[1] == "  0201 003003 0102 0302 1011 01");
    }

    SECTION("trim_unix_file")
    {
        auto result = trim_unix_windows_file_by_line("line 1\nline 2");

        REQUIRE(result.size() == 2);
        REQUIRE(result[0] == "line 1");
        REQUIRE(result[1] == "line 2");
    }

    SECTION("trim_windows_file")
    {
        auto result = trim_unix_windows_file_by_line("line 1\rline 2");

        REQUIRE(result.size() == 2);
        REQUIRE(result[0] == "line 1");
        REQUIRE(result[1] == "line 2");
    }

    SECTION("trim_by_string")
    {
        auto result = trim_str_vec_by_string({"mtDNA", "Pop", "Grange des Peres  ,  0201 003003 0102 0302 1011 01"}, "Pop");

        REQUIRE(result.size() == 2);
        REQUIRE(result[0][0] == "mtDNA");
        REQUIRE(result[1][0] == "Grange des Peres  ,  0201 003003 0102 0302 1011 01");
    }

    SECTION("read_file")
    {
        auto result = read_file("./read_test.txt");
        REQUIRE(result == "OK");

        REQUIRE_THROWS(read_file("./unknow"));
    }
}

TEST_CASE("gen_pop_input_test")
{
    SECTION("empty_constructor")
    {
        genepop_input_c<2> input;

        REQUIRE(input.Genotype.size() == 0);
    }

    SECTION("trim_locus")
    {
        genepop_input_c<2> input;

        REQUIRE(input.trim_locus("003003") == std::array{3, 3});
        REQUIRE(input.trim_locus("0201") == std::array{2, 1});
        REQUIRE(input.trim_locus("01") == std::array{1, 0});
    }

    SECTION("constructor")
    {
        genepop_input_c<2> input("input_test.txt");
        //Verify ADH-4 , ADH-5 \n mtDNA have been well separate
        REQUIRE(input.Locus_name.size() == 6);
        REQUIRE(input.Locus_name[4] == "adh-5");
        REQUIRE(input.Locus_name[5] == "mtdna");
        //pop name
        REQUIRE(input.Pop_name.size() == 4);
        REQUIRE(input.Pop_name[2][0] == "0");
        REQUIRE(input.Pop_name[2][1] == "2.00");
        REQUIRE(input.Pop_name[3][0] == "-0.00");
        REQUIRE(input.Pop_name[3][1] == "3");
        //Name of indiv
        REQUIRE(input.Indiv_name.size() == 4);
        REQUIRE(input.Indiv_name[0][0] == "grange des peres");
        REQUIRE(input.Indiv_name[2][1] == "bonneau 02");
        REQUIRE(input.Indiv_name[3][2] == "");
        //Number of pop
        REQUIRE(input.Genotype.size() == 4);
        //Pop size
        REQUIRE(input.Genotype[0].size() == 5);
        //First indiv
        REQUIRE(input.Genotype[0][0].at(0) == std::array{2, 1});
        REQUIRE(input.Genotype[0][0].at(5) == std::array{1, 0});
        //Last indiv -1
        REQUIRE(input.Genotype[3][2].at(0) == std::array{0, 10});
        REQUIRE(input.Genotype[3][2].at(4) == std::array{8, 7});
    }

    SECTION("constructor_win_format")
    {
        genepop_input_c<2> input("input_test_win_form.txt");
        //Verify ADH-4 , ADH-5 \n mtDNA have been well separate
        REQUIRE(input.Locus_name.size() == 6);
        REQUIRE(input.Locus_name[4] == "adh-5");
        REQUIRE(input.Locus_name[5] == "mtdna");
        //pop name
        REQUIRE(input.Pop_name.size() == 4);
        REQUIRE(input.Pop_name[2][0] == "bonneau");
        REQUIRE(input.Pop_name[3][0] == "last");
        //Name of indiv
        REQUIRE(input.Indiv_name.size() == 4);
        REQUIRE(input.Indiv_name[0][0] == "grange des peres");
        REQUIRE(input.Indiv_name[2][1] == "bonneau 02");
        REQUIRE(input.Indiv_name[3][2] == "");
        //Number of pop
        REQUIRE(input.Genotype.size() == 4);
        //Pop size
        REQUIRE(input.Genotype[0].size() == 5);
        //First indiv
        REQUIRE(input.Genotype[0][0].at(0) == std::array{2, 1});
        REQUIRE(input.Genotype[0][0].at(5) == std::array{1, 0});
        //Last indiv -1
        REQUIRE(input.Genotype[3][2].at(0) == std::array{0, 10});
        REQUIRE(input.Genotype[3][2].at(4) == std::array{8, 7});
    }
    SECTION("dist_btw_pop")
    {
        genepop_input_c<2> input("input_test.txt");
        REQUIRE(input.Dist_btw_pop[0] == std::vector<float>{0, 1, 2, 3});
        REQUIRE(input.Dist_btw_pop[1] == std::vector<float>{1, 0, 1, 2});
        REQUIRE(input.Dist_btw_pop[2] == std::vector<float>{2, 1, 0, 1});
        REQUIRE(input.Dist_btw_pop[3] == std::vector<float>{3, 2, 1, 0});
    }
}