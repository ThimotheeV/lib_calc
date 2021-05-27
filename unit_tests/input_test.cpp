#define CATCH_CONFIG_MAIN InputTest
#include "catch.hpp"

#include "input.hpp"

using namespace gss;
TEST_CASE("input_handler_test")
{
    SECTION("remove_spaces_tab_in_range")
    {
        std::string str = "Grange des Peres  ";
        auto result = remove_spaces_tab_in_range(str, 15, str.size());
        result = str_tolower(result);

        REQUIRE(result == "grange des peres");
    }

    SECTION("remove_spaces_tab_underscores_andlower")
    {
        auto result = remove_spaces_tab_underscores("Grange_des Peres  ");
        result = str_tolower(result);

        REQUIRE(result == "grangedesperes");
    }

    SECTION("sep_by_space")
    {
        auto result = slice_by_char(" 0201  003003 0102   0302 1011 01 ", ' ');

        REQUIRE(result.size() == 6);
        REQUIRE(result[0] == "0201");
        REQUIRE(result[1] == "003003");
        REQUIRE(result[2] == "0102");
        REQUIRE(result[3] == "0302");
        REQUIRE(result[4] == "1011");
        REQUIRE(result[5] == "01");
    }

    SECTION("sep_by_comma")
    {
        auto result = slice_by_char(" Grange des Peres  ,  0201 003003 0102 0302 1011 01", ',');

        REQUIRE(result.size() == 2);
        REQUIRE(result[0] == " Grange des Peres  ");
        REQUIRE(result[1] == "  0201 003003 0102 0302 1011 01");
    }

    SECTION("trim_unix_file")
    {
        auto result = slice_unix_windows_file_by_line("line 1\n\nline 2\n\n");

        REQUIRE(result.size() == 2);
        REQUIRE(result[0] == "line 1");
        REQUIRE(result[1] == "line 2");
    }

    SECTION("trim_windows_file")
    {
        auto result = slice_unix_windows_file_by_line("line 1\r\rline 2\r\r");

        REQUIRE(result.size() == 2);
        REQUIRE(result[0] == "line 1");
        REQUIRE(result[1] == "line 2");
    }

    SECTION("trim_by_string")
    {
        auto result = slice_str_vec_by_string({"mtDNA", "Pop", "Grange des Peres  ,  0201 003003 0102 0302 1011 01"}, "Pop");

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

TEST_CASE("genepop_input_test")
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
        genepop_input_c<2> input("input_test.txt", 10);
        REQUIRE(input.Nbr_geo_dist_class == 10);
        REQUIRE(input.Nbr_chr_dist_class == 0);
        REQUIRE(input.Dist_btw_loc == std::vector<std::vector<std::vector<double>>>{});
        REQUIRE(input.Dist_class_btw_loc == std::vector<std::vector<std::vector<int>>>{});
        //Verify ADH-4 , ADH-5 \n mtDNA have been well separate
        REQUIRE(input.Locus_name.size() == 6);
        REQUIRE(input.Locus_name[0] == "adhlocus1");
        REQUIRE(input.Locus_name[4] == "adh-5");
        REQUIRE(input.Locus_name[5] == "mtdna");
        //deme name
        REQUIRE(input.Pop_name.size() == 4);
        REQUIRE(input.Pop_name[2][0] == "0");
        REQUIRE(input.Pop_name[2][1] == "2.00");
        REQUIRE(input.Pop_name[3][0] == "-0.00");
        REQUIRE(input.Pop_name[3][1] == "3");
        //Name of indiv
        REQUIRE(input.Indiv_name.size() == 4);
        REQUIRE(input.Indiv_name[0][0] == "grange des peres  ");
        REQUIRE(input.Indiv_name[2][1] == "bonneau 02   ");
        REQUIRE(input.Indiv_name[3][2] == "NaN");
        //Number of deme
        REQUIRE(input.Genotype.size() == 4);
        //Pop size
        REQUIRE(input.Genotype[0].size() == 3);
        //First indiv
        REQUIRE(input.Genotype[0][0].at(0) == std::array{2, 1});
        REQUIRE(input.Genotype[0][0].at(5) == std::array{1, 0});
        //Last indiv
        REQUIRE(input.Genotype[3][3].at(0) == std::array{01, 01});
        REQUIRE(input.Genotype[3][3].at(5) == std::array{2, 0});
    }

    SECTION("constructor_win_format")
    {
        genepop_input_c<2> input("input_test_win_form.txt");
        //Verify ADH-4 , ADH-5 \n mtDNA have been well separate
        REQUIRE(input.Locus_name.size() == 6);
        REQUIRE(input.Locus_name[4] == "adh-5");
        REQUIRE(input.Locus_name[5] == "mtdna");
        //deme name
        REQUIRE(input.Pop_name.size() == 4);
        REQUIRE(input.Pop_name[2][0] == "bonneau");
        REQUIRE(input.Pop_name[3][0] == "last");
        //Name of indiv
        REQUIRE(input.Indiv_name.size() == 4);
        REQUIRE(input.Indiv_name[0][0] == "grange des peres  ");
        REQUIRE(input.Indiv_name[2][1] == "bonneau 02   ");
        REQUIRE(input.Indiv_name[3][2] == "NaN");
        //Number of deme
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
    SECTION("geo_dist_btw_gene for int dist")
    {
        genepop_input_c<2> input("input_test.txt", 10);
        REQUIRE(input.Dist_btw_deme[0] == std::vector<double>{0, 1, 2, 3});
        REQUIRE(input.Dist_btw_deme[1] == std::vector<double>{1, 0, 1, 2});
        REQUIRE(input.Dist_btw_deme[2] == std::vector<double>{2, 1, 0, 1});
        REQUIRE(input.Dist_btw_deme[3] == std::vector<double>{3, 2, 1, 0});
        //10 classes
        REQUIRE(input.Dist_class_btw_deme[0] == std::vector<int>{0, 3, 6, 9});
        REQUIRE(input.Dist_class_btw_deme[1] == std::vector<int>{3, 0, 3, 6});
        REQUIRE(input.Dist_class_btw_deme[2] == std::vector<int>{6, 3, 0, 3});
        REQUIRE(input.Dist_class_btw_deme[3] == std::vector<int>{9, 6, 3, 0});
    }

    SECTION("geo_dist_btw_gene for int dist")
    {
        genepop_input_c<2> input("float_dist_test.txt", 3);
        REQUIRE(input.Nbr_geo_dist_class == 3);
        REQUIRE(input.Dist_btw_deme[0][0] == Approx(0).margin(0.01));
        REQUIRE(input.Dist_btw_deme[0][1] == Approx(1.5).margin(0.01));
        REQUIRE(input.Dist_btw_deme[0][2] == Approx(6.8).margin(0.01));
        REQUIRE(input.Dist_btw_deme[0][3] == Approx(8.1).margin(0.01));

        REQUIRE(input.Dist_btw_deme[1][0] == Approx(1.5).margin(0.01));
        REQUIRE(input.Dist_btw_deme[1][1] == Approx(0).margin(0.01));
        REQUIRE(input.Dist_btw_deme[1][2] == Approx(5.3).margin(0.01));
        REQUIRE(input.Dist_btw_deme[1][3] == Approx(6.6).margin(0.01));

        REQUIRE(input.Dist_btw_deme[2][0] == Approx(6.8).margin(0.01));
        REQUIRE(input.Dist_btw_deme[2][1] == Approx(5.3).margin(0.01));
        REQUIRE(input.Dist_btw_deme[2][2] == Approx(0).margin(0.01));
        REQUIRE(input.Dist_btw_deme[2][3] == Approx(1.3).margin(0.01));

        REQUIRE(input.Dist_btw_deme[3][0] == Approx(8.1).margin(0.01));
        REQUIRE(input.Dist_btw_deme[3][1] == Approx(6.6).margin(0.01));
        REQUIRE(input.Dist_btw_deme[3][2] == Approx(1.3).margin(0.01));
        REQUIRE(input.Dist_btw_deme[3][3] == Approx(0).margin(0.01));

        REQUIRE(input.Dist_class_btw_deme[0] == std::vector<int>{0, 0, 2, 2});
        REQUIRE(input.Dist_class_btw_deme[1] == std::vector<int>{0, 0, 1, 2});
        REQUIRE(input.Dist_class_btw_deme[2] == std::vector<int>{2, 1, 0, 0});
        REQUIRE(input.Dist_class_btw_deme[3] == std::vector<int>{2, 2, 0, 0});
    }

    SECTION("Dist_btw_loc")
    {
        genepop_input_c<2> input("input_test.txt", 0, "input_test.map", 0);

        REQUIRE(input.Dist_btw_loc[0][0] == std::vector<double>{0, 1, 2, 3});
        REQUIRE(input.Dist_btw_loc[0][1] == std::vector<double>{1, 0, 1, 2});
        REQUIRE(input.Dist_btw_loc[0][2] == std::vector<double>{2, 1, 0, 1});
        REQUIRE(input.Dist_btw_loc[0][3] == std::vector<double>{3, 2, 1, 0});
        REQUIRE(input.Dist_btw_loc[1][0] == std::vector<double>{0, 1});
        REQUIRE(input.Dist_btw_loc[1][1] == std::vector<double>{1, 0});
    }

    SECTION("Dist_class_btw_loc")
    {
        genepop_input_c<2> input("input_test.txt", 0, "input_test.map", 4);

        REQUIRE(input.Dist_class_btw_loc[0][0] == std::vector<int>{0, 1, 2, 3});
        REQUIRE(input.Dist_class_btw_loc[0][1] == std::vector<int>{1, 0, 1, 2});
        REQUIRE(input.Dist_class_btw_loc[0][2] == std::vector<int>{2, 1, 0, 1});
        REQUIRE(input.Dist_class_btw_loc[0][3] == std::vector<int>{3, 2, 1, 0});
        REQUIRE(input.Dist_class_btw_loc[1][0] == std::vector<int>{0, 1});
        REQUIRE(input.Dist_class_btw_loc[1][1] == std::vector<int>{1, 0});

        genepop_input_c<2> input1("input_test.txt", 0, "input_test.map", 2);

        REQUIRE(input1.Dist_class_btw_loc[0][0] == std::vector<int>{0, 0, 1, 1});
        REQUIRE(input1.Dist_class_btw_loc[0][1] == std::vector<int>{0, 0, 0, 1});
        REQUIRE(input1.Dist_class_btw_loc[0][2] == std::vector<int>{1, 0, 0, 0});
        REQUIRE(input1.Dist_class_btw_loc[0][3] == std::vector<int>{1, 1, 0, 0});
        REQUIRE(input1.Dist_class_btw_loc[1][0] == std::vector<int>{0, 0});
        REQUIRE(input1.Dist_class_btw_loc[1][1] == std::vector<int>{0, 0});
    }
}

TEST_CASE("selector_input_test")
{
    SECTION("constructor")
    {
        selector_input_c setting("test_selector.txt");

        REQUIRE(setting.Data_filename == "./test.txt");
        REQUIRE(setting.Genetic_map_name == "./test.map");

        REQUIRE(setting.Hobs == true);
        REQUIRE(setting.Hexp == true);
        REQUIRE(setting.Nb_allele == false);
        REQUIRE(setting.Var == false);
        REQUIRE(setting.MGW == false);

        REQUIRE(setting.F_stat == true);
        REQUIRE(setting.Q_stat == false);

        REQUIRE(setting.Qr == true);
        REQUIRE(setting.Ar == false);
        REQUIRE(setting.Er == false);

        REQUIRE(setting.Eta == true);

        REQUIRE(setting.Missing_data == false);
        REQUIRE(setting.Nbr_geo_dist_class == 100);
    }
}

TEST_CASE("map_input_test")
{
    SECTION("read")
    {
        auto const map_file = gss::str_tolower(gss::read_file("input_test.map"));
        auto const map_file_vec = gss::slice_unix_windows_file_by_line(map_file);

        std::vector<std::vector<std::string>> crude_map_vec(map_file_vec.size());
        for (int i = 0; i < crude_map_vec.size(); ++i)
        {
            crude_map_vec[i] = gss::slice_by_char(map_file_vec[i], '\t');
            REQUIRE(crude_map_vec[i].size() == 4);
        }
    }
}