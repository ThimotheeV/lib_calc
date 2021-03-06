#define CATCH_CONFIG_MAIN CalcStatTest
#include "catch.hpp"

#include "calc_stat.hpp"

TEST_CASE("Q_intra_indiv Q_inter_indiv_intra_deme Q_inter_deme calc_stat_test")
{
    SECTION("Q_intra_indiv without missing value")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        REQUIRE(calc_Q_intra_indiv(data_plane_vec) == Approx(5.0 / 18).margin(0.000001));
    }
    SECTION("Q_intra_indiv with missing value")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 0}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {0, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {0, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        REQUIRE(calc_Q_intra_indiv(data_plane_vec) == Approx(5.0 / 15).margin(0.000001));
    }
    SECTION("Q_inter_indiv_intra_deme without missing value")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 3 indiv, 2 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}}, {{1, 1}, {1, 1}}, {{1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}}, {{1, 3}, {2, 3}}, {{2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}}, {{1, 3}, {2, 1}}, {{2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        auto result = calc_Q_inter_indiv_intra_deme(data_plane_vec);
        REQUIRE(result == Approx(28.0 / 72).margin(0.000001));
    }
    SECTION("Q_inter_indiv_intra_deme with missing value")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 3 indiv, 2 locus
        genepop_input.Genotype = {{{{1, 0}, {1, 1}}, {{1, 1}, {1, 1}}, {{0, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}}, {{0, 3}, {2, 3}}, {{2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}}, {{1, 3}, {2, 1}}, {{2, 2}, {2, 0}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        auto result = calc_Q_inter_indiv_intra_deme(data_plane_vec);
        double mean = ((2.0 / 5 + 3.0 / 8 + 2.0 / 12) / 3 + (8.0 / 12 + 4.0 / 12 + 4.0 / 8) / 3) / 2;
        REQUIRE(result == mean);
    }
    SECTION("calc_Q_inter_deme_per_chr_per_locus without missing value")
    {
        genepop_input_c<2> genepop_input;
        //2 deme, 4 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {2, 2}, {1, 1}}, {{1, 2}, {1, 2}, {1, 3}}, {{1, 1}, {1, 1}, {1, 1}}, {{1, 1}, {2, 2}, {3, 3}}},
                                  {{{2, 3}, {1, 2}, {2, 3}}, {{3, 4}, {3, 1}, {4, 3}}, {{1, 2}, {2, 5}, {5, 4}}, {{1, 1}, {4, 5}, {1, 1}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        auto result = calc_Q_inter_deme_per_chr_per_locus(data_plane_vec, 0, 0);
        REQUIRE(result == 22.0 / 64);
        result = calc_Q_inter_deme_per_chr_per_locus(data_plane_vec, 1, 0);
        REQUIRE(result == 16.0 / 64);
        result = calc_Q_inter_deme_per_chr_per_locus(data_plane_vec, 2, 0);
        REQUIRE(result == 16.0 / 64);
    }
    // SECTION("Q_inter_deme with missing value")
    // {
    //     genepop_input_c<2> genepop_input;
    //     //3 deme, 3 indiv, 2 locus
    //     genepop_input.Genotype = {{{{1, 0}, {1, 1}}, {{1, 1}, {1, 1}}, {{0, 2}, {1, 3}}},
    //                               {{{1, 3}, {1, 3}}, {{0, 3}, {2, 3}}, {{2, 3}, {2, 3}}},
    //                               {{{1, 1}, {1, 2}}, {{1, 3}, {2, 1}}, {{2, 2}, {2, 0}}}};

    //     data_plane_vec_c data_plane_vec(genepop_input);
    //     auto result = calc_Q_inter_deme(data_plane_vec);
    //     REQUIRE(result == 53.0 / 170);
    // }
    // SECTION("Q_inter_deme with diff deme size")
    // {
    //     genepop_input_c<2> genepop_input;
    //     //3 deme, 2-1-3 indiv, 3 locus
    //     genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
    //                               {{{1, 3}, {1, 3}, {1, 3}}},
    //                               {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}, {{2, 3}, {2, 3}, {2, 3}}}};

    //     data_plane_vec_c data_plane_vec(genepop_input);
    //     auto result = calc_Q_inter_deme(data_plane_vec);
    //     REQUIRE(result == 44.0 / 132);
    // }
}
TEST_CASE("calc_Hnei_calc_stat_test")
{
    //Compare with IBDSim
    SECTION("calc_Hnei one locus")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //3 deme, 3 indiv, 1 locus
        genepop_input.Genotype = {{{{1, 2}}, {{3, 1}}, {{4, 5}}},
                                  {{{6, 1}}, {{5, 7}}, {{8, 9}}},
                                  {{{10, 11}}, {{12, 1}}, {{5, 13}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        auto result = calc_Hnei_per_chr_per_loc(data_plane_vec, 0, 0);
        REQUIRE(result == Approx(0.941176).margin(0.0001));
    }

    SECTION("calc_Hnei multi locus")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //3 deme, 3 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 2}, {1, 2}}, {{3, 1}, {3, 2}, {3, 4}}, {{4, 5}, {4, 5}, {5, 6}}},
                                  {{{6, 1}, {6, 7}, {7, 8}}, {{5, 7}, {1, 8}, {8, 9}}, {{8, 9}, {9, 1}, {10, 11}}},
                                  {{{10, 11}, {4, 1}, {1, 11}}, {{12, 1}, {8, 5}, {6, 12}}, {{5, 13}, {10, 1}, {13, 14}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        auto result = calc_Hnei(data_plane_vec);
        REQUIRE(result == Approx(0.941176).margin(0.0001));
    }
}

TEST_CASE("calc_Var_calc_stat_test")
{
    SECTION("calc_Var simple case")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2 indiv, 1 locus
        genepop_input.Genotype = {{{{10, 12}}, {{10, 8}}},
                                  {{{10, 8}}, {{12, 12}}},
                                  {{{10, 10}}, {{8, 10}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        REQUIRE(calc_Var_per_chr_per_loc(data_plane_vec, 0, 0) == (24.0 / 12) * (12.0 / 11));
    }

    SECTION("calc_Var complexe case")
    {
        genepop_input_c<2> genepop_input;
        //2 deme, 2 indiv, 1 locus
        genepop_input.Genotype = {{{{1, 2}}, {{2, 1}}},
                                  {{{5, 4}}, {{9, 7}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        REQUIRE(calc_Var_per_chr_per_loc(data_plane_vec, 0, 0) == Approx(8.6964).margin(0.0001));
    }
}

TEST_CASE("calc_Hobs_calc_stat_test")
{
    SECTION("calc_Hobs without missing value")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        REQUIRE(calc_Hobs_per_chr_per_loc(data_plane_vec, 0, 0) == 4.0 / 6);
        REQUIRE(calc_Hobs_per_chr_per_loc(data_plane_vec, 1, 0) == 4.0 / 6);
        REQUIRE(calc_Hobs_per_chr_per_loc(data_plane_vec, 2, 0) == 5.0 / 6);
    }
}

TEST_CASE("calc_MGW_calc_stat_test")
{
    //Compare with IBDSim
    SECTION("calc_MGW one locus")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //3 deme, 3 indiv, 1 locus
        genepop_input.Genotype = {{{{1, 2}}, {{3, 1}}, {{4, 5}}},
                                  {{{6, 1}}, {{5, 7}}, {{8, 6}}},
                                  {{{10, 8}}, {{10, 1}}, {{5, 10}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        auto result = calc_MGW_per_chr_per_loc(data_plane_vec, 0, 0);
        REQUIRE(result == 0.9);
    }

    SECTION("calc_MGW multi locus")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //3 deme, 3 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 2}, {1, 2}}, {{3, 1}, {3, 2}, {3, 4}}, {{4, 5}, {4, 5}, {5, 6}}},
                                  {{{6, 1}, {6, 7}, {7, 8}}, {{5, 7}, {1, 8}, {8, 9}}, {{8, 9}, {9, 1}, {10, 11}}},
                                  {{{10, 11}, {4, 1}, {1, 11}}, {{12, 1}, {8, 5}, {6, 12}}, {{5, 13}, {10, 1}, {13, 14}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        auto result = calc_MGW(data_plane_vec);
        REQUIRE(result == 1);
    }
}

TEST_CASE("qr_all_loc_calc_stat_test")
{
    SECTION("qr by locus with one indiv haploid per deme")
    {
        genepop_input_c<1> genepop_input;
        genepop_input.Geo_dist_class_nbr = 3;
        genepop_input.Geo_dist_class_btw_deme = {{0, 1, 2},
                                                 {1, 0, 1},
                                                 {2, 1, 0}};
        //3 deme, 1 indiv, 3 locus
        genepop_input.Genotype = {{{{2}, {1}, {1}}},
                                  {{{1}, {3}, {1}}},
                                  {{{1}, {2}, {1}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        auto result = calc_qr_per_chr_by_loc(data_plane_vec, 0, 0);
        REQUIRE(result[1] == 1.0 / 2);
        REQUIRE(result[2] == 0.0 / 1);

        result = calc_qr_per_chr_by_loc(data_plane_vec, 1, 0);
        REQUIRE(result[1] == 0.0 / 2);
        REQUIRE(result[2] == 0.0 / 1);

        result = calc_qr_per_chr_by_loc(data_plane_vec, 2, 0);
        REQUIRE(result[1] == 2.0 / 2);
        REQUIRE(result[2] == 1.0 / 1);
    }

    SECTION("qr by locus without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_class_nbr = 3;
        genepop_input.Geo_dist_class_btw_deme = {{0, 1, 2},
                                                 {1, 0, 1},
                                                 {2, 1, 0}};
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        auto result = calc_qr_per_chr_by_loc(data_plane_vec, 0, 0);
        REQUIRE(result[0] == 5.0 / 12);
        REQUIRE(result[1] == 8.0 / 32);
        REQUIRE(result[2] == 10.0 / 16);

        result = calc_qr_per_chr_by_loc(data_plane_vec, 1, 0);
        REQUIRE(result[0] == 5.0 / 12);
        REQUIRE(result[1] == 8.0 / 32);
        REQUIRE(result[2] == 6.0 / 16);

        result = calc_qr_per_chr_by_loc(data_plane_vec, 2, 0);
        REQUIRE(result[0] == 4.0 / 12);
        REQUIRE(result[1] == 11.0 / 32);
        REQUIRE(result[2] == 5.0 / 16);
    }

    SECTION("qr without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_class_nbr = 3;
        genepop_input.Geo_dist_class_btw_deme = {{0, 1, 2},
                                                 {1, 0, 1},
                                                 {2, 1, 0}};
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        auto result = calc_qr(data_plane_vec);
        REQUIRE(result[0] == 14.0 / 36);
        REQUIRE(result[1] == 27.0 / 96);
        REQUIRE(result[2] == 21.0 / 48);
    }

    SECTION("haploid qr with diff deme size")
    {
        genepop_input_c<1> genepop_input;
        genepop_input.Geo_dist_class_nbr = 3;
        genepop_input.Geo_dist_class_btw_deme = {{0, 1, 2},
                                                 {1, 0, 1},
                                                 {2, 1, 0}};
        //3 deme, 2-1-3 indiv, 3 locus
        genepop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {1}, {1}}},
                                  {{{1}, {1}, {1}}},
                                  {{{1}, {1}, {1}}, {{2}, {2}, {2}}, {{2}, {2}, {2}}}};

        data_plane_vec_c data_plane_vec(genepop_input);
        auto result = calc_qr(data_plane_vec);
        REQUIRE(result[0] == 6.0 / 12);
        REQUIRE(result[1] == 9.0 / 15);
        REQUIRE(result[2] == 6.0 / 18);
    }
}

TEST_CASE("ar_calc_stat_test")
{
    SECTION("diploid ar without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        std::vector<std::array<double, 2>> expectation = {{0, -0.1538}, {1, -0.1538}, {1, 0.3076}, {2, -0.1538}, {2, 0.3076}, {1, -0.2692}, {1, 0.0769}, {2, -0.5}, {2, -0.1538}, {0, -0.1538}, {1, -0.2692}, {1, 0.0769}, {1, 0.0769}, {1, -0.2692}, {0, -0.1538}};

        auto result = ar_by_indiv_pair(data_plane_vec);
        auto result_itr = result.begin();
        for (auto const &exp_value : expectation)
        {
            REQUIRE(exp_value.at(0) == result_itr->at(0));
            REQUIRE(exp_value.at(1) == Approx(result_itr->at(1)).margin(0.0001));
            ++result_itr;
        }
    }

    SECTION("diploid ar with missing value")
    {
        genepop_input_c<2> genepop_input; //3 deme, 2-1-3 indiv, 3 locus, 3 chr
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}, {{2, 3}, {2, 3}, {2, 3}}}};
        genepop_input.Geo_dist_btw_deme = {{0}};
        //1 deme, 4 indiv, 3 locus
        genepop_input.Genotype = {{{{0, 2}, {1, 1}, {1, 2}}, {{1, 2}, {1, 0}, {1, 3}}, {{1, 1}, {2, 3}, {0, 3}}, {{2, 2}, {3, 3}, {3, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        std::vector<std::array<double, 2>> expectation = {{0, -0.125}, {0, 1.75}, {0, 1.25}, {0, 0.25}, {0, 0}, {0, 1.375}};

        auto result = ar_by_indiv_pair(data_plane_vec);
        auto result_itr = result.begin();
        for (auto const &exp_value : expectation)
        {
            REQUIRE(exp_value.at(0) == result_itr->at(0));
            REQUIRE(exp_value.at(1) == Approx(result_itr->at(1)).margin(0.0001));
            ++result_itr;
        }
    }
}

TEST_CASE("er_calc_stat_test")
{
    SECTION("simple er without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //3 deme, 1 indiv, 1 locus
        genepop_input.Genotype = {{{{1, 1}}},
                                  {{{1, 2}}},
                                  {{{2, 2}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        std::vector<std::array<double, 2>> expectation = {{1, -0.25}, {2, 1.25}, {1, -0.25}};

        auto result = er_by_indiv_pair(data_plane_vec);
        auto result_itr = result.begin();
        for (auto const &exp_value : expectation)
        {
            REQUIRE(exp_value.at(0) == result_itr->at(0));
            REQUIRE(exp_value.at(1) == Approx(result_itr->at(1)).margin(0.0001));
            ++result_itr;
        }
    }
    //Resul from IBDSim
    SECTION("er without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2, 3},
                                           {1, 0, 1, 2},
                                           {2, 1, 0, 1},
                                           {3, 2, 1, 0}};
        //4 deme, 1 indiv, 2 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}}},
                                  {{{1, 3}, {1, 1}}},
                                  {{{3, 3}, {2, 1}}},
                                  {{{1, 1}, {1, 1}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        std::vector<std::array<double, 2>> expectation = {{1, -0.033334}, {2, 0.383333}, {3, -0.283334}, {1, -0.116667}, {2, -0.116667}, {1, 0.633333}};

        auto result = er_by_indiv_pair(data_plane_vec);
        auto result_itr = result.begin();
        for (auto const &exp_value : expectation)
        {
            REQUIRE(exp_value.at(0) == result_itr->at(0));
            REQUIRE(exp_value.at(1) == Approx(result_itr->at(1)).margin(0.0001));
            ++result_itr;
        }
    }

    SECTION("er without missing value, bis")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 2}, {1, 1}}, {{3, 3}, {3, 4}, {1, 2}}},
                                  {{{3, 4}, {5, 6}, {2, 3}}, {{1, 5}, {5, 7}, {1, 1}}},
                                  {{{6, 7}, {6, 3}, {1, 2}}, {{2, 8}, {5, 8}, {1, 1}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        std::vector<std::array<double, 2>> expectation = {{0, 0.1224}, {1, 0.2117}, {1, -0.1454}, {2, 0.0688}, {2, -0.1454}, {1, -0.1454}, {1, 0.1403}, {2, -0.0739}, {2, 0.1403}, {0, 0.1224}, {1, -0.0918}, {1, 0.1224}, {1, 0.0867}, {1, -0.1275}, {0, 0.0867}};

        auto result = er_by_indiv_pair(data_plane_vec);
        auto result_itr = result.begin();
        for (auto const &exp_value : expectation)
        {
            REQUIRE(exp_value.at(0) == result_itr->at(0));
            REQUIRE(exp_value.at(1) == Approx(result_itr->at(1)).margin(0.0001));
            ++result_itr;
        }
    }

    SECTION("diploid er with missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0}};
        //1 deme, 4 indiv, 3 locus
        genepop_input.Genotype = {{{{0, 2}, {1, 1}, {1, 2}},
                                   {{1, 2}, {1, 0}, {1, 3}},
                                   {{1, 1}, {2, 3}, {0, 3}},
                                   {{2, 2}, {3, 3}, {3, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        std::vector<std::array<double, 2>> expectation = {{0, -0.0625}, {0, 0.5}, {0, 0.625}, {0, -0.25}, {0, -0.125}, {0, 0.375}};

        auto result = er_by_indiv_pair(data_plane_vec);
        auto result_itr = result.begin();
        for (auto const &exp_value : expectation)
        {
            REQUIRE(exp_value.at(0) == result_itr->at(0));
            REQUIRE(exp_value.at(1) == Approx(result_itr->at(1)).margin(0.0001));
            ++result_itr;
        }
    }
}

TEST_CASE("Fstat_per_chr_by_loc_with_probid_calc_stat_test")
{
    SECTION("Fstat_per_chr_by_loc_with_probid simple case")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //3 deme, 2 indiv, 1 locus
        genepop_input.Genotype = {{{{1, 1}}, {{1, 1}}},
                                  {{{1, 2}}, {{1, 2}}},
                                  {{{2, 2}}, {{2, 2}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result = Fstat_per_chr_by_loc_with_probid(data_plane_vec, 0, 0);
        //FiS
        REQUIRE(result.at(0).at(0) / result.at(0).at(1) == Approx(-1).margin(0.1));
        //FST
        REQUIRE(result.at(1).at(0) / result.at(1).at(1) == Approx(0.75).margin(0.01));
    }

    SECTION("Fstat_per_chr_by_loc_with_probid without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2, 3},
                                           {1, 0, 1, 2},
                                           {2, 1, 0, 1},
                                           {3, 2, 1, 0}};
        //4 deme, 2 indiv, 1 locus
        genepop_input.Genotype = {{{{1, 1}}, {{1, 2}}},
                                  {{{1, 2}}, {{2, 2}}},
                                  {{{2, 3}}, {{2, 3}}},
                                  {{{3, 3}}, {{2, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result = Fstat_per_chr_by_loc_with_probid(data_plane_vec, 0, 0);
        //FiS
        REQUIRE(result.at(0).at(0) / result.at(0).at(1) == -0.250000);
        //FST
        REQUIRE(result.at(1).at(0) / result.at(1).at(1) == Approx(0.314286).margin(0.000001));
    }
}

TEST_CASE("Fstat_per_chr_by_loc_with_indic_calc_stat_test")
{
    SECTION("Fstat_per_chr_by_loc_with_indic simple case")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //3 deme, 2 indiv, 1 locus, 1 chr
        genepop_input.Genotype = {{{{1, 1}}, {{1, 1}}},
                                  {{{1, 2}}, {{1, 2}}},
                                  {{{2, 2}}, {{2, 2}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result = Fstat_per_chr_by_loc_with_indic(data_plane_vec, 0, 0);
        //FiS
        REQUIRE(result.at(0).at(0) / result.at(0).at(1) == Approx(-1).margin(0.1));
        //FST
        REQUIRE(result.at(1).at(0) / result.at(1).at(1) == Approx(0.75).margin(0.01));
    }

    SECTION("Fstat_per_chr_by_loc_with_indic without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2, 3},
                                           {1, 0, 1, 2},
                                           {2, 1, 0, 1},
                                           {3, 2, 1, 0}};
        //4 deme, 2 indiv, 1 locus, 1 chr
        genepop_input.Genotype = {{{{1, 1}}, {{1, 2}}},
                                  {{{1, 2}}, {{2, 2}}},
                                  {{{2, 3}}, {{2, 3}}},
                                  {{{3, 3}}, {{2, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result = Fstat_per_chr_by_loc_with_indic(data_plane_vec, 0, 0);
        //FiS
        REQUIRE(result.at(0).at(0) / result.at(0).at(1) == -0.250000);
        //FST
        REQUIRE(result.at(1).at(0) / result.at(1).at(1) == Approx(0.314286).margin(0.000001));
    }

    SECTION("Fstat_per_chr_by_loc_with_indic simple case with missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //3 deme, 2 indiv, 1 locus
        genepop_input.Genotype = {{{{1, 1}}, {{1, 1}}},
                                  {{{1, 2}}, {{1, 2}}},
                                  {{{0, 2}}, {{2, 2}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result = Fstat_per_chr_by_loc_with_indic(data_plane_vec, 0, 0);
        //FiS
        REQUIRE(result.at(0).at(0) / result.at(0).at(1) == Approx(-1).margin(0.1));
        //FST
        REQUIRE(result.at(1).at(0) / result.at(1).at(1) == Approx(0.686275).margin(0.000001));
    }

    SECTION("Fstat_per_chr_by_loc_with_indic with missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2, 3},
                                           {1, 0, 1, 2},
                                           {2, 1, 0, 1},
                                           {3, 2, 1, 0}};
        //4 deme, 2 indiv, 1 locus
        genepop_input.Genotype = {{{{1, 0}}, {{1, 2}}},
                                  {{{1, 2}}, {{2, 2}}},
                                  {{{2, 3}}, {{0, 3}}},
                                  {{{3, 3}}, {{2, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result = Fstat_per_chr_by_loc_with_indic(data_plane_vec, 0, 0);
        //FiS
        REQUIRE(result.at(0).at(0) / result.at(0).at(1) == Approx(-0.142857).margin(0.000001));
        //FST
        REQUIRE(result.at(1).at(0) / result.at(1).at(1) == Approx(0.141509).margin(0.000001));
    }
}

TEST_CASE("Fstat_genepop")
{
    SECTION("Fstat without missing value with_probid")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result_0 = Fstat_per_chr_by_loc_with_probid(data_plane_vec, 0, 0);
        REQUIRE(result_0.at(0).at(0) / result_0.at(0).at(1) == Approx(-0.1429).margin(0.0001));
        REQUIRE(result_0.at(1).at(0) / result_0.at(1).at(1) == Approx(0.0667).margin(0.0001));

        auto result_1 = Fstat_per_chr_by_loc_with_probid(data_plane_vec, 1, 0);
        REQUIRE(result_1.at(0).at(0) / result_1.at(0).at(1) == Approx(-0.1429).margin(0.0001));
        REQUIRE(result_1.at(1).at(0) / result_1.at(1).at(1) == Approx(0.1765).margin(0.0001));

        auto result_2 = Fstat_per_chr_by_loc_with_probid(data_plane_vec, 2, 0);
        REQUIRE(result_2.at(0).at(0) / result_2.at(0).at(1) == Approx(-0.2500).margin(0.0001));
        REQUIRE(result_2.at(1).at(0) / result_2.at(1).at(1) == Approx(0).margin(0.0001));

        auto result = Fstat_genepop(data_plane_vec, false);
        REQUIRE(result.at(0) == Approx(-0.1818).margin(0.0001));
        REQUIRE(result.at(1) == Approx(0.0833).margin(0.0001));
    }

    SECTION("Fstat without missing value with_indic")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result_0 = Fstat_per_chr_by_loc_with_indic(data_plane_vec, 0, 0);
        REQUIRE(result_0.at(0).at(0) / result_0.at(0).at(1) == Approx(-0.1429).margin(0.0001));
        REQUIRE(result_0.at(1).at(0) / result_0.at(1).at(1) == Approx(0.0667).margin(0.0001));

        auto result_1 = Fstat_per_chr_by_loc_with_indic(data_plane_vec, 1, 0);
        REQUIRE(result_1.at(0).at(0) / result_1.at(0).at(1) == Approx(-0.1429).margin(0.0001));
        REQUIRE(result_1.at(1).at(0) / result_1.at(1).at(1) == Approx(0.1765).margin(0.0001));

        auto result_2 = Fstat_per_chr_by_loc_with_indic(data_plane_vec, 2, 0);
        REQUIRE(result_2.at(0).at(0) / result_2.at(0).at(1) == Approx(-0.2500).margin(0.0001));
        REQUIRE(result_2.at(1).at(0) / result_2.at(1).at(1) == Approx(0).margin(0.0001));

        auto result = Fstat_genepop(data_plane_vec, true);
        REQUIRE(result.at(0) == Approx(-0.1818).margin(0.0001));
        REQUIRE(result.at(1) == Approx(0.0833).margin(0.0001));
    }

    SECTION("Fstat with missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{0, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {0, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 0}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 0}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result_0 = Fstat_per_chr_by_loc_with_indic(data_plane_vec, 0, 0);
        REQUIRE(result_0.at(0).at(0) / result_0.at(0).at(1) == Approx(-0.0909).margin(0.0001));
        REQUIRE(result_0.at(1).at(0) / result_0.at(1).at(1) == Approx(0.1456).margin(0.0001));

        auto result_1 = Fstat_per_chr_by_loc_with_indic(data_plane_vec, 1, 0);
        REQUIRE(result_1.at(0).at(0) / result_1.at(0).at(1) == Approx(-0.0909).margin(0.0001));
        REQUIRE(result_1.at(1).at(0) / result_1.at(1).at(1) == Approx(0.1852).margin(0.0001));

        auto result_2 = Fstat_per_chr_by_loc_with_indic(data_plane_vec, 2, 0);
        REQUIRE(result_2.at(0).at(0) / result_2.at(0).at(1) == Approx(-0.2000).margin(0.0001));
        REQUIRE(result_2.at(1).at(0) / result_2.at(1).at(1) == Approx(-0.3158).margin(0.0001));

        auto result = Fstat_genepop(data_plane_vec, true);
        REQUIRE(result.at(0) == Approx(-0.1244).margin(0.0001));
        REQUIRE(result.at(1) == Approx(0.0601).margin(0.0001));
    }
}

TEST_CASE("AFS")
{
    SECTION("AFS 2 states without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //1 deme, 6 indiv, 3 locus; 1 ancestry state, 2 derived state
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 2}}, {{1, 2}, {1, 2}, {1, 2}}, {{2, 2}, {2, 2}, {2, 2}}, {{1, 1}, {1, 2}, {1, 2}}, {{2, 1}, {2, 2}, {2, 2}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result = calc_AFS(data_plane_vec, 0);

        //for frequency : 5,7 /12
        std::vector<int> target{1, 2};

        auto target_itr = target.begin();

        for (auto const &pair : result.at(2))
        {
            REQUIRE(pair.second == *target_itr);
            ++target_itr;
        }
    }

    SECTION("AFS 2 states with missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};
        //1 deme, 3 indiv, 3 locus; 1 ancestry state, 2 derived state
        genepop_input.Genotype = {{{{0, 2}, {1, 1}, {1, 2}}, {{1, 2}, {1, 0}, {1, 2}}, {{0, 1}, {2, 2}, {2, 2}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result = calc_AFS(data_plane_vec, 0);

        //for frequency : 5,7 /12
        std::vector<double> target{0.066, 0.533, 2, 0.4};

        auto target_itr = target.begin();
        auto const &map_state = result.at(1);

        for (auto const &pair : map_state)
        {
            REQUIRE(pair.second == Approx(*target_itr).margin(0.001));
            ++target_itr;
        }
    }
}

TEST_CASE("lin_Fst")
{
    SECTION("lin_Fst diploide without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};

        //3 deme, 4 indiv, 3 locus, 3 chr
        genepop_input.Genotype = {
            {{{1, 2}, {2, 2}, {1, 1}}, {{1, 2}, {1, 2}, {1, 3}}, {{1, 1}, {1, 1}, {1, 1}}, {{1, 1}, {2, 2}, {3, 3}}},
            {{{1, 3}, {2, 2}, {1, 3}}, {{2, 3}, {4, 3}, {2, 3}}, {{1, 3}, {4, 3}, {1, 3}}, {{1, 3}, {1, 1}, {1, 3}}},
            {{{2, 3}, {1, 2}, {2, 3}}, {{3, 4}, {3, 1}, {4, 3}}, {{1, 2}, {2, 5}, {5, 4}}, {{1, 1}, {4, 5}, {1, 1}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result = lin_Fst_by_deme_pair(data_plane_vec);

        //See lint_Fst in golden master file
        std::vector<std::array<double, 2>> target{{1, 0.051470588235294}, {2, -0.014705882352941}, {1, -0.073529411764706}};

        auto target_itr = target.begin();

        for (auto const &value : result)
        {
            REQUIRE(value.at(1) == Approx(target_itr->at(1)).margin(0.0000001));
            ++target_itr;
        }
    }

    SECTION("lin_Fst haploide without missing value")
    {
        genepop_input_c<1> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};

        //3 deme, 4 indiv, 3 locus, 3 chr
        genepop_input.Genotype = {
            {{{1}, {2}, {1}}, {{1}, {1}, {1}}, {{1}, {1}, {1}}, {{1}, {2}, {3}}},
            {{{1}, {2}, {1}}, {{2}, {4}, {2}}, {{1}, {4}, {1}}, {{1}, {1}, {1}}},
            {{{2}, {1}, {2}}, {{3}, {3}, {4}}, {{1}, {2}, {5}}, {{1}, {4}, {1}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result = lin_Fst_by_deme_pair(data_plane_vec);

        //See lint_Fst in golden master file
        std::vector<std::array<double, 2>> target{{1, 0.260714285714286}, {2, -0.060714285714286}, {1, -0.060714285714286}};

        auto target_itr = target.begin();

        for (auto const &value : result)
        {
            REQUIRE(value.at(1) == Approx(target_itr->at(1)).margin(0.0000001));
            ++target_itr;
        }
    }

    SECTION("lin_Fst diploide with missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};

        //3 deme, 4 indiv, 3 locus, 3 chr
        genepop_input.Genotype = {
            {{{0, 2}, {2, 2}, {1, 1}}, {{1, 2}, {0, 2}, {1, 3}}, {{1, 1}, {1, 1}, {1, 1}}, {{1, 1}, {2, 2}, {3, 3}}},
            {{{1, 3}, {2, 2}, {1, 3}}, {{2, 3}, {4, 3}, {2, 3}}, {{1, 3}, {4, 3}, {0, 3}}, {{1, 3}, {1, 1}, {1, 3}}},
            {{{2, 3}, {1, 2}, {2, 3}}, {{0, 4}, {3, 1}, {4, 3}}, {{1, 2}, {2, 5}, {5, 4}}, {{1, 1}, {4, 5}, {1, 1}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result = lin_Fst_by_deme_pair(data_plane_vec);

        //See lint_Fst in golden master file
        std::vector<std::array<double, 2>> target{{1, 0.039888682745826}, {2, 0.012456930824278}, {1, -0.081168831168831}};

        auto target_itr = target.begin();

        for (auto const &value : result)
        {
            REQUIRE(value.at(1) == Approx(target_itr->at(1)).margin(0.0000001));
            ++target_itr;
        }
    }

    SECTION("lin_Fst haploide with missing value")
    {
        genepop_input_c<1> genepop_input;
        genepop_input.Geo_dist_btw_deme = {{0, 1, 2},
                                           {1, 0, 1},
                                           {2, 1, 0}};

        //3 deme, 4 indiv, 3 locus, 3 chr
        genepop_input.Genotype = {
            {{{0}, {2}, {1}}, {{1}, {0}, {1}}, {{1}, {1}, {1}}, {{1}, {2}, {3}}},
            {{{1}, {2}, {1}}, {{2}, {4}, {2}}, {{1}, {4}, {0}}, {{1}, {1}, {1}}},
            {{{2}, {1}, {2}}, {{0}, {3}, {4}}, {{1}, {2}, {5}}, {{1}, {4}, {1}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        auto result = lin_Fst_by_deme_pair(data_plane_vec);

        //See lint_Fst in golden master file
        std::vector<std::array<double, 2>> target{{1, 0.228571428571429}, {2, 0.025}, {1, 0.014285714285714}};

        auto target_itr = target.begin();

        for (auto const &value : result)
        {
            REQUIRE(value.at(1) == Approx(target_itr->at(1)).margin(0.0000001));
            ++target_itr;
        }
    }
}