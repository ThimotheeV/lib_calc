#define CATCH_CONFIG_MAIN CalcStatTest
#include <catch2/catch.hpp>

#include "calc_stat.hpp"

TEST_CASE("Q_intra_indiv Q_inter_indiv_intra_pop Q_inter_pop calc_stat_test")
{
    SECTION("Q_intra_indiv without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        REQUIRE(calc_Q_intra_indiv(data_plane_vec) == 5.0 / 18);
    }
    SECTION("Q_intra_indiv with missing value")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 0}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {0, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {0, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        REQUIRE(calc_Q_intra_indiv(data_plane_vec) == 5.0 / 15);
    }
    SECTION("Q_inter_indiv_intra_pop without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 3 indiv, 2 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}}, {{1, 1}, {1, 1}}, {{1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}}, {{1, 3}, {2, 3}}, {{2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}}, {{1, 3}, {2, 1}}, {{2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_Q_inter_indiv_intra_pop(data_plane_vec);
        REQUIRE(result == 28.0 / 72);
    }
    SECTION("Q_inter_indiv_intra_pop with missing value")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 3 indiv, 2 locus
        gen_pop_input.Genotype = {{{{1, 0}, {1, 1}}, {{1, 1}, {1, 1}}, {{0, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}}, {{0, 3}, {2, 3}}, {{2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}}, {{1, 3}, {2, 1}}, {{2, 2}, {2, 0}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_Q_inter_indiv_intra_pop(data_plane_vec);
        REQUIRE(result == 23.0 / 57);
    }
    // SECTION("Q_inter_pop without missing value")
    // {
    //     genepop_input_c<2> gen_pop_input;
    //     //3 pop, 3 indiv, 2 locus
    //     gen_pop_input.Genotype = {{{{1, 2}, {1, 1}}, {{1, 1}, {1, 1}}, {{1, 2}, {1, 3}}},
    //                               {{{1, 3}, {1, 3}}, {{1, 3}, {2, 3}}, {{2, 3}, {2, 3}}},
    //                               {{{1, 1}, {1, 2}}, {{1, 3}, {2, 1}}, {{2, 2}, {2, 3}}}};

    //     data_plane_vec_c data_plane_vec(gen_pop_input);
    //     auto result = calc_Q_inter_pop(data_plane_vec);
    //     REQUIRE(result == 67.0 / 216);
    // }
    SECTION("Q_inter_pop with missing value")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 3 indiv, 2 locus
        gen_pop_input.Genotype = {{{{1, 0}, {1, 1}}, {{1, 1}, {1, 1}}, {{0, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}}, {{0, 3}, {2, 3}}, {{2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}}, {{1, 3}, {2, 1}}, {{2, 2}, {2, 0}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_Q_inter_pop(data_plane_vec);
        REQUIRE(result == 53.0 / 170);
    }
    SECTION("Q_inter_pop with diff pop size")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2-1-3 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}, {{2, 3}, {2, 3}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_Q_inter_pop(data_plane_vec);
        REQUIRE(result == 44.0 / 132);
    }
}

TEST_CASE("qr_all_loc_test")
{
    SECTION("qr by locus without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2},
                                      {1, 0, 1},
                                      {2, 1, 0}};
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_qr_loc_by_loc(data_plane_vec, 2, 0);
        REQUIRE(result[0] == std::array<double, 2>{0, 7.0 / 18});
        REQUIRE(result[1] == std::array<double, 2>{1, 8.0 / 32});
        REQUIRE(result[2] == std::array<double, 2>{2, 10.0 / 16});

        result = calc_qr_loc_by_loc(data_plane_vec, 2, 1);
        REQUIRE(result[0] == std::array<double, 2>{0, 7.0 / 18});
        REQUIRE(result[1] == std::array<double, 2>{1, 8.0 / 32});
        REQUIRE(result[2] == std::array<double, 2>{2, 6.0 / 16});

        result = calc_qr_loc_by_loc(data_plane_vec, 2, 2);
        REQUIRE(result[0] == std::array<double, 2>{0, 5.0 / 18});
        REQUIRE(result[1] == std::array<double, 2>{1, 11.0 / 32});
        REQUIRE(result[2] == std::array<double, 2>{2, 5.0 / 16});
    }

    SECTION("qr without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2},
                                      {1, 0, 1},
                                      {2, 1, 0}};
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_qr_all_loc(data_plane_vec, 2);
        REQUIRE(result[0] == std::array<double, 2>{0, 19.0 / 54});
        REQUIRE(result[1] == std::array<double, 2>{1, 27.0 / 96});
        REQUIRE(result[2] == std::array<double, 2>{2, 21.0 / 48});
    }

    SECTION("haploid qr with diff pop size")
    {
        genepop_input_c<1> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2},
                                      {1, 0, 1},
                                      {2, 1, 0}};
        //3 pop, 2-1-3 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {1}, {1}}},
                                  {{{1}, {1}, {1}}},
                                  {{{1}, {1}, {1}}, {{2}, {2}, {2}}, {{2}, {2}, {2}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_qr_all_loc(data_plane_vec, 2);
        REQUIRE(result[0] == std::array<double, 2>{0, 6.0 / 12});
        REQUIRE(result[1] == std::array<double, 2>{1, 9.0 / 15});
        REQUIRE(result[2] == std::array<double, 2>{2, 6.0 / 18});
    }
}

TEST_CASE("ar_test")
{
    SECTION("diploid ar without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2},
                                      {1, 0, 1},
                                      {2, 1, 0}};
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);

        std::vector<std::array<double, 2>> expectation = {{0, -0.1538}, {1, -0.1538}, {1, 0.3076}, {2, -0.1538}, {2, 0.3076}, {1, -0.2692}, {1, 0.0769}, {2, -0.5}, {2, -0.1538}, {0, -0.1538}, {1, -0.2692}, {1, 0.0769}, {1, 0.0769}, {1, -0.2692}, {0, -0.1538}};

        auto result = ar_by_pair(data_plane_vec);
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
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0}};
        //1 pop, 4 indiv, 3 locus
        gen_pop_input.Genotype = {{{{0, 2}, {1, 1}, {1, 2}}, {{1, 2}, {1, 0}, {1, 3}}, {{1, 1}, {2, 3}, {0, 3}}, {{2, 2}, {3, 3}, {3, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);

        std::vector<std::array<double, 2>> expectation = {{0, -0.125}, {0, 1.75}, {0, 1.25}, {0, 0.25}, {0, 0}, {0, 1.375}};

        auto result = ar_by_pair(data_plane_vec);
        auto result_itr = result.begin();
        for (auto const &exp_value : expectation)
        {
            REQUIRE(exp_value.at(0) == result_itr->at(0));
            REQUIRE(exp_value.at(1) == Approx(result_itr->at(1)).margin(0.0001));
            ++result_itr;
        }
    }
}

TEST_CASE("er_test")
{
    SECTION("simple er without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2},
                                      {1, 0, 1},
                                      {2, 1, 0}};
        //3 pop, 1 indiv, 1 locus
        gen_pop_input.Genotype = {{{{1, 1}}},
                                  {{{1, 2}}},
                                  {{{2, 2}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);

        std::vector<std::array<double, 2>> expectation = {{1, -0.25}, {2, 1.25}, {1, -0.25}};

        auto result = er_by_pair(data_plane_vec);
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
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2, 3},
                                      {1, 0, 1, 2},
                                      {2, 1, 0, 1},
                                      {3, 2, 1, 0}};
        //4 pop, 1 indiv, 2 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}}},
                                  {{{1, 3}, {1, 1}}},
                                  {{{3, 3}, {2, 1}}},
                                  {{{1, 1}, {1, 1}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);

        std::vector<std::array<double, 2>> expectation = {{1, -0.033334}, {2, 0.383333}, {3, -0.283334}, {1, -0.116667}, {2, -0.116667}, {1, 0.633333}};

        auto result = er_by_pair(data_plane_vec);
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
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2},
                                      {1, 0, 1},
                                      {2, 1, 0}};
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 2}, {1, 1}}, {{3, 3}, {3, 4}, {1, 2}}},
                                  {{{3, 4}, {5, 6}, {2, 3}}, {{1, 5}, {5, 7}, {1, 1}}},
                                  {{{6, 7}, {6, 3}, {1, 2}}, {{2, 8}, {5, 8}, {1, 1}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);

        std::vector<std::array<double, 2>> expectation = {{0, 0.1224}, {1, 0.2117}, {1, -0.1454}, {2, 0.0688}, {2, -0.1454}, {1, -0.1454}, {1, 0.1403}, {2, -0.0739}, {2, 0.1403}, {0, 0.1224}, {1, -0.0918}, {1, 0.1224}, {1, 0.0867}, {1, -0.1275}, {0, 0.0867}};

        auto result = er_by_pair(data_plane_vec);
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
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0}};
        //1 pop, 4 indiv, 3 locus
        gen_pop_input.Genotype = {{{{0, 2}, {1, 1}, {1, 2}},
                                   {{1, 2}, {1, 0}, {1, 3}},
                                   {{1, 1}, {2, 3}, {0, 3}},
                                   {{2, 2}, {3, 3}, {3, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);

        std::vector<std::array<double, 2>> expectation = {{0, -0.0625}, {0, 0.5}, {0, 0.625}, {0, -0.25}, {0, -0.125}, {0, 0.375}};

        auto result = er_by_pair(data_plane_vec);
        auto result_itr = result.begin();
        for (auto const &exp_value : expectation)
        {
            REQUIRE(exp_value.at(0) == result_itr->at(0));
            REQUIRE(exp_value.at(1) == Approx(result_itr->at(1)).margin(0.0001));
            ++result_itr;
        }
    }
}

TEST_CASE("Fstat_by_loc_with_probid_test")
{
    SECTION("Fstat_by_loc_with_probid simple case")
    {
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2},
                                      {1, 0, 1},
                                      {2, 1, 0}};
        //3 pop, 2 indiv, 1 locus
        gen_pop_input.Genotype = {{{{1, 1}}, {{1, 1}}},
                                  {{{1, 2}}, {{1, 2}}},
                                  {{{2, 2}}, {{2, 2}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);

        auto result = Fstat_by_loc_with_probid(data_plane_vec, 0);
        //FiS
        REQUIRE(result.at(0).at(0) / result.at(0).at(1) == Approx(-1).margin(0.1));
        //FST
        REQUIRE(result.at(1).at(0) / result.at(1).at(1) == Approx(0.75).margin(0.01));
    }

    SECTION("Fstat_by_loc_with_probid without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2, 3},
                                      {1, 0, 1, 2},
                                      {2, 1, 0, 1},
                                      {3, 2, 1, 0}};
        //4 pop, 2 indiv, 1 locus
        gen_pop_input.Genotype = {{{{1, 1}}, {{1, 2}}},
                                  {{{1, 2}}, {{2, 2}}},
                                  {{{2, 3}}, {{2, 3}}},
                                  {{{3, 3}}, {{2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);

        auto result = Fstat_by_loc_with_probid(data_plane_vec, 0);
        //FiS
        REQUIRE(result.at(0).at(0) / result.at(0).at(1) == -0.250000);
        //FST
        REQUIRE(result.at(1).at(0) / result.at(1).at(1) == Approx(0.314286).margin(0.000001));
    }
}

TEST_CASE("Fstat_by_loc_with_indic_test")
{
    SECTION("Fstat_by_loc_with_indic simple case")
    {
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2},
                                      {1, 0, 1},
                                      {2, 1, 0}};
        //3 pop, 2 indiv, 1 locus
        gen_pop_input.Genotype = {{{{1, 1}}, {{1, 1}}},
                                  {{{1, 2}}, {{1, 2}}},
                                  {{{2, 2}}, {{2, 2}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);

        auto result = Fstat_by_loc_with_indic(data_plane_vec, 0);
        //FiS
        REQUIRE(result.at(0).at(0) / result.at(0).at(1) == Approx(-1).margin(0.1));
        //FST
        REQUIRE(result.at(1).at(0) / result.at(1).at(1) == Approx(0.75).margin(0.01));
    }

    SECTION("Fstat_by_loc_with_indic without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2, 3},
                                      {1, 0, 1, 2},
                                      {2, 1, 0, 1},
                                      {3, 2, 1, 0}};
        //4 pop, 2 indiv, 1 locus
        gen_pop_input.Genotype = {{{{1, 1}}, {{1, 2}}},
                                  {{{1, 2}}, {{2, 2}}},
                                  {{{2, 3}}, {{2, 3}}},
                                  {{{3, 3}}, {{2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);

        auto result = Fstat_by_loc_with_indic(data_plane_vec, 0);
        //FiS
        REQUIRE(result.at(0).at(0) / result.at(0).at(1) == -0.250000);
        //FST
        REQUIRE(result.at(1).at(0) / result.at(1).at(1) == Approx(0.314286).margin(0.000001));
    }

    SECTION("Fstat_by_loc_with_indic simple case with missing value")
    {
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2},
                                      {1, 0, 1},
                                      {2, 1, 0}};
        //3 pop, 2 indiv, 1 locus
        gen_pop_input.Genotype = {{{{1, 1}}, {{1, 1}}},
                                  {{{1, 2}}, {{1, 2}}},
                                  {{{0, 2}}, {{2, 2}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);

        auto result = Fstat_by_loc_with_indic(data_plane_vec, 0);
        //FiS
        REQUIRE(result.at(0).at(0) / result.at(0).at(1) == Approx(-1).margin(0.1));
        //FST
        REQUIRE(result.at(1).at(0) / result.at(1).at(1) == Approx(0.686275).margin(0.000001));
    }

    SECTION("Fstat_by_loc_with_indic with missing value")
    {
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2, 3},
                                      {1, 0, 1, 2},
                                      {2, 1, 0, 1},
                                      {3, 2, 1, 0}};
        //4 pop, 2 indiv, 1 locus
        gen_pop_input.Genotype = {{{{1, 0}}, {{1, 2}}},
                                  {{{1, 2}}, {{2, 2}}},
                                  {{{2, 3}}, {{0, 3}}},
                                  {{{3, 3}}, {{2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);

        auto result = Fstat_by_loc_with_indic(data_plane_vec, 0);
        //FiS
        REQUIRE(result.at(0).at(0) / result.at(0).at(1) == Approx(-0.142857).margin(0.000001));
        //FST
        REQUIRE(result.at(1).at(0) / result.at(1).at(1) == Approx(0.141509).margin(0.000001));
    }
}

TEST_CASE("Fstat_genepop")
{
    SECTION("Fstat without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2},
                                      {1, 0, 1},
                                      {2, 1, 0}};
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);

        auto result_0 = Fstat_by_loc_with_indic(data_plane_vec, 0);
        REQUIRE(result_0.at(0).at(0) / result_0.at(0).at(1) == Approx(-0.1429).margin(0.0001));
        REQUIRE(result_0.at(1).at(0) / result_0.at(1).at(1) == Approx(0.0667).margin(0.0001));

        auto result_1 = Fstat_by_loc_with_indic(data_plane_vec, 1);
        REQUIRE(result_1.at(0).at(0) / result_1.at(0).at(1) == Approx(-0.1429).margin(0.0001));
        REQUIRE(result_1.at(1).at(0) / result_1.at(1).at(1) == Approx(0.1765).margin(0.0001));

        auto result_2 = Fstat_by_loc_with_indic(data_plane_vec, 2);
        REQUIRE(result_2.at(0).at(0) / result_2.at(0).at(1) == Approx(-0.2500).margin(0.0001));
        REQUIRE(result_2.at(1).at(0) / result_2.at(1).at(1) == Approx(0).margin(0.0001));

        auto result = Fstat_genepop(data_plane_vec);
        REQUIRE(result.at(0) == Approx(-0.1818).margin(0.0001));
        REQUIRE(result.at(1) == Approx(0.0833).margin(0.0001));
    }

    SECTION("Fstat with missing value")
    {
        genepop_input_c<2> gen_pop_input;
        gen_pop_input.Dist_btw_pop = {{0, 1, 2},
                                      {1, 0, 1},
                                      {2, 1, 0}};
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{0, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {0, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 0}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 0}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);

        auto result_0 = Fstat_by_loc_with_indic(data_plane_vec, 0);
        REQUIRE(result_0.at(0).at(0) / result_0.at(0).at(1) == Approx(-0.0909).margin(0.0001));
        REQUIRE(result_0.at(1).at(0) / result_0.at(1).at(1) == Approx(0.1456).margin(0.0001));

        auto result_1 = Fstat_by_loc_with_indic(data_plane_vec, 1);
        REQUIRE(result_1.at(0).at(0) / result_1.at(0).at(1) == Approx(-0.0909).margin(0.0001));
        REQUIRE(result_1.at(1).at(0) / result_1.at(1).at(1) == Approx(0.1852).margin(0.0001));

        auto result_2 = Fstat_by_loc_with_indic(data_plane_vec, 2);
        REQUIRE(result_2.at(0).at(0) / result_2.at(0).at(1) == Approx(-0.2000).margin(0.0001));
        REQUIRE(result_2.at(1).at(0) / result_2.at(1).at(1) == Approx(-0.3158).margin(0.0001));

        auto result = Fstat_genepop(data_plane_vec);
        REQUIRE(result.at(0) == Approx(-0.1244).margin(0.0001));
        REQUIRE(result.at(1) == Approx(0.0601).margin(0.0001));
    }
}