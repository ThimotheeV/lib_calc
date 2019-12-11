#define CATCH_CONFIG_MAIN CustomVecTest
#include <catch2/catch.hpp>

#include "calc_stat.hpp"

TEST_CASE("Q0_Q1_Q2_calc_stat_test")
{
    SECTION("Q0 without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        REQUIRE(calc_Q0(data_plane_vec) == 5.0 / 18);
    }
    SECTION("Q1 without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_Q1(data_plane_vec);
        REQUIRE(result == 19.0 / 54);
    }
    SECTION("Q2 without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_Q2(data_plane_vec);
        REQUIRE(result == 48.0 / 144);
    }
    SECTION("Q2 with diff pop size")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2-1-3 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}, {{2, 3}, {2, 3}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_Q2(data_plane_vec);
        REQUIRE(result == 44.0 / 132);
    }

    SECTION("Fst without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_Fstat(calc_Q1(data_plane_vec), calc_Q2(data_plane_vec));
        REQUIRE(result == Approx(0.02777777).margin(0.00000001));
    }
}

TEST_CASE("calc_qstat_test")
{
    SECTION("qstat one locus without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_qstat_loc_by_loc(data_plane_vec, 0);
        REQUIRE(result.Q0_intra_ind == Approx(2.0 / 6).margin(0.00000001));
        REQUIRE(result.Q1_intra_pop == 19.0 / 54);
        REQUIRE(result.Q2_inter_pop == 48.0 / 144);

        result = calc_qstat_loc_by_loc(data_plane_vec, 1);
        REQUIRE(result.Q0_intra_ind == Approx(6.0 / 18).margin(0.00000001));
        REQUIRE(result.Q1_intra_pop == 19.0 / 54);
        REQUIRE(result.Q2_inter_pop == 48.0 / 144);

        result = calc_qstat_loc_by_loc(data_plane_vec, 2);
        REQUIRE(result.Q0_intra_ind == Approx(3.0 / 18).margin(0.00000001));
        REQUIRE(result.Q1_intra_pop == 19.0 / 54);
        REQUIRE(result.Q2_inter_pop == 48.0 / 144);
    }

    SECTION("qstat without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_qstat_all_loc(data_plane_vec);
        REQUIRE(result.Q0_intra_ind == Approx(5.0 / 18).margin(0.00000001));
        REQUIRE(result.Q1_intra_pop == 19.0 / 54);
        REQUIRE(result.Q2_inter_pop == 48.0 / 144);
    }

    SECTION("qstat with diff pop size")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2-1-3 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}, {{2, 3}, {2, 3}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_qstat_all_loc(data_plane_vec);
        REQUIRE(result.Q2_inter_pop == 44.0 / 132);
    }

    SECTION("haploid qstat with diff pop size")
    {
        genepop_input_c<1> gen_pop_input;
        //3 pop, 2-1-3 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {1}, {1}}},
                                  {{{1}, {1}, {1}}},
                                  {{{1}, {1}, {1}}, {{2}, {2}, {2}}, {{2}, {2}, {2}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_qstat_all_loc(data_plane_vec);
        REQUIRE(result.Q2_inter_pop == 15.0 / 33);
    }
}

TEST_CASE("calc_qr_all_loc_test")
{
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
        REQUIRE(result[0] == 19.0 / 54);
        REQUIRE(result[1] == 27.0 / 96);
        REQUIRE(result[2] == 21.0 / 48);
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
        REQUIRE(result[0] == 6.0 / 12);
        REQUIRE(result[1] == 9.0 / 15);
        REQUIRE(result[2] == 6.0 / 18);
    }

    SECTION("Fst_qstat_all_loc without missing value")
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
        auto Fst = Fst_qstat_all_loc(result);
        REQUIRE(Fst[0] == 0);
        REQUIRE(Fst[1] == Approx(0.098228).margin(0.000001));
        REQUIRE(Fst[2] == Approx(-0.152263).margin(0.000001));
    }

    SECTION("haploid Fst_qstat_all_loc with diff pop size")
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
        auto Fst = Fst_qstat_all_loc(result);
        REQUIRE(Fst[0] == 0);
        REQUIRE(Fst[1] == Approx(-0.25).margin(0.000001));
        REQUIRE(Fst[2] == Approx(0.25).margin(0.000001));
    }
}

TEST_CASE("mean_var_cov_regr_test")
{
    SECTION("double mean(std::vector<value> X)")
    {
        std::vector<int> vec_i = {479, 720, 63, 435, 234, 967, 608, 1, 270, 991,
                                  220, 279, 836, 117, 498, 300, 101, 269, 983, 575,
                                  923, 401, 37, 197, 398, 355, 228, 865, 307, 897};
        auto res = mean(vec_i);
        REQUIRE(res == 451.8);

        std::vector<double> vec_d = {0.971, 0.564, 0.329, 0.962, 0.803, 0.894, 0.178, 0.709, 0.633, 0.607,
                                     0.917, 0.810, 0.410, 0.952, 0.234, 0.708, 0.168, 0.648, 0.066, 0.902,
                                     0.616, 0.426, 0.839, 0.397, 0.682, 0.398, 0.940, 0.190, 0.852, 0.822};
        res = mean(vec_d);
        REQUIRE(res == Approx(0.6209).margin(0.0001));
    }

    SECTION("double var(std::vector<value> X, double meanX)")
    {
        std::vector<int> vec_i = {479, 720, 63, 435, 234, 967, 608, 1, 270, 991,
                                  220, 279, 836, 117, 498, 300, 101, 269, 983, 575,
                                  923, 401, 37, 197, 398, 355, 228, 865, 307, 897};
        auto meanX = mean(vec_i);
        auto res = var(vec_i, meanX);
        REQUIRE(res == Approx(94151.293333).margin(0.000001));

        std::vector<double> vec_d = {0.971, 0.564, 0.329, 0.962, 0.803, 0.894, 0.178, 0.709, 0.633, 0.607,
                                     0.917, 0.810, 0.410, 0.952, 0.234, 0.708, 0.168, 0.648, 0.066, 0.902,
                                     0.616, 0.426, 0.839, 0.397, 0.682, 0.398, 0.940, 0.190, 0.852, 0.822};
        meanX = mean(vec_d);
        res = var(vec_d, meanX);
        REQUIRE(res == 0.07458149);
    }

    SECTION("double cov_X_Y(std::vector<value1> X_vec, double meanX, std::vector<value2> Y_vec, double meanY)")
    {
        std::vector<int> vec_i = {479, 720, 63, 435, 234, 967, 608, 1, 270, 991,
                                  220, 279, 836, 117, 498, 300, 101, 269, 983, 575,
                                  923, 401, 37, 197, 398, 355, 228, 865, 307, 897};

        std::vector<double> vec_d = {0.971, 0.564, 0.329, 0.962, 0.803, 0.894, 0.178, 0.709, 0.633, 0.607,
                                     0.917, 0.810, 0.410, 0.952, 0.234, 0.708, 0.168, 0.648, 0.066, 0.902,
                                     0.616, 0.426, 0.839, 0.397, 0.682, 0.398, 0.940, 0.190, 0.852, 0.822};
        auto meanX = mean(vec_i);
        auto meanY = mean(vec_d);
        auto res = cov_X_Y(vec_i, meanX, vec_d, meanY);
        REQUIRE(res == Approx(-18.231953).margin(0.000001));
    }

    SECTION("double linear_regres_X_Y(std::vector<value1> X_vec, std::vector<value2> Y_vec)")
    {
        std::vector<int> vec_i = {479, 720, 63, 435, 234, 967, 608, 1, 270, 991,
                                  220, 279, 836, 117, 498, 300, 101, 269, 983, 575,
                                  923, 401, 37, 197, 398, 355, 228, 865, 307, 897};

        std::vector<double> vec_d = {0.971, 0.564, 0.329, 0.962, 0.803, 0.894, 0.178, 0.709, 0.633, 0.607,
                                     0.917, 0.810, 0.410, 0.952, 0.234, 0.708, 0.168, 0.648, 0.066, 0.902,
                                     0.616, 0.426, 0.839, 0.397, 0.682, 0.398, 0.940, 0.190, 0.852, 0.822};
        auto res = linear_regres_X_Y(vec_i, vec_d);
        REQUIRE(res.at(0) == Approx(-0.0001936).margin(0.0000001));
        REQUIRE(res.at(1) == Approx(0.7083889).margin(0.0000001));
    }
}
