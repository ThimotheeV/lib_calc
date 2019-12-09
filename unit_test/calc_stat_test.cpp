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
    SECTION("qstat without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_qstat(data_plane_vec);
        REQUIRE(result.Q0_intra_ind == 5.0 / 18);
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
        auto result = calc_qstat(data_plane_vec);
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
        auto result = calc_qstat(data_plane_vec);
        REQUIRE(result.Q2_inter_pop == 15.0 / 33);
    }
}

TEST_CASE("calc_qr_test")
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
        auto result = calc_qr(data_plane_vec, 2);
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
        auto result = calc_qr(data_plane_vec, 2);
        REQUIRE(result[0] == 6.0 / 12);
        REQUIRE(result[1] == 9.0 / 15);
        REQUIRE(result[2] == 6.0 / 18);
    }
}