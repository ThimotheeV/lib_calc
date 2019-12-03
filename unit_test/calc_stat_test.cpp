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
        REQUIRE(calc_Q0(data_plane_vec) == std::array<int, 2>{5, 18});
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
        REQUIRE(result == std::array<int, 2>{19, 54});
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
        REQUIRE(result == std::array<int, 2>{48, 144});
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
        REQUIRE(result == std::array<int, 2>{44, 132});
    }

    SECTION("Fst without missing value")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(gen_pop_input);
        auto result = calc_Fst(calc_Q1(data_plane_vec), calc_Q2(data_plane_vec));
        REQUIRE(static_cast<double>(result.at(0)) / result.at(1) == Approx(0.02777777).margin(0.00000001));
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
        REQUIRE(result.Q0 == std::array<int, 2>{5, 18});
        REQUIRE(result.Q1 == std::array<int, 2>{19, 54});
        REQUIRE(result.Q2 == std::array<int, 2>{48, 144});
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
        REQUIRE(result.Q2 == std::array<int, 2>{44, 132});
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
        REQUIRE(result.Q2 == std::array<int, 2>{15, 33});
    }
}