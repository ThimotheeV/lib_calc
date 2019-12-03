#define CATCH_CONFIG_MAIN CustomVecTest
#include <catch2/catch.hpp>

#include "custom_vec.hpp"

TEST_CASE("haploid_custom_vec_test")
{
    SECTION("constructor")
    {
        genepop_input_c<1> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}, {{3}, {2}, {2}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(gen_pop_input);
        //3 locus, 3 pop, 2 indiv
        //{{{{1}, {1}}, {{3}, {3}}, {{1}, {2}}},
        // {{{1}, {2}}, {{1}, {2}}, {{2}, {2}}},
        // {{{1}, {3}}, {{3}, {2}}, {{3}, {2}}}};
        REQUIRE(plane_3d_vec.pop_nbr() == 3);
        REQUIRE(plane_3d_vec.indiv_nbr_per_pop(0) == 2);
        REQUIRE(plane_3d_vec.locus_nbr() == 3);
        REQUIRE(plane_3d_vec.indiv_nbr() == 6);

        REQUIRE(plane_3d_vec.cumul_indiv_nbr_per_pop() == std::vector<int>{0, 2, 4});

        std::vector<int> result = {1, 1, 3, 3, 1, 2, 1, 2, 1, 2, 2, 2, 1, 3, 3, 2, 3, 2};
        REQUIRE(plane_3d_vec.get_plane_vec() == result);
    }

    SECTION("operator()(int pop, int indiv, int locus)")
    {
        genepop_input_c<1> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}, {{3}, {2}, {2}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(gen_pop_input);
        //3 locus, 3 pop, 2 indiv
        REQUIRE(plane_3d_vec(0, 0, 0, 0) == 1);
        REQUIRE(plane_3d_vec(1, 1, 0, 0) == 1);
        REQUIRE(plane_3d_vec(2, 2, 1, 0) == 2);
    }

    SECTION("index_begin_locus && index_end_locus")
    {
        genepop_input_c<1> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}, {{3}, {2}, {2}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(gen_pop_input);
        //3 locus, 3 pop, 2 indiv
        REQUIRE(plane_3d_vec.index_begin_locus(0) == 0);
        REQUIRE(plane_3d_vec.index_end_locus(0) == 6);
        REQUIRE(plane_3d_vec.index_begin_locus(2) == 12);
        REQUIRE(plane_3d_vec.index_end_locus(2) == 18);
    }

    SECTION("same_indiv with same size pop")
    {
        genepop_input_c<1> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}, {{3}, {2}, {2}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(gen_pop_input);
        //3 locus, 3 pop, 2 indiv
        //{{{{1}, {1}}, {{3}, {3}}, {{1}, {2}}},
        // {{{1}, {2}}, {{1}, {2}}, {{2}, {2}}},
        // {{{1}, {3}}, {{3}, {2}}, {{3}, {2}}}};
        //Same locus && same indiv ?
        REQUIRE(plane_3d_vec.same_indiv(0, 1) == false);
        REQUIRE(plane_3d_vec.same_indiv(0, 2) == false);
        REQUIRE(plane_3d_vec.same_indiv(0, 4) == false);
        REQUIRE(plane_3d_vec.same_indiv(0, 12) == false);

        REQUIRE(plane_3d_vec.same_indiv(4, 1) == false);
        REQUIRE(plane_3d_vec.same_indiv(4, 5) == false);
        REQUIRE(plane_3d_vec.same_indiv(4, 6) == false);
        REQUIRE(plane_3d_vec.same_indiv(4, 16) == false);
    }
    SECTION("same_pop with same size pop")
    {
        genepop_input_c<1> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}, {{3}, {2}, {2}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(gen_pop_input);
        //3 locus, 3 pop, 2 indiv
        //{{{{1}, {1}}, {{3}, {3}}, {{1}, {2}}},
        // {{{1}, {2}}, {{1}, {2}}, {{2}, {2}}},
        // {{{1}, {3}}, {{3}, {2}}, {{3}, {2}}}};

        //Same locus && same pop ?
        REQUIRE(plane_3d_vec.same_pop(0, 0, 1) == true);
        REQUIRE(plane_3d_vec.same_pop(0, 0, 2) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 0, 4) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 0, 12) == false);

        REQUIRE(plane_3d_vec.same_pop(0, 4, 1) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 4, 5) == true);
        REQUIRE(plane_3d_vec.same_pop(0, 4, 6) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 4, 8) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 4, 16) == false);

        REQUIRE(plane_3d_vec.same_pop(0, 12, 1) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 12, 5) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 12, 6) == false);
        REQUIRE(plane_3d_vec.same_pop(2, 12, 13) == true);
        REQUIRE(plane_3d_vec.same_pop(2, 12, 15) == false);
    }

    SECTION("same_indiv with dif size pop")
    {
        genepop_input_c<1> gen_pop_input;
        //3 pop, 2-1-3 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}, {{3}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(gen_pop_input);
        //3 locus, 3 pop, 2-1-3 indiv
        //{{{{1}, {1}}, {{3}}, {{1}, {2}, {3}}},
        // {{{1}, {2}}, {{1}}, {{2}, {2}, {2}}},
        // {{{1}, {3}}, {{3}}, {{3}, {2}, {2}}}};
        //Same locus && same indiv ?
        REQUIRE(plane_3d_vec.same_indiv(0, 1) == false);
        REQUIRE(plane_3d_vec.same_indiv(0, 2) == false);
        REQUIRE(plane_3d_vec.same_indiv(0, 4) == false);
        REQUIRE(plane_3d_vec.same_indiv(0, 12) == false);

        REQUIRE(plane_3d_vec.same_indiv(4, 1) == false);
        REQUIRE(plane_3d_vec.same_indiv(4, 5) == false);
        REQUIRE(plane_3d_vec.same_indiv(4, 6) == false);
        REQUIRE(plane_3d_vec.same_indiv(4, 16) == false);
    }

    SECTION("same_pop with dif size pop")
    {
        genepop_input_c<1> gen_pop_input;
        //3 pop, 2-1-3 indiv, 3 locus
        //3 pop, 2-1-3 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}, {{3}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(gen_pop_input);
        //3 locus, 3 pop, 2-1-3 indiv
        //{{{{1}, {1}}, {{3}}, {{1}, {2}, {3}}},
        // {{{1}, {2}}, {{1}}, {{2}, {2}, {2}}},
        // {{{1}, {3}}, {{3}}, {{3}, {2}, {2}}}};

        //Same locus && same pop ?
        REQUIRE(plane_3d_vec.same_pop(0, 0, 1) == true);
        REQUIRE(plane_3d_vec.same_pop(0, 0, 2) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 0, 4) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 0, 12) == false);

        REQUIRE(plane_3d_vec.same_pop(0, 2, 1) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 2, 3) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 2, 4) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 2, 12) == false);

        REQUIRE(plane_3d_vec.same_pop(0, 3, 1) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 3, 4) == true);
        REQUIRE(plane_3d_vec.same_pop(0, 3, 5) == true);
        REQUIRE(plane_3d_vec.same_pop(0, 3, 8) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 3, 16) == false);
    }
}

TEST_CASE("diploid_custom_vec_test")
{
    SECTION("constructor")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(gen_pop_input);
        //3 locus, 3 pop, 2 indiv
        //{{{{1, 2}, {1, 1}}, {{1, 3}, {2, 3}}, {{1, 1}, {2, 1}}},
        // {{{1, 1}, {1, 2}}, {{1, 3}, {2, 3}}, {{1, 2}, {2, 2}}},
        // {{{1, 1}, {1, 3}}, {{1, 3}, {2, 3}}, {{1, 3}, {2, 3}}}};
        REQUIRE(plane_3d_vec.pop_nbr() == 3);
        REQUIRE(plane_3d_vec.indiv_nbr_per_pop(0) == 2);
        REQUIRE(plane_3d_vec.locus_nbr() == 3);
        REQUIRE(plane_3d_vec.indiv_nbr() == 6);

        REQUIRE(plane_3d_vec.cumul_indiv_nbr_per_pop() == std::vector<int>{0, 2, 4});

        std::vector<int> result = {1, 2, 1, 1, 1, 3, 2, 3, 1, 1, 2, 1, 1, 1, 1, 2, 1, 3, 2, 3, 1, 2, 2, 2, 1, 1, 1, 3, 1, 3, 2, 3, 1, 3, 2, 3};
        REQUIRE(plane_3d_vec.get_plane_vec() == result);
    }

    SECTION("operator()(int pop, int indiv, int locus)")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(gen_pop_input);
        //3 locus, 3 pop, 2 indiv
        REQUIRE(plane_3d_vec(0, 0, 0, 0) == 1);
        REQUIRE(plane_3d_vec(1, 1, 0, 1) == 3);
        REQUIRE(plane_3d_vec(2, 2, 1, 0) == 2);
    }

    SECTION("index_begin_locus && index_end_locus")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(gen_pop_input);
        //3 locus, 3 pop, 2 indiv
        REQUIRE(plane_3d_vec.index_begin_locus(0) == 0);
        REQUIRE(plane_3d_vec.index_end_locus(0) == 12);
        REQUIRE(plane_3d_vec.index_begin_locus(2) == 24);
        REQUIRE(plane_3d_vec.index_end_locus(2) == 36);
    }

    SECTION("same_indiv with same size pop")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(gen_pop_input);
        //3 locus, 3 pop, 2 indiv
        //{{{{1, 2}, {1, 1}}, {{1, 3}, {2, 3}}, {{1, 1}, {2, 1}}},
        // {{{1, 1}, {1, 2}}, {{1, 3}, {2, 3}}, {{1, 2}, {2, 2}}},
        // {{{1, 1}, {1, 3}}, {{1, 3}, {2, 3}}, {{1, 3}, {2, 3}}}};
        //Same locus && same indiv ?
        REQUIRE(plane_3d_vec.same_indiv(0, 1) == true);
        REQUIRE(plane_3d_vec.same_indiv(0, 2) == false);
        REQUIRE(plane_3d_vec.same_indiv(0, 4) == false);
        REQUIRE(plane_3d_vec.same_indiv(0, 12) == false);

        REQUIRE(plane_3d_vec.same_indiv(4, 1) == false);
        REQUIRE(plane_3d_vec.same_indiv(4, 5) == true);
        REQUIRE(plane_3d_vec.same_indiv(4, 6) == false);
        REQUIRE(plane_3d_vec.same_indiv(4, 16) == false);
    }
    SECTION("same_pop with same size pop")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(gen_pop_input);
        //3 locus, 3 pop, 2 indiv
        //{{{{1, 2}, {1, 1}}, {{1, 3}, {2, 3}}, {{1, 1}, {2, 1}}},
        // {{{1, 1}, {1, 2}}, {{1, 3}, {2, 3}}, {{1, 2}, {2, 2}}},
        // {{{1, 1}, {1, 3}}, {{1, 3}, {2, 3}}, {{1, 3}, {2, 3}}}};

        //Same locus && same pop ?
        REQUIRE(plane_3d_vec.same_pop(0, 0, 1) == true);
        REQUIRE(plane_3d_vec.same_pop(0, 0, 2) == true);
        REQUIRE(plane_3d_vec.same_pop(0, 0, 4) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 0, 12) == false);

        REQUIRE(plane_3d_vec.same_pop(0, 4, 1) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 4, 5) == true);
        REQUIRE(plane_3d_vec.same_pop(0, 4, 6) == true);
        REQUIRE(plane_3d_vec.same_pop(0, 4, 8) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 4, 16) == false);

        REQUIRE(plane_3d_vec.same_pop(0, 12, 1) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 12, 5) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 12, 6) == false);
        REQUIRE(plane_3d_vec.same_pop(1, 12, 14) == true);
        REQUIRE(plane_3d_vec.same_pop(1, 12, 15) == true);
        REQUIRE(plane_3d_vec.same_pop(1, 12, 16) == false);
    }

    SECTION("same_indiv with dif size pop")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2-1-3 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}, {{2, 3}, {2, 3}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(gen_pop_input);
        //3 locus, 3 pop, 2-1-3 indiv
        //{{{{1, 2}, {1, 1}}, {{1, 3}}, {{1, 1}, {2, 1}, {2, 3}}},
        // {{{1, 1}, {1, 2}}, {{1, 3}}, {{1, 2}, {2, 2}, {2, 3}}},
        // {{{1, 1}, {1, 3}}, {{1, 3}}, {{1, 3}, {2, 3}, {2, 3}}}};
        //Same locus && same indiv ?
        REQUIRE(plane_3d_vec.same_indiv(0, 1) == true);
        REQUIRE(plane_3d_vec.same_indiv(0, 2) == false);
        REQUIRE(plane_3d_vec.same_indiv(0, 4) == false);
        REQUIRE(plane_3d_vec.same_indiv(0, 12) == false);

        REQUIRE(plane_3d_vec.same_indiv(4, 1) == false);
        REQUIRE(plane_3d_vec.same_indiv(4, 5) == true);
        REQUIRE(plane_3d_vec.same_indiv(4, 6) == false);
        REQUIRE(plane_3d_vec.same_indiv(4, 16) == false);
    }

    SECTION("same_pop with dif size pop")
    {
        genepop_input_c<2> gen_pop_input;
        //3 pop, 2-1-3 indiv, 3 locus
        gen_pop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}, {{2, 3}, {2, 3}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(gen_pop_input);
        //3 locus, 3 pop, 2-1-3 indiv
        //{{{{1, 2}, {1, 1}}, {{1, 3}}, {{1, 1}, {2, 1}, {2, 3}}},
        // {{{1, 1}, {1, 2}}, {{1, 3}}, {{1, 2}, {2, 2}, {2, 3}}},
        // {{{1, 1}, {1, 3}}, {{1, 3}}, {{1, 3}, {2, 3}, {2, 3}}}};

        //Same locus && same pop ?
        REQUIRE(plane_3d_vec.same_pop(0, 0, 1) == true);
        REQUIRE(plane_3d_vec.same_pop(0, 0, 2) == true);
        REQUIRE(plane_3d_vec.same_pop(0, 0, 4) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 0, 12) == false);

        REQUIRE(plane_3d_vec.same_pop(0, 4, 1) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 4, 5) == true);
        auto result = plane_3d_vec.same_pop(0, 4, 6);
        REQUIRE(result == false);
        REQUIRE(plane_3d_vec.same_pop(0, 4, 8) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 4, 16) == false);

        REQUIRE(plane_3d_vec.same_pop(0, 12, 1) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 12, 5) == false);
        REQUIRE(plane_3d_vec.same_pop(0, 12, 6) == false);
        REQUIRE(plane_3d_vec.same_pop(1, 12, 14) == true);
        REQUIRE(plane_3d_vec.same_pop(1, 12, 15) == true);
        REQUIRE(plane_3d_vec.same_pop(1, 12, 16) == false);
    }
}