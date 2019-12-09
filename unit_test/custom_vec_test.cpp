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
        REQUIRE(plane_3d_vec.get_feature(0).Pop == 0);
        REQUIRE(plane_3d_vec.get_feature(3).Pop == 1);
        REQUIRE(plane_3d_vec.get_feature(4).Pop == 2);
    }

    SECTION("get_indiv(int gene_index)")
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
        REQUIRE(plane_3d_vec.get_indiv(0) == 0);
        REQUIRE(plane_3d_vec.get_indiv(2) == 2);
        REQUIRE(plane_3d_vec.get_indiv(3) == 3);
        REQUIRE(plane_3d_vec.get_indiv(7) == 1);
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

    SECTION("index_begin_locus() && index_end_locus(")
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

    SECTION("same_indiv(int gene_index1, int gene_index2) with same size pop")
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
    SECTION("pop_at_dist(int gene_index1, int gene_index2, int dist) with same size pop")
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

        //same pop ?
        REQUIRE(plane_3d_vec.pop_at_dist(0, 1, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(0, 2, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(0, 4, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(0, 12, 0) == true);

        REQUIRE(plane_3d_vec.pop_at_dist(4, 1, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(4, 5, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(4, 6, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(4, 8, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(4, 16, 0) == true);

        REQUIRE(plane_3d_vec.pop_at_dist(12, 1, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(12, 5, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(12, 6, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(12, 13, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(12, 15, 0) == false);
    }

    SECTION("same_indiv(int gene_index1, int gene_index2) with dif size pop")
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

    SECTION("pop_at_dist(int gene_index1, int gene_index2, int dist) with dif size pop")
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

        //same pop ?
        REQUIRE(plane_3d_vec.pop_at_dist(0, 1, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(0, 2, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(0, 4, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(0, 12, 0) == true);

        REQUIRE(plane_3d_vec.pop_at_dist(2, 1, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(2, 3, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(2, 4, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(2, 12, 0) == false);

        REQUIRE(plane_3d_vec.pop_at_dist(3, 1, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(3, 4, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(3, 5, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(3, 8, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(3, 16, 0) == true);
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
        REQUIRE(plane_3d_vec.get_feature(0).Pop == 0);
        REQUIRE(plane_3d_vec.get_feature(3).Pop == 1);
        REQUIRE(plane_3d_vec.get_feature(4).Pop == 2);
    }

    SECTION("get_indiv(int gene_index)")
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
        REQUIRE(plane_3d_vec.get_indiv(0) == 0);
        REQUIRE(plane_3d_vec.get_indiv(2) == 1);
        REQUIRE(plane_3d_vec.get_indiv(4) == 2);
        REQUIRE(plane_3d_vec.get_indiv(12) == 0);
        REQUIRE(plane_3d_vec.get_indiv(14) == 1);
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

    SECTION("index_begin_locus() && index_end_locus()")
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

    SECTION("same_indiv(int gene_index1, int gene_index2) with same size pop")
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
    SECTION("pop_at_dist(int gene_index1, int gene_index2, int dist) with same size pop")
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

        //same pop ?
        REQUIRE(plane_3d_vec.pop_at_dist(0, 1, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(0, 2, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(0, 4, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(0, 12, 0) == true);

        REQUIRE(plane_3d_vec.pop_at_dist(4, 1, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(4, 5, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(4, 6, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(4, 8, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(4, 16, 0) == true);

        REQUIRE(plane_3d_vec.pop_at_dist(12, 1, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(12, 5, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(12, 6, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(12, 14, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(12, 15, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(12, 16, 0) == false);
    }

    SECTION("same_indiv(int gene_index1, int gene_index2) with dif size pop")
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

    SECTION("pop_at_dist(int gene_index1, int gene_index2, int dist) with dif size pop")
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

        //Same pop ?
        REQUIRE(plane_3d_vec.pop_at_dist(0, 1, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(0, 2, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(0, 4, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(0, 12, 0) == true);

        REQUIRE(plane_3d_vec.pop_at_dist(4, 1, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(4, 5, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(4, 6, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(4, 8, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(4, 16, 0) == true);

        REQUIRE(plane_3d_vec.pop_at_dist(12, 1, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(12, 5, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(12, 6, 0) == false);
        REQUIRE(plane_3d_vec.pop_at_dist(12, 14, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(12, 15, 0) == true);
        REQUIRE(plane_3d_vec.pop_at_dist(12, 16, 0) == false);
    }
}