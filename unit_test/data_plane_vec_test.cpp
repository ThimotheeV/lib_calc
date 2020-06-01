#define CATCH_CONFIG_MAIN DataPlaneVec
#include <catch2/catch.hpp>

#include "data_plane_vec.hpp"

TEST_CASE("haploid_data_plane_vec_test")
{
    SECTION("constructor")
    {
        genepop_input_c<1> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}, {{3}, {2}, {2}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
        //{{{{1}, {1}}, {{3}, {3}}, {{1}, {2}}},
        // {{{1}, {2}}, {{1}, {2}}, {{2}, {2}}},
        // {{{1}, {3}}, {{3}, {2}}, {{3}, {2}}}};
        REQUIRE(plane_3d_vec.nbr_of_deme() == 3);
        REQUIRE(plane_3d_vec.nbr_of_indiv_per_deme(0) == 2);
        REQUIRE(plane_3d_vec.nbr_of_locus() == 3);
        REQUIRE(plane_3d_vec.nbr_of_indiv() == 6);

        REQUIRE(plane_3d_vec.cumul_nbr_of_indiv_per_deme() == std::vector<int>{0, 2, 4});
        REQUIRE(plane_3d_vec.nbr_allele_per_loc(1) == 2);
        REQUIRE(plane_3d_vec.allele_state_per_loc(1) == std::vector<std::array<int, 2>>{{1, 2}, {2, 4}});

        std::vector<int> result = {1, 1, 3, 3, 1, 2, 1, 2, 1, 2, 2, 2, 1, 3, 3, 2, 3, 2};
        REQUIRE(plane_3d_vec.get_plane_vec() == result);
        REQUIRE(plane_3d_vec.get_feature(0).Pop == 0);
        REQUIRE(plane_3d_vec.get_feature(3).Pop == 1);
        REQUIRE(plane_3d_vec.get_feature(4).Pop == 2);
    }

    SECTION("constructor_with_missing_value")
    {
        genepop_input_c<1> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1}, {1}, {0}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}, {{3}, {2}, {2}}},
                                  {{{1}, {0}, {3}}, {{2}, {0}, {2}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
        //{{{{1}, {1}}, {{3}, {3}}, {{1}, {2}}},
        // {{{1}, {2}}, {{1}, {2}}, {{0}, {0}}},
        // {{{0}, {3}}, {{3}, {2}}, {{3}, {2}}}};
        REQUIRE(plane_3d_vec.nbr_allele_per_loc(1) == 2);
        REQUIRE(plane_3d_vec.allele_state_per_loc(1) == std::vector<std::array<int, 2>>{{1, 2}, {2, 2}});

        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc(0) == 6);
        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc(1) == 4);
        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc(2) == 5);

        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc_per_deme(0) == std::vector<int>{2, 2, 2});
        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc_per_deme(0, 0) == 2);

        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc_per_deme(1) == std::vector<int>{2, 2});
        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc_per_deme(2) == std::vector<int>{1, 2, 2});
        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc_per_deme(2, 2) == 2);

        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc(0) == 6);
        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc(1) == 4);
        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc(2) == 5);

        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc_per_deme(0) == std::vector<int>{2, 2, 2});
        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc_per_deme(0, 0) == 2);

        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc_per_deme(1) == std::vector<int>{2, 2});
        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc_per_deme(2) == std::vector<int>{1, 2, 2});
        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc_per_deme(2, 2) == 2);

        REQUIRE(plane_3d_vec.nomiss_nbr_of_deme_per_loc(0) == 3);
    }

    SECTION("get_indiv(int gene_index)")
    {
        genepop_input_c<1> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}, {{3}, {2}, {2}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
        //{{{{1}, {1}}, {{3}, {3}}, {{1}, {2}}},
        // {{{1}, {2}}, {{1}, {2}}, {{2}, {2}}},
        // {{{1}, {3}}, {{3}, {2}}, {{3}, {2}}}};
        REQUIRE(plane_3d_vec.get_indiv(0) == 0);
        REQUIRE(plane_3d_vec.get_indiv(2) == 2);
        REQUIRE(plane_3d_vec.get_indiv(3) == 3);
        REQUIRE(plane_3d_vec.get_indiv(7) == 1);
    }

    SECTION("operator()(int locus, int deme, int indiv, int gene)")
    {
        genepop_input_c<1> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}, {{3}, {2}, {2}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
        REQUIRE(plane_3d_vec(0, 0, 0, 0) == 1);
        REQUIRE(plane_3d_vec(1, 1, 0, 0) == 1);
        REQUIRE(plane_3d_vec(2, 2, 1, 0) == 2);
    }

    SECTION("operator()(int locus, int indiv, int gene)")
    {
        genepop_input_c<1> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}, {{3}, {2}, {2}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
        REQUIRE(plane_3d_vec(0, 0, 0) == 1);
        REQUIRE(plane_3d_vec(0, 2, 0) == 3);
        REQUIRE(plane_3d_vec(2, 4, 0) == 3);
    }

    SECTION("index_begin_locus() && index_end_locus(")
    {
        genepop_input_c<1> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}, {{3}, {2}, {2}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
        REQUIRE(plane_3d_vec.index_begin_locus(0) == 0);
        REQUIRE(plane_3d_vec.index_end_locus(0) == 6);
        REQUIRE(plane_3d_vec.index_begin_locus(2) == 12);
        REQUIRE(plane_3d_vec.index_end_locus(2) == 18);
    }

    SECTION("same_indiv(int gene_index1, int gene_index2) with same size deme")
    {
        genepop_input_c<1> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}, {{3}, {2}, {2}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
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

    SECTION("deme_at_dist(int gene_index1, int gene_index2, int dist) with same size deme")
    {
        genepop_input_c<1> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}, {{3}, {2}, {2}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
        //{{{{1}, {1}}, {{3}, {3}}, {{1}, {2}}},
        // {{{1}, {2}}, {{1}, {2}}, {{2}, {2}}},
        // {{{1}, {3}}, {{3}, {2}}, {{3}, {2}}}};

        //same deme ?
        REQUIRE(plane_3d_vec.same_deme(0, 1) == true);
        REQUIRE(plane_3d_vec.same_deme(0, 2) == false);
        REQUIRE(plane_3d_vec.same_deme(0, 4) == false);
        REQUIRE(plane_3d_vec.same_deme(0, 12) == true);

        REQUIRE(plane_3d_vec.same_deme(4, 1) == false);
        REQUIRE(plane_3d_vec.same_deme(4, 5) == true);
        REQUIRE(plane_3d_vec.same_deme(4, 6) == false);
        REQUIRE(plane_3d_vec.same_deme(4, 8) == false);
        REQUIRE(plane_3d_vec.same_deme(4, 16) == true);

        REQUIRE(plane_3d_vec.same_deme(12, 1) == true);
        REQUIRE(plane_3d_vec.same_deme(12, 5) == false);
        REQUIRE(plane_3d_vec.same_deme(12, 6) == true);
        REQUIRE(plane_3d_vec.same_deme(12, 13) == true);
        REQUIRE(plane_3d_vec.same_deme(12, 15) == false);
    }

    SECTION("same_indiv(int gene_index1, int gene_index2) with dif size deme")
    {
        genepop_input_c<1> genepop_input;
        //3 deme, 2-1-3 indiv, 3 locus
        genepop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}, {{3}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2-1-3 indiv
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

    SECTION("bool same_deme(int gene_index1, int gene_index2) with dif size deme")
    {
        genepop_input_c<1> genepop_input;
        //3 deme, 2-1-3 indiv, 3 locus
        genepop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{3}, {1}, {3}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}, {{3}, {2}, {2}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2-1-3 indiv
        //{{{{1}, {1}}, {{3}}, {{1}, {2}, {3}}},
        // {{{1}, {2}}, {{1}}, {{2}, {2}, {2}}},
        // {{{1}, {3}}, {{3}}, {{3}, {2}, {2}}}};

        //same deme ?
        REQUIRE(plane_3d_vec.same_deme(0, 1) == true);
        REQUIRE(plane_3d_vec.same_deme(0, 2) == false);
        REQUIRE(plane_3d_vec.same_deme(0, 4) == false);
        REQUIRE(plane_3d_vec.same_deme(0, 12) == true);

        REQUIRE(plane_3d_vec.same_deme(2, 1) == false);
        REQUIRE(plane_3d_vec.same_deme(2, 3) == false);
        REQUIRE(plane_3d_vec.same_deme(2, 4) == false);
        REQUIRE(plane_3d_vec.same_deme(2, 12) == false);

        REQUIRE(plane_3d_vec.same_deme(3, 1) == false);
        REQUIRE(plane_3d_vec.same_deme(3, 4) == true);
        REQUIRE(plane_3d_vec.same_deme(3, 5) == true);
        REQUIRE(plane_3d_vec.same_deme(3, 8) == false);
        REQUIRE(plane_3d_vec.same_deme(3, 16) == true);
    }

    SECTION("bool nomiss_data_indiv_per_loc(int indiv, int locus) with missing data")
    {
        genepop_input_c<1> genepop_input;
        //3 deme, 2-1-3 indiv, 3 locus
        genepop_input.Genotype = {{{{1}, {1}, {1}}, {{1}, {2}, {3}}},
                                  {{{0}, {1}, {3}}},
                                  {{{1}, {2}, {3}}, {{2}, {2}, {2}}, {{3}, {2}, {0}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2-1-3 indiv
        //{{{{1}, {1}}, {{0}}, {{1}, {2}, {3}}},
        // {{{1}, {2}}, {{1}}, {{2}, {2}, {2}}},
        // {{{1}, {3}}, {{3}}, {{3}, {2}, {0}}}};

        REQUIRE(plane_3d_vec.nomiss_data_indiv_per_loc(0, 0) == true);
        REQUIRE(plane_3d_vec.nomiss_data_indiv_per_loc(2, 0) == false);
        REQUIRE(plane_3d_vec.nomiss_data_indiv_per_loc(5, 2) == false);
    }
}

TEST_CASE("diploid_data_plane_vec_test")
{
    SECTION("constructor")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
        //{{{{1, 2}, {1, 1}}, {{1, 3}, {2, 3}}, {{1, 1}, {2, 1}}},
        // {{{1, 1}, {1, 2}}, {{1, 3}, {2, 3}}, {{1, 2}, {2, 2}}},
        // {{{1, 1}, {1, 3}}, {{1, 3}, {2, 3}}, {{1, 3}, {2, 3}}}};
        REQUIRE(plane_3d_vec.nbr_of_deme() == 3);
        REQUIRE(plane_3d_vec.nbr_of_indiv_per_deme(0) == 2);
        REQUIRE(plane_3d_vec.nbr_of_locus() == 3);
        REQUIRE(plane_3d_vec.nbr_of_indiv() == 6);

        REQUIRE(plane_3d_vec.cumul_nbr_of_indiv_per_deme() == std::vector<int>{0, 2, 4});
        REQUIRE(plane_3d_vec.nbr_allele_per_loc(1) == 3);
        REQUIRE(plane_3d_vec.allele_state_per_loc(1) == std::vector<std::array<int, 2>>{{1, 5}, {2, 5}, {3, 2}});

        std::vector<int> result = {1, 2, 1, 1, 1, 3, 2, 3, 1, 1, 2, 1, 1, 1, 1, 2, 1, 3, 2, 3, 1, 2, 2, 2, 1, 1, 1, 3, 1, 3, 2, 3, 1, 3, 2, 3};
        REQUIRE(plane_3d_vec.get_plane_vec() == result);
        REQUIRE(plane_3d_vec.get_feature(0).Pop == 0);
        REQUIRE(plane_3d_vec.get_feature(3).Pop == 1);
        REQUIRE(plane_3d_vec.get_feature(4).Pop == 2);
    }

    SECTION("constructor_with_missing_value")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {0, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {0, 0}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
        //{{{{1, 2}, {1, 1}}, {{1, 3}, {2, 3}}, {{1, 1}, {2, 1}}},
        // {{{1, 1}, {1, 2}}, {{1, 3}, {0, 3}}, {{1, 2}, {2, 2}}},
        // {{{1, 1}, {1, 3}}, {{1, 3}, {2, 3}}, {{0, 0}, {2, 3}}}};

        REQUIRE(plane_3d_vec.nbr_allele_per_loc(1) == 3);
        REQUIRE(plane_3d_vec.allele_state_per_loc(1) == std::vector<std::array<int, 2>>{{1, 5}, {2, 4}, {3, 2}});

        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc(0) == 12);
        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc(1) == 11);
        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc(2) == 10);

        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc_per_deme(0) == std::vector<int>{4, 4, 4});
        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc_per_deme(0, 0) == 4);

        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc_per_deme(1) == std::vector<int>{4, 3, 4});
        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc_per_deme(2) == std::vector<int>{4, 4, 2});
        REQUIRE(plane_3d_vec.nomiss_nbr_of_gene_per_loc_per_deme(2, 2) == 2);

        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc(0) == 6);
        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc(1) == 5);
        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc(2) == 5);

        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc_per_deme(0) == std::vector<int>{2, 2, 2});
        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc_per_deme(0, 0) == 2);

        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc_per_deme(1) == std::vector<int>{2, 1, 2});
        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc_per_deme(2) == std::vector<int>{2, 2, 1});
        REQUIRE(plane_3d_vec.nomiss_nbr_of_indiv_per_loc_per_deme(2, 2) == 1);

        REQUIRE(plane_3d_vec.nomiss_nbr_of_deme_per_loc(0) == 3);
    }

    SECTION("get_indiv(int gene_index)")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
        //{{{{1, 2}, {1, 1}}, {{1, 3}, {2, 3}}, {{1, 1}, {2, 1}}},
        // {{{1, 1}, {1, 2}}, {{1, 3}, {2, 3}}, {{1, 2}, {2, 2}}},
        // {{{1, 1}, {1, 3}}, {{1, 3}, {2, 3}}, {{1, 3}, {2, 3}}}};
        REQUIRE(plane_3d_vec.get_indiv(0) == 0);
        REQUIRE(plane_3d_vec.get_indiv(2) == 1);
        REQUIRE(plane_3d_vec.get_indiv(4) == 2);
        REQUIRE(plane_3d_vec.get_indiv(12) == 0);
        REQUIRE(plane_3d_vec.get_indiv(14) == 1);
    }

    SECTION("operator()(int locus, int deme, int indiv, int gene)")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
        REQUIRE(plane_3d_vec(0, 0, 0, 0) == 1);
        REQUIRE(plane_3d_vec(1, 1, 0, 1) == 3);
        REQUIRE(plane_3d_vec(2, 2, 1, 0) == 2);
    }

    SECTION("operator()(int locus, int deme, int indiv, int gene)")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
        REQUIRE(plane_3d_vec(0, 0, 0) == 1);
        REQUIRE(plane_3d_vec(1, 2, 1) == 3);
        REQUIRE(plane_3d_vec(2, 5, 0) == 2);
    }

    SECTION("index_begin_locus() && index_end_locus()")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
        REQUIRE(plane_3d_vec.index_begin_locus(0) == 0);
        REQUIRE(plane_3d_vec.index_end_locus(0) == 12);
        REQUIRE(plane_3d_vec.index_begin_locus(2) == 24);
        REQUIRE(plane_3d_vec.index_end_locus(2) == 36);
    }

    SECTION("same_indiv(int gene_index1, int gene_index2) with same size deme")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
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
    SECTION("bool same_deme(int gene_index1, int gene_index2) with same size deme")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}, {{2, 3}, {2, 3}, {2, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2 indiv
        //{{{{1, 2}, {1, 1}}, {{1, 3}, {2, 3}}, {{1, 1}, {2, 1}}},
        // {{{1, 1}, {1, 2}}, {{1, 3}, {2, 3}}, {{1, 2}, {2, 2}}},
        // {{{1, 1}, {1, 3}}, {{1, 3}, {2, 3}}, {{1, 3}, {2, 3}}}};

        //same deme ?
        REQUIRE(plane_3d_vec.same_deme(0, 1) == true);
        REQUIRE(plane_3d_vec.same_deme(0, 2) == true);
        REQUIRE(plane_3d_vec.same_deme(0, 4) == false);
        REQUIRE(plane_3d_vec.same_deme(0, 12) == true);

        REQUIRE(plane_3d_vec.same_deme(4, 1) == false);
        REQUIRE(plane_3d_vec.same_deme(4, 5) == true);
        REQUIRE(plane_3d_vec.same_deme(4, 6) == true);
        REQUIRE(plane_3d_vec.same_deme(4, 8) == false);
        REQUIRE(plane_3d_vec.same_deme(4, 16) == true);

        REQUIRE(plane_3d_vec.same_deme(12, 1) == true);
        REQUIRE(plane_3d_vec.same_deme(12, 5) == false);
        REQUIRE(plane_3d_vec.same_deme(12, 6) == false);
        REQUIRE(plane_3d_vec.same_deme(12, 14) == true);
        REQUIRE(plane_3d_vec.same_deme(12, 15) == true);
        REQUIRE(plane_3d_vec.same_deme(12, 16) == false);
    }

    SECTION("same_indiv(int gene_index1, int gene_index2) with dif size deme")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2-1-3 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}, {{2, 3}, {2, 3}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2-1-3 indiv
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

    SECTION("bool same_deme(int gene_index1, int gene_index2) with dif size deme")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2-1-3 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}, {{2, 3}, {2, 3}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2-1-3 indiv
        //{{{{1, 2}, {1, 1}}, {{1, 3}}, {{1, 1}, {2, 1}, {2, 3}}},
        // {{{1, 1}, {1, 2}}, {{1, 3}}, {{1, 2}, {2, 2}, {2, 3}}},
        // {{{1, 1}, {1, 3}}, {{1, 3}}, {{1, 3}, {2, 3}, {2, 3}}}};

        //Same deme ?
        REQUIRE(plane_3d_vec.same_deme(0, 1) == true);
        REQUIRE(plane_3d_vec.same_deme(0, 2) == true);
        REQUIRE(plane_3d_vec.same_deme(0, 4) == false);
        REQUIRE(plane_3d_vec.same_deme(0, 12) == true);

        REQUIRE(plane_3d_vec.same_deme(4, 1) == false);
        REQUIRE(plane_3d_vec.same_deme(4, 5) == true);
        REQUIRE(plane_3d_vec.same_deme(4, 6) == false);
        REQUIRE(plane_3d_vec.same_deme(4, 8) == false);
        REQUIRE(plane_3d_vec.same_deme(4, 16) == true);

        REQUIRE(plane_3d_vec.same_deme(12, 1) == true);
        REQUIRE(plane_3d_vec.same_deme(12, 5) == false);
        REQUIRE(plane_3d_vec.same_deme(12, 6) == false);
        REQUIRE(plane_3d_vec.same_deme(12, 14) == true);
        REQUIRE(plane_3d_vec.same_deme(12, 15) == true);
        REQUIRE(plane_3d_vec.same_deme(12, 16) == false);
    }

    SECTION("bool nomiss_data_indiv_per_loc(int indiv, int locus) with missing data")
    {
        genepop_input_c<2> genepop_input;
        //3 deme, 2-1-3 indiv, 3 locus
        genepop_input.Genotype = {{{{0, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}, {{2, 3}, {2, 3}, {2, 0}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2-1-3 indiv
        //{{{{0, 2}, {1, 1}}, {{1, 3}}, {{1, 1}, {2, 1}, {2, 3}}},
        // {{{1, 1}, {1, 2}}, {{1, 3}}, {{1, 2}, {2, 2}, {2, 3}}},
        // {{{1, 1}, {1, 3}}, {{1, 3}}, {{1, 3}, {2, 3}, {2, 0}}}};

        REQUIRE(plane_3d_vec.nomiss_data_indiv_per_loc(0, 0) == false);
        REQUIRE(plane_3d_vec.nomiss_data_indiv_per_loc(2, 0) == true);
        REQUIRE(plane_3d_vec.nomiss_data_indiv_per_loc(5, 2) == false);
    }
}

TEST_CASE("Indiv_feat_test")
{
    SECTION("double dist_btw_deme(int gene_index1, int gene_index2)")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Nbr_dist_class = 3;
        genepop_input.Dist_btw_deme = {{0, 1, 2},
                                      {1, 0, 1},
                                      {2, 1, 0}};
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}, {{2, 3}, {2, 3}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2-1-3 indiv
        //{{{{1, 2}, {1, 1}}, {{1, 3}}, {{1, 1}, {2, 1}, {2, 3}}},
        // {{{1, 1}, {1, 2}}, {{1, 3}}, {{1, 2}, {2, 2}, {2, 3}}},
        // {{{1, 1}, {1, 3}}, {{1, 3}}, {{1, 3}, {2, 3}, {2, 3}}}};

        REQUIRE(plane_3d_vec.nbr_of_dist_class() == 3);

        REQUIRE(plane_3d_vec.dist_btw_deme(0, 1) == 0);
        REQUIRE(plane_3d_vec.dist_btw_deme(0, 2) == 0);
        REQUIRE(plane_3d_vec.dist_btw_deme(0, 4) == 1);
        REQUIRE(plane_3d_vec.dist_btw_deme(0, 7) == 2);
        REQUIRE(plane_3d_vec.dist_btw_deme(0, 12) == 0);

        REQUIRE(plane_3d_vec.dist_btw_deme(4, 1) == 1);
        REQUIRE(plane_3d_vec.dist_btw_deme(4, 5) == 0);
        REQUIRE(plane_3d_vec.dist_btw_deme(4, 6) == 1);
        REQUIRE(plane_3d_vec.dist_btw_deme(4, 8) == 1);
        REQUIRE(plane_3d_vec.dist_btw_deme(4, 16) == 0);

        REQUIRE(plane_3d_vec.dist_btw_deme(12, 1) == 0);
        REQUIRE(plane_3d_vec.dist_btw_deme(12, 5) == 1);
        REQUIRE(plane_3d_vec.dist_btw_deme(12, 6) == 2);
        REQUIRE(plane_3d_vec.dist_btw_deme(12, 14) == 0);
        REQUIRE(plane_3d_vec.dist_btw_deme(12, 15) == 0);
        REQUIRE(plane_3d_vec.dist_btw_deme(12, 16) == 1);
    }

    SECTION("int dist_class_btw_deme(int gene_index1, int gene_index2)")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Nbr_dist_class = 3;
        genepop_input.Dist_class_btw_deme = {{0, 1, 2},
                                            {1, 0, 1},
                                            {2, 1, 0}};
        //3 deme, 2 indiv, 3 locus
        genepop_input.Genotype = {{{{1, 2}, {1, 1}, {1, 1}}, {{1, 1}, {1, 2}, {1, 3}}},
                                  {{{1, 3}, {1, 3}, {1, 3}}},
                                  {{{1, 1}, {1, 2}, {1, 3}}, {{2, 1}, {2, 2}, {2, 3}}, {{2, 3}, {2, 3}, {2, 3}}}};

        data_plane_vec_c plane_3d_vec(genepop_input);
        //3 locus, 3 deme, 2-1-3 indiv
        //{{{{1, 2}, {1, 1}}, {{1, 3}}, {{1, 1}, {2, 1}, {2, 3}}},
        // {{{1, 1}, {1, 2}}, {{1, 3}}, {{1, 2}, {2, 2}, {2, 3}}},
        // {{{1, 1}, {1, 3}}, {{1, 3}}, {{1, 3}, {2, 3}, {2, 3}}}};
        REQUIRE(plane_3d_vec.nbr_of_dist_class() == 3);

        REQUIRE(plane_3d_vec.dist_class_btw_deme(0, 1) == 0);
        REQUIRE(plane_3d_vec.dist_class_btw_deme(0, 2) == 0);
        REQUIRE(plane_3d_vec.dist_class_btw_deme(0, 4) == 1);
        REQUIRE(plane_3d_vec.dist_class_btw_deme(0, 7) == 2);
        REQUIRE(plane_3d_vec.dist_class_btw_deme(0, 12) == 0);

        REQUIRE(plane_3d_vec.dist_class_btw_deme(4, 1) == 1);
        REQUIRE(plane_3d_vec.dist_class_btw_deme(4, 5) == 0);
        REQUIRE(plane_3d_vec.dist_class_btw_deme(4, 6) == 1);
        REQUIRE(plane_3d_vec.dist_class_btw_deme(4, 8) == 1);
        REQUIRE(plane_3d_vec.dist_class_btw_deme(4, 16) == 0);

        REQUIRE(plane_3d_vec.dist_class_btw_deme(12, 1) == 0);
        REQUIRE(plane_3d_vec.dist_class_btw_deme(12, 5) == 1);
        REQUIRE(plane_3d_vec.dist_class_btw_deme(12, 6) == 2);
        REQUIRE(plane_3d_vec.dist_class_btw_deme(12, 14) == 0);
        REQUIRE(plane_3d_vec.dist_class_btw_deme(12, 15) == 0);
        REQUIRE(plane_3d_vec.dist_class_btw_deme(12, 16) == 1);
    }
}