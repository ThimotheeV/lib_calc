#define CATCH_CONFIG_MAIN CalcStatDLTest
#include "catch.hpp"

#include "calc_stat_dl.hpp"

TEST_CASE("calc_phi_calc_stat_dl")
{
    SECTION("Simple case controle")
    {
        //2 allelles max
        REQUIRE(calc_phi_ij_xy({{{1, 1, 1, 1}, {1, 1, 1, 1}}}) == 1);
        REQUIRE(calc_phi_ij_xy({{{1, 2, 2, 1}, {1, 1, 1, 1}}}) == 0.5);
        REQUIRE(calc_phi_ij_xy({{{1, 2, 2, 1}, {1, 2, 1, 2}}}) == 0.25);
        REQUIRE(calc_phi_ij_xy({{{1, 1, 2, 2}, {1, 1, 1, 1}}}) == 0);
        REQUIRE(calc_phi_ij_xy({{{1, 1, 2, 2}, {1, 2, 1, 2}}}) == 0);

        REQUIRE(calc_phi_ij_xy({{{1, 1, 1, 1}, {1, 2, 2, 1}}}) == 0.5);
        REQUIRE(calc_phi_ij_xy({{{1, 2, 1, 2}, {1, 2, 2, 1}}}) == 0.25);
        REQUIRE(calc_phi_ij_xy({{{1, 1, 1, 1}, {1, 1, 2, 2}}}) == 0);
        REQUIRE(calc_phi_ij_xy({{{1, 2, 1, 2}, {1, 1, 2, 2}}}) == 0);

        //3 alleles max
        REQUIRE(calc_phi_ij_xy({{{1, 1, 1, 1}, {1, 3, 3, 2}}}) == 0.25);
        REQUIRE(calc_phi_ij_xy({{{1, 2, 3, 2}, {1, 2, 1, 2}}}) == 0.125);
        REQUIRE(calc_phi_ij_xy({{{1, 3, 3, 2}, {1, 1, 1, 2}}}) == 0.125);
        REQUIRE(calc_phi_ij_xy({{{1, 2, 3, 2}, {3, 1, 3, 2}}}) == 0.0625);
        REQUIRE(calc_phi_ij_xy({{{1, 2, 3, 2}, {1, 3, 2, 2}}}) == 0);

        REQUIRE(calc_phi_ij_xy({{{1, 3, 3, 2}, {1, 1, 1, 1}}}) == 0.25);
        REQUIRE(calc_phi_ij_xy({{{1, 2, 1, 2}, {1, 2, 3, 2}}}) == 0.125);
        REQUIRE(calc_phi_ij_xy({{{1, 1, 1, 2}, {1, 3, 3, 2}}}) == 0.125);
        REQUIRE(calc_phi_ij_xy({{{3, 1, 3, 2}, {1, 2, 3, 2}}}) == 0.0625);
        REQUIRE(calc_phi_ij_xy({{{1, 3, 2, 2}, {1, 2, 3, 2}}}) == 0);

        //4 alleles max
        REQUIRE(calc_phi_ij_xy({{{1, 2, 3, 4}, {1, 1, 1, 1}}}) == 0);
        REQUIRE(calc_phi_ij_xy({{{1, 2, 1, 2}, {1, 2, 3, 4}}}) == 0);
    }

    SECTION("phi 2 diploid pop without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Dist_btw_deme = {{0, 1},
                                       {1, 0}};
        //2 deme, 4 indiv, 3 locus
        genepop_input.Genotype = {
            {{{1, 2}, {2, 2}, {1, 1}}, {{1, 2}, {1, 2}, {1, 3}}, {{1, 1}, {1, 1}, {1, 1}}, {{1, 1}, {2, 2}, {3, 3}}},
            {{{1, 3}, {2, 2}, {1, 3}}, {{2, 3}, {4, 3}, {2, 3}}, {{1, 3}, {4, 3}, {1, 3}}, {{1, 3}, {1, 1}, {1, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        //Loc pair 0-1
        //value come from gmf.ods (golden master file)
        std::array<double, 16> expect = {0.25, 0, 0, 0, 0.125, 0, 0, 0.125, 0, 0, 0, 0.5, 0.5, 0, 0, 0};
        auto expect_itr = expect.begin();
        int debug = 0;

        for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv_per_deme(0); ++indiv)
        {
            int locus_i_indiv1_gene1 = data_plane_vec(0, 0, indiv, 0);
            int locus_i_indiv1_gene2 = data_plane_vec(0, 0, indiv, 1);

            int locus_j_indiv1_gene1 = data_plane_vec(1, 0, indiv, 0);
            int locus_j_indiv1_gene2 = data_plane_vec(1, 0, indiv, 1);

            for (int indiv_other_deme = 0; indiv_other_deme < data_plane_vec.nbr_of_indiv_per_deme(1); ++indiv_other_deme)
            {
                int locus_i_indiv2_gene1 = data_plane_vec(0, 1, indiv_other_deme, 0);
                int locus_i_indiv2_gene2 = data_plane_vec(0, 1, indiv_other_deme, 1);

                int locus_j_indiv2_gene1 = data_plane_vec(1, 1, indiv_other_deme, 0);
                int locus_j_indiv2_gene2 = data_plane_vec(1, 1, indiv_other_deme, 1);

                REQUIRE(calc_phi_ij_xy({{{locus_i_indiv1_gene1, locus_i_indiv1_gene2, locus_i_indiv2_gene1, locus_i_indiv2_gene2}, {locus_j_indiv1_gene1, locus_j_indiv1_gene2, locus_j_indiv2_gene1, locus_j_indiv2_gene2}}}) == *expect_itr);
                ++expect_itr;
            }
        }

        //Loc pair 1-2
        expect = {0.5, 0, 0, 0, 0.25, 0, 0, 0.25, 0, 0, 0, 0.5, 0.5, 0, 0, 0};
        expect_itr = expect.begin();
        debug = 0;

        for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv_per_deme(0); ++indiv)
        {
            int locus_i_indiv1_gene1 = data_plane_vec(1, 0, indiv, 0);
            int locus_i_indiv1_gene2 = data_plane_vec(1, 0, indiv, 1);

            int locus_j_indiv1_gene1 = data_plane_vec(2, 0, indiv, 0);
            int locus_j_indiv1_gene2 = data_plane_vec(2, 0, indiv, 1);

            for (int indiv_other_deme = 0; indiv_other_deme < data_plane_vec.nbr_of_indiv_per_deme(1); ++indiv_other_deme)
            {
                int locus_i_indiv2_gene1 = data_plane_vec(1, 1, indiv_other_deme, 0);
                int locus_i_indiv2_gene2 = data_plane_vec(1, 1, indiv_other_deme, 1);

                int locus_j_indiv2_gene1 = data_plane_vec(2, 1, indiv_other_deme, 0);
                int locus_j_indiv2_gene2 = data_plane_vec(2, 1, indiv_other_deme, 1);

                REQUIRE(calc_phi_ij_xy({{{locus_i_indiv1_gene1, locus_i_indiv1_gene2, locus_i_indiv2_gene1, locus_i_indiv2_gene2}, {locus_j_indiv1_gene1, locus_j_indiv1_gene2, locus_j_indiv2_gene1, locus_j_indiv2_gene2}}}) == *expect_itr);
                ++expect_itr;
            }
        }

        //Loc pair 0-2
        expect = {0.125, 0, 0.125, 0.125, 0.125, 0.0625, 0.125, 0.125, 0.25, 0, 0.25, 0.25, 0.25, 0, 0.25, 0.25};
        expect_itr = expect.begin();
        debug = 0;

        for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv_per_deme(0); ++indiv)
        {
            int locus_i_indiv1_gene1 = data_plane_vec(0, 0, indiv, 0);
            int locus_i_indiv1_gene2 = data_plane_vec(0, 0, indiv, 1);

            int locus_j_indiv1_gene1 = data_plane_vec(2, 0, indiv, 0);
            int locus_j_indiv1_gene2 = data_plane_vec(2, 0, indiv, 1);

            for (int indiv_other_deme = 0; indiv_other_deme < data_plane_vec.nbr_of_indiv_per_deme(1); ++indiv_other_deme)
            {
                int locus_i_indiv2_gene1 = data_plane_vec(0, 1, indiv_other_deme, 0);
                int locus_i_indiv2_gene2 = data_plane_vec(0, 1, indiv_other_deme, 1);

                int locus_j_indiv2_gene1 = data_plane_vec(2, 1, indiv_other_deme, 0);
                int locus_j_indiv2_gene2 = data_plane_vec(2, 1, indiv_other_deme, 1);

                REQUIRE(calc_phi_ij_xy({{{locus_i_indiv1_gene1, locus_i_indiv1_gene2, locus_i_indiv2_gene1, locus_i_indiv2_gene2}, {locus_j_indiv1_gene1, locus_j_indiv1_gene2, locus_j_indiv2_gene1, locus_j_indiv2_gene2}}}) == *expect_itr);
                ++expect_itr;
            }
        }
    }

    SECTION("calc_phi_ij 2 diploid pop without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Dist_btw_deme = {{0, 1},
                                       {1, 0}};
        //2 deme, 3 indiv, 3 locus
        genepop_input.Genotype = {
            {{{1, 2}, {2, 2}, {1, 1}}, {{1, 2}, {1, 2}, {1, 3}}, {{1, 1}, {1, 1}, {1, 1}}},
            {{{1, 3}, {2, 2}, {1, 3}}, {{2, 3}, {4, 3}, {2, 3}}, {{1, 3}, {4, 3}, {1, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        //Loc pair 0-1
        //value come from gmf.ods (golden master file)
        std::array<double, 3> expect = {0.066666, 0.1625, 0.091666};

        auto result = calc_phi_ij(data_plane_vec, 2);

        for (int index = 0; index < expect.size(); ++index)
        {
            REQUIRE(expect[index] == Approx(result[index]).margin(0.000001));
        }
    }
}

TEST_CASE("calc_eta diploid without missing data calc_stat_dl")
{
    SECTION("eta 2 diploid pop without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Dist_btw_deme = {{0, 1},
                                       {1, 0}};

        //2 deme, 4 indiv, 3 locus
        genepop_input.Genotype = {
            {{{1, 2}, {2, 2}, {1, 1}}, {{1, 2}, {1, 2}, {1, 3}}, {{1, 1}, {1, 1}, {1, 1}}, {{1, 1}, {2, 2}, {3, 3}}},
            {{{1, 3}, {2, 2}, {1, 3}}, {{2, 3}, {4, 3}, {2, 3}}, {{1, 3}, {4, 3}, {1, 3}}, {{1, 3}, {1, 1}, {1, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        //value come from gmf.ods (golden master file)
        std::array<std::array<double, 3>, 3> expect = {{{1, 1, -0.0505051}, {1, 2, -0.1135681}, {1, 1, 0.0600601}}};

        auto result = calc_eta(data_plane_vec);
        auto result_itr = result.begin();

        for (auto value : expect)
        {
            REQUIRE(result_itr->at(0) == value.at(0));
            REQUIRE(result_itr->at(1) == value.at(1));
            REQUIRE(result_itr->at(2) == Approx(value.at(2)).margin(0.0000001));
            ++result_itr;
        }
    }

    SECTION("eta 3 diploid pop without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Dist_btw_deme = {{0, 1, 1},
                                       {1, 0, 1},
                                       {1, 1, 0}};

        //3 deme, 4 indiv, 3 locus
        genepop_input.Genotype = {
            {{{1, 2}, {2, 2}, {1, 1}}, {{1, 2}, {1, 2}, {1, 3}}, {{1, 1}, {1, 1}, {1, 1}}, {{1, 1}, {2, 2}, {3, 3}}},
            {{{1, 3}, {2, 2}, {1, 3}}, {{2, 3}, {4, 3}, {2, 3}}, {{1, 3}, {4, 3}, {1, 3}}, {{1, 3}, {1, 1}, {1, 3}}},
            {{{2, 3}, {1, 2}, {2, 3}}, {{3, 4}, {3, 1}, {4, 3}}, {{1, 2}, {2, 5}, {5, 4}}, {{1, 1}, {4, 5}, {1, 1}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        //value come from gmf.ods (golden master file)
        std::array<std::array<double, 3>, 9> expect = {{{1, 1, -0.0495151}, {1, 1, -0.0891272}, {1, 1, 0.0354858}, {1, 2, -0.0947932}, {1, 2, 0.0410162}, {1, 2, 0.0583342}, {1, 1, 0.0484066}, {1, 1, -0.0750302}, {1, 1, 0.0282371}}};

        auto result = calc_eta(data_plane_vec);
        auto result_itr = result.begin();

        for (auto value : expect)
        {
            REQUIRE(result_itr->at(0) == value.at(0));
            REQUIRE(result_itr->at(1) == value.at(1));
            REQUIRE(result_itr->at(2) == Approx(value.at(2)).margin(0.0000001));
            ++result_itr;
        }
    }

    SECTION("eta_q1_v 3 diploid pop without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Dist_btw_deme = {{0, 1, 1},
                                       {1, 0, 1},
                                       {1, 1, 0}};

        //3 deme, 4 indiv, 3 locus
        genepop_input.Genotype = {
            {{{1, 2}, {2, 2}, {1, 1}}, {{1, 2}, {1, 2}, {1, 3}}, {{1, 1}, {1, 1}, {1, 1}}, {{1, 1}, {2, 2}, {3, 3}}},
            {{{1, 3}, {2, 2}, {1, 3}}, {{2, 3}, {4, 3}, {2, 3}}, {{1, 3}, {4, 3}, {1, 3}}, {{1, 3}, {1, 1}, {1, 3}}},
            {{{2, 3}, {1, 2}, {2, 3}}, {{3, 4}, {3, 1}, {4, 3}}, {{1, 2}, {2, 5}, {5, 4}}, {{1, 1}, {4, 5}, {1, 1}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        //value come from gmf.ods (golden master file)
        std::array<std::array<double, 3>, 9> expect = {{{1, 1, -0.0666666}, {1, 1, -0.1028571}, {1, 1, 0.0285714}, {1, 2, -0.1434482}, {1, 2, 0.0416666}, {1, 2, 0.0494208}, {1, 1, 0.0574712}, {1, 1, -0.0738095}, {1, 1, 0.0219987}}};

        auto result = calc_eta_q1_version(data_plane_vec);
        auto result_itr = result.begin();

        for (auto value : expect)
        {
            REQUIRE(result_itr->at(0) == value.at(0));
            REQUIRE(result_itr->at(1) == value.at(1));
            REQUIRE(result_itr->at(2) == Approx(value.at(2)).margin(0.0000001));
            ++result_itr;
        }
    }

    SECTION("eta 3 diploid 1 indiv_pop without missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Dist_btw_deme = {{0, 1, 1},
                                       {1, 0, 1},
                                       {1, 1, 0}};

        //3 deme, 1 indiv, 3 locus
        genepop_input.Genotype = {
            {{{1, 2}, {2, 2}, {1, 1}}},
            {{{1, 3}, {2, 2}, {1, 3}}},
            {{{2, 3}, {1, 2}, {2, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        //value come from gmf.ods (golden master file)
        std::array<std::array<double, 3>, 9> expect = {{{1, 1, 1}, {1, 1, 0.5}, {1, 1, 0.5}, {1, 2, 0.222222}, {1, 2, 0}, {1, 2, 0.111111}, {1, 1, 1.111111}, {1, 1, -0.888888}, {1, 1, -0.388888}}};

        auto result = calc_eta_1_indiv_deme_v(data_plane_vec);
        auto result_itr = result.begin();

        for (auto value : expect)
        {
            REQUIRE(result_itr->at(0) == value.at(0));
            REQUIRE(result_itr->at(1) == value.at(1));
            REQUIRE(result_itr->at(2) == Approx(value.at(2)).margin(0.0000001));
            ++result_itr;
        }
    }
}

TEST_CASE("calc_eta haploid without missing data calc_stat_dl")
{
    SECTION("eta 2 haploid pop without missing value")
    {
        genepop_input_c<1> genepop_input;
        genepop_input.Dist_btw_deme = {{0, 1},
                                       {1, 0}};

        //2 deme, 4 indiv, 3 locus
        genepop_input.Genotype = {
            {{{1}, {2}, {1}}, {{1}, {1}, {1}}, {{1}, {1}, {1}}, {{1}, {2}, {3}}},
            {{{1}, {2}, {1}}, {{2}, {4}, {2}}, {{1}, {4}, {1}}, {{1}, {1}, {1}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        //value come from gmf.ods (golden master file)
        std::array<std::array<double, 3>, 3> expect = {{{1, 1, 0.3333333}, {1, 2, 1.7142857}, {1, 1, 0.1904761}}};

        auto result = calc_eta(data_plane_vec);
        auto result_itr = result.begin();

        for (auto value : expect)
        {
            REQUIRE(result_itr->at(0) == value.at(0));
            REQUIRE(result_itr->at(1) == value.at(1));
            REQUIRE(result_itr->at(2) == Approx(value.at(2)).margin(0.0000001));
            ++result_itr;
        }
    }

    SECTION("eta 3 haploid pop without missing value")
    {
        genepop_input_c<1> genepop_input;
        genepop_input.Dist_btw_deme = {{0, 1, 1},
                                       {1, 0, 1},
                                       {1, 1, 0}};

        //3 deme, 4 indiv, 3 locus
        genepop_input.Genotype = {
            {{{1}, {2}, {1}}, {{1}, {1}, {1}}, {{1}, {1}, {1}}, {{1}, {2}, {3}}},
            {{{1}, {2}, {1}}, {{2}, {4}, {2}}, {{1}, {4}, {1}}, {{1}, {1}, {1}}},
            {{{2}, {1}, {2}}, {{3}, {3}, {4}}, {{1}, {2}, {5}}, {{1}, {4}, {1}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        //value come from gmf.ods (golden master file)
        std::array<std::array<double, 3>, 9> expect = {
            {{1, 1, 0.1904761}, {1, 1, 0.0846560}, {1, 1, 0.2962962}, {1, 2, 0.6428571}, {1, 2, 0.1428571}, {1, 2, 0.5714285}, {1, 1, 0.125}, {1, 1, -0.0833333}, {1, 1, 0.0833333}}};

        auto result = calc_eta(data_plane_vec);
        auto result_itr = result.begin();

        for (auto value : expect)
        {
            REQUIRE(result_itr->at(0) == value.at(0));
            REQUIRE(result_itr->at(1) == value.at(1));
            REQUIRE(result_itr->at(2) == Approx(value.at(2)).margin(0.0000001));
            ++result_itr;
        }
    }

    SECTION("eta 3 haploid 1 indiv_pop without missing value")
    {
        genepop_input_c<1> genepop_input;
        genepop_input.Dist_btw_deme = {{0, 1, 1},
                                       {1, 0, 1},
                                       {1, 1, 0}};

        //3 deme, 1 indiv, 3 locus
        genepop_input.Genotype = {
            {{{1}, {2}, {1}}},
            {{{1}, {2}, {1}}},
            {{{2}, {1}, {2}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        //value come from gmf.ods (golden master file)
        std::array<std::array<double, 3>, 9> expect = {{{1, 1, 2}, {1, 1, -0.25}, {1, 1, -0.25}, {1, 2, 2}, {1, 2, -0.25}, {1, 2, -0.25}, {1, 1, 2}, {1, 1, -0.25}, {1, 1, -0.25}}};

        auto result = calc_eta_1_indiv_deme_v(data_plane_vec);
        auto result_itr = result.begin();

        for (auto value : expect)
        {
            REQUIRE(result_itr->at(0) == value.at(0));
            REQUIRE(result_itr->at(1) == value.at(1));
            REQUIRE(result_itr->at(2) == Approx(value.at(2)).margin(0.0000001));
            ++result_itr;
        }
    }
}

TEST_CASE("calc_eta diploid with missing data calc_stat_dl")
{
    SECTION("eta 2 diploid pop with missing value")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Dist_btw_deme = {{0, 1},
                                       {1, 0}};

        //2 deme, 4 indiv, 3 locus
        genepop_input.Genotype = {
            {{{1, 0}, {2, 2}, {1, 1}}, {{1, 2}, {1, 2}, {1, 3}}, {{1, 0}, {1, 0}, {1, 1}}, {{1, 1}, {2, 2}, {3, 3}}},
            {{{1, 3}, {2, 2}, {1, 3}}, {{2, 3}, {0, 3}, {2, 3}}, {{1, 3}, {4, 3}, {1, 3}}, {{1, 3}, {1, 0}, {1, 3}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        //value come from gmf.ods (golden master file)
        std::array<std::array<double, 3>, 3> expect = {{{1, 1, -0.1287030}, {1, 2, -0.1617969}, {1, 1, 0.0006974}}};

        auto result = calc_eta(data_plane_vec);
        auto result_itr = result.begin();

        for (auto value : expect)
        {
            REQUIRE(result_itr->at(0) == value.at(0));
            REQUIRE(result_itr->at(1) == value.at(1));
            REQUIRE(result_itr->at(2) == Approx(value.at(2)).margin(0.0000001));
            ++result_itr;
        }
    }
}

TEST_CASE("calc_eta corner case calc_stat_dl")
{
    SECTION("eta without difference")
    {
        genepop_input_c<2> genepop_input;
        genepop_input.Dist_btw_deme = {{0, 1},
                                       {1, 0}};

        //2 deme, 1 indiv, 2 locus
        genepop_input.Genotype = {
            {{{7, 7}, {1, 1}}},
            {{{6, 6}, {1, 1}}}};

        data_plane_vec_c data_plane_vec(genepop_input);

        REQUIRE_THROWS(calc_eta(data_plane_vec));
    }
}