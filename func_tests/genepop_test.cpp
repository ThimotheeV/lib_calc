#define CATCH_CONFIG_MAIN GenedemeTest
#include "catch.hpp"
#include <iostream>

#include "input.hpp"
#include "common_tools.hpp"
#include "calc_stat.hpp"

TEST_CASE("haplo_genepop_input_test")
{
    SECTION("read_test")
    {
        genepop_input_c<1> input("genotype_genepop_format.txt");
        //Verify ADH-4 , ADH-5 \n mtDNA have been well separate
        REQUIRE(input.Locus_name.size() == 10);
        REQUIRE(input.Locus_name[4] == "4");
        REQUIRE(input.Locus_name[5] == "5");
        //deme name
        REQUIRE(input.Pop_name.size() == 5);
        REQUIRE(input.Indiv_name.size() == 5);
        REQUIRE(input.Pop_name[0][0] == "-0.1");
        REQUIRE(input.Pop_name[0][1] == "0.1");
        REQUIRE(input.Pop_name[3][0] == "-0.1");
        REQUIRE(input.Pop_name[3][1] == "3.1");
        //Name of indiv
        REQUIRE(input.Indiv_name.size() == 5);
        REQUIRE(input.Indiv_name[0][0] == "-0.1 0.1 ");
        REQUIRE(input.Indiv_name[2][1] == "-0.1 2.1 ");
        REQUIRE(input.Indiv_name[3][2] == "-0.1 3.1 ");
        REQUIRE(input.Indiv_name[4][2] == "-0.1 4.1 ");
        //Number of deme
        REQUIRE(input.Genotype.size() == 5);
        //Pop size
        REQUIRE(input.Genotype[0].size() == 10);
        //First indiv
        REQUIRE(input.Genotype[0][0].at(0) == std::array{37});
        REQUIRE(input.Genotype[0][0].at(5) == std::array{39});
        //Last indiv
        REQUIRE(input.Genotype[4][9].at(0) == std::array{1});
        REQUIRE(input.Genotype[4][9].at(9) == std::array{46});
    }

    SECTION("Geo_dist_btw_deme")
    {
        genepop_input_c<1> input("genotype_genepop_format.txt");
        REQUIRE(input.Geo_dist_btw_deme[0][0] == Approx(0).margin(0.0001));
        REQUIRE(input.Geo_dist_btw_deme[0][1] == Approx(1).margin(0.0001));
        REQUIRE(input.Geo_dist_btw_deme[0][2] == Approx(2).margin(0.0001));
        REQUIRE(input.Geo_dist_btw_deme[0][3] == Approx(3).margin(0.0001));
        REQUIRE(input.Geo_dist_btw_deme[0][4] == Approx(4).margin(0.0001));
    }

    SECTION("regres_ar_without_missing_value")
    {
        genepop_input_c<2> input("test_ar.txt");
        data_plane_vec_c data_plane_vec(input);

        auto ar = ar_by_pair(data_plane_vec);
        auto regr_ar_loc = linear_regres_X_Y(ar);

        REQUIRE(Approx(regr_ar_loc.at(0)).margin(0.000001) == -0.0184381);
        REQUIRE(Approx(regr_ar_loc.at(1)).margin(0.000001) == 0.238358);
    }

    SECTION("regres_ar_with_missing_value")
    {
        genepop_input_c<2> input("test_ar_mv.txt");
        data_plane_vec_c data_plane_vec(input);

        auto ar = ar_by_pair(data_plane_vec);
        auto regr_ar_loc = linear_regres_X_Y(ar);

        REQUIRE(Approx(regr_ar_loc.at(0)).margin(0.000001) == -0.0130383);
        REQUIRE(Approx(regr_ar_loc.at(1)).margin(0.00001) == 0.23328);
    }
}

TEST_CASE("func_test")
{
    SECTION("stat_test")
    {
        genepop_input_c<2> input("Test_func.txt", 15);
        data_plane_vec_c data_plane_vec(input);
        result_c result(15);

        std::cout << "\n######Hobs calculation######" << std::endl;
        std::vector<double> Vec_value(data_plane_vec.nbr_locus());
        double Hobs_mean = 0;
        int loc_abs = 0;
        for (std::size_t chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
        {
            for (std::size_t loc = 0; loc < data_plane_vec.nbr_locus(chr); ++loc)
            {
                Vec_value[loc_abs] = calc_Hobs_per_chr_per_loc(data_plane_vec, chr, loc);
                Hobs_mean += Vec_value[loc];
                ++loc_abs;
            }
        }

        result.Hobs_mean = Hobs_mean / data_plane_vec.nbr_locus();
        result.Hobs_var = var(Vec_value, result.Hobs_mean);

        /*******************************************/

        std::cout << "\n######Hexp calculation######" << std::endl;
        Vec_value.clear();
        Vec_value.resize(data_plane_vec.nbr_locus());
        double Hexp_mean = 0;
        loc_abs = 0;
        for (std::size_t chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
        {
            for (std::size_t loc = 0; loc < data_plane_vec.nbr_locus(chr); ++loc)
            {
                Vec_value[loc_abs] = calc_Hnei_per_chr_per_loc(data_plane_vec, chr, loc);
                Hexp_mean += Vec_value[loc];
                ++loc_abs;
            }
        }

        result.Hexp_mean = Hexp_mean / data_plane_vec.nbr_locus();
        result.Hexp_var = var(Vec_value, result.Hexp_mean);

        /*******************************************/

        std::cout << "\n######Var calculation######" << std::endl;
        Vec_value.clear();
        Vec_value.resize(data_plane_vec.nbr_locus());
        double var_mean = 0;
        loc_abs = 0;
        for (std::size_t chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
        {
            for (std::size_t loc = 0; loc < data_plane_vec.nbr_locus(chr); ++loc)
            {
                Vec_value[loc_abs] = calc_Var_per_chr_per_loc(data_plane_vec, chr, loc);
                var_mean += Vec_value[loc];
                ++loc_abs;
            }
        }
        result.Var_mean = var_mean / data_plane_vec.nbr_locus();
        result.Var_var = var(Vec_value, result.Var_mean);

        /*******************************************/

        std::cout << "\n######Nb_allele calculation######" << std::endl;
        Vec_value.clear();
        Vec_value.resize(data_plane_vec.nbr_locus());
        double nb_allele_mean = 0;
        loc_abs = 0;
        for (std::size_t chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
        {
            for (std::size_t loc = 0; loc < data_plane_vec.nbr_locus(chr); ++loc)
            {
                Vec_value[loc_abs] = data_plane_vec.nbr_allele(chr, loc);
                nb_allele_mean += Vec_value[loc];
                ++loc_abs;
            }
        }

        result.Nb_allele_mean = nb_allele_mean / data_plane_vec.nbr_locus();
        result.Nb_allele_var = var(Vec_value, result.Nb_allele_mean);

        /*******************************************/

        std::cout << "\n######MGW calculation######" << std::endl;
        Vec_value.clear();
        Vec_value.resize(data_plane_vec.nbr_locus());
        double MGW_mean = 0;
        loc_abs = 0;
        for (std::size_t chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
        {
            for (std::size_t loc = 0; loc < data_plane_vec.nbr_locus(chr); ++loc)
            {
                Vec_value[loc_abs] = calc_MGW_per_chr_per_loc(data_plane_vec, chr, loc);
                MGW_mean += Vec_value[loc];
                ++loc_abs;
            }
        }
        result.MGW_mean = MGW_mean / data_plane_vec.nbr_locus();
        result.MGW_var = var(Vec_value, result.MGW_mean);

        /*******************************************/

        std::cout << "\n######F_stat calculation######" << std::endl;

        auto temp = Fstat_genepop(data_plane_vec, false);
        result.Fis = temp.at(0);
        result.Fst = temp.at(1);

        /*******************************************/

        std::cout << "\n######Qr calculation######" << std::endl;
        result.Qr = calc_qr_all_loc(data_plane_vec);

        /*******************************************/

        std::cout << "\n######Ar calculation######" << std::endl;
        auto Ar = ar_by_pair(data_plane_vec);
        result.Ar_reg = linear_regres_X_Y(Ar);

        /*******************************************/

        std::cout << "\n######Er calculation######" << std::endl;
        auto er = er_by_pair(data_plane_vec);
        result.Er_reg = linear_regres_X_Y(er);

        /*******************************************/
        // Value from IBDSim
        REQUIRE(result.Hobs_mean == Approx(0.5555556).margin(0.0000001));
        REQUIRE(result.Hexp_mean == Approx(0.7237660).margin(0.0000001));
        REQUIRE(result.Var_mean == Approx(3.158315).margin(0.000001));
        REQUIRE(result.MGW_mean == Approx(0.9297619).margin(0.0000001));
        REQUIRE(result.Nb_allele_mean == 6.1);
        REQUIRE(result.Qr[0] == Approx(0.37333).margin(0.00001));
        REQUIRE(result.Qr[1] == Approx(0.31281).margin(0.00001));
        REQUIRE(result.Qr[2] == Approx(0.29392).margin(0.00001));
        REQUIRE(result.Qr[3] == Approx(0.28254).margin(0.00001));
        REQUIRE(result.Qr[4] == Approx(0.27965).margin(0.00001));
        REQUIRE(result.Qr[5] == Approx(0.27442).margin(0.00001));
        REQUIRE(result.Qr[6] == Approx(0.27813).margin(0.00001));
        REQUIRE(result.Qr[7] == Approx(0.26475).margin(0.00001));
        REQUIRE(result.Qr[8] == Approx(0.26013).margin(0.00001));
        REQUIRE(result.Qr[9] == Approx(0.26194).margin(0.00001));
        REQUIRE(result.Qr[10] == Approx(0.24749).margin(0.00001));
        REQUIRE(result.Qr[11] == Approx(0.26760).margin(0.00001));
        REQUIRE(result.Qr[12] == Approx(0.26709).margin(0.00001));
        // Value from Genepop
        REQUIRE(result.Ar_reg.at(0) == Approx(0.0108752).margin(0.0000001));  //Ar_slope
        REQUIRE(result.Ar_reg.at(1) == Approx(0.231215).margin(0.0000001));   //Ar_intercept
        REQUIRE(result.Er_reg.at(0) == Approx(0.00857859).margin(0.0000001)); //Er_slope
        REQUIRE(result.Er_reg.at(1) == Approx(-0.0572093).margin(0.0000001)); //Er_intercept
    }
}