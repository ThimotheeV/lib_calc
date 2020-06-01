#define CATCH_CONFIG_MAIN GenedemeTest
#include <catch2/catch.hpp>
#include <iostream>

#include "input.hpp"
#include "common_tools.hpp"
#include "calc_stat.hpp"

TEST_CASE("haplo_genepop_input_test")
{
    SECTION("read_test")
    {
        genepop_input_c<1> input("genotype_genedeme_format.txt");
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

    SECTION("Dist_btw_deme")
    {
        genepop_input_c<1> input("genotype_genedeme_format.txt");
        REQUIRE(input.Dist_btw_deme[0][0] == Approx(0).margin(0.0001));
        REQUIRE(input.Dist_btw_deme[0][1] == Approx(1).margin(0.0001));
        REQUIRE(input.Dist_btw_deme[0][2] == Approx(2).margin(0.0001));
        REQUIRE(input.Dist_btw_deme[0][3] == Approx(3).margin(0.0001));
        REQUIRE(input.Dist_btw_deme[0][4] == Approx(4).margin(0.0001));
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

        std::vector<double> Vec_value(data_plane_vec.nbr_of_locus());
        double Hexp_moy = 0;
        for (int loc = 0; loc < Vec_value.size(); ++loc)
        {
            Vec_value[loc] = calc_Hnei_per_loc(data_plane_vec, loc);
            Hexp_moy += Vec_value[loc];
        }

        Hexp_moy /= data_plane_vec.nbr_of_locus();

        /*******************************************/

        double var_moy = 0;
        for (int loc = 0; loc < Vec_value.size(); ++loc)
        {
            Vec_value[loc] = calc_Var_per_loc(data_plane_vec, loc);
            var_moy += Vec_value[loc];
        }
        var_moy /= data_plane_vec.nbr_of_locus();

        /*******************************************/
        double nb_allele_moy = 0;
        for (int loc = 0; loc < Vec_value.size(); ++loc)
        {
            Vec_value[loc] = data_plane_vec.nbr_allele_per_loc(loc);
            nb_allele_moy += Vec_value[loc];
        }

        nb_allele_moy /= data_plane_vec.nbr_of_locus();

        /*******************************************/

        double MGW_moy = 0;
        for (int loc = 0; loc < Vec_value.size(); ++loc)
        {
            Vec_value[loc] = calc_MGW_per_loc(data_plane_vec, loc);
            MGW_moy += Vec_value[loc];
        }

        MGW_moy /= data_plane_vec.nbr_of_locus();

        /*******************************************/

        double Hobs_moy = 0;
        for (int loc = 0; loc < Vec_value.size(); ++loc)
        {
            Vec_value[loc] = calc_Hobs_per_loc(data_plane_vec, loc);
            Hobs_moy += Vec_value[loc];
        }

        Hobs_moy /= data_plane_vec.nbr_of_locus();

        /*******************************************/

        //Fis or Fst
        double F_moy = 0;
        for (int loc = 0; loc < Vec_value.size(); ++loc)
        {
            //    //WARNING : Array <Fis, Fst> => <.at(0), .at(1)>
            auto fract_Fst = Fstat_by_loc_with_indic(data_plane_vec, loc).at(1);
            Vec_value[loc] = fract_Fst.at(0) / fract_Fst.at(1);
            F_moy += Vec_value[loc];
        }

        F_moy /= data_plane_vec.nbr_of_locus();

        /*******************************************/

        auto Qr = calc_qr_all_loc(data_plane_vec);

        /*******************************************/

        auto Ar = ar_by_pair(data_plane_vec);
        auto Ar_regr = linear_regres_X_Y(Ar);

        /*******************************************/

        auto er = er_by_pair(data_plane_vec);
        auto er_regr = linear_regres_X_Y(er);

        /*******************************************/
        // Value from IBDSim
        REQUIRE(Hobs_moy == Approx(0.5555556).margin(0.0000001));
        REQUIRE(Hexp_moy == Approx(0.7237660).margin(0.0000001));
        REQUIRE(var_moy == Approx(3.158315).margin(0.000001));
        REQUIRE(MGW_moy == Approx(0.9297619).margin(0.0000001));
        REQUIRE(F_moy == Approx(0.2452756).margin(0.0000001)); //Le Fst de genedeme => moyenne des Fst des locus dans genedeme (et non pas le Fst moyen de Genedeme)
        REQUIRE(nb_allele_moy == 6.1);
        REQUIRE(Qr[0] == Approx(0.37333).margin(0.00001));
        REQUIRE(Qr[1] == Approx(0.31281).margin(0.00001));
        REQUIRE(Qr[2] == Approx(0.29392).margin(0.00001));
        REQUIRE(Qr[3] == Approx(0.28254).margin(0.00001));
        REQUIRE(Qr[4] == Approx(0.27965).margin(0.00001));
        REQUIRE(Qr[5] == Approx(0.27442).margin(0.00001));
        REQUIRE(Qr[6] == Approx(0.27813).margin(0.00001));
        REQUIRE(Qr[7] == Approx(0.26475).margin(0.00001));
        REQUIRE(Qr[8] == Approx(0.26013).margin(0.00001));
        REQUIRE(Qr[9] == Approx(0.26194).margin(0.00001));
        REQUIRE(Qr[10] == Approx(0.24749).margin(0.00001));
        REQUIRE(Qr[11] == Approx(0.26760).margin(0.00001));
        REQUIRE(Qr[12] == Approx(0.26709).margin(0.00001));
        REQUIRE(Ar_regr.at(0) == Approx(0.0108752).margin(0.0000001));  //Ar_slope
        REQUIRE(Ar_regr.at(1) == Approx(0.231215).margin(0.0000001));   //Ar_intercept
        REQUIRE(er_regr.at(0) == Approx(0.00857859).margin(0.0000001)); //Er_slope
        REQUIRE(er_regr.at(1) == Approx(-0.0572093).margin(0.0000001)); //Er_intercept
    }
}