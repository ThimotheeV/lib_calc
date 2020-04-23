#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>

#include "calc_stat.hpp"
#include "output_file.hpp"

//TODO : Test dif opti float/double jeux de données de lézard avec 13 000 SNP
int main(int argc, char *argv[])
{
    std::string file_in_name(argv[1]);
    std::string file_out_name(argv[2]);

    genepop_input_c<2> input(file_in_name, 3);
    data_plane_vec_c data_plane_vec(input);

    std::vector<double> Vec_value(data_plane_vec.nbr_of_locus());
    // double Hexp_moy = 0;
    // for (int loc = 0; loc < Vec_value.size(); ++loc)
    // {
    //     Vec_value[loc] = calc_Hnei_per_loc(data_plane_vec, loc);
    //     Hexp_moy += Vec_value[loc];
    // }

    // Hexp_moy /= data_plane_vec.nbr_of_locus();
    // // double Hexp_V = var(Vec_value, Hexp_moy);

    // /*******************************************/

    // double var_moy = 0;
    // for (int loc = 0; loc < Vec_value.size(); ++loc)
    // {
    //     Vec_value[loc] = calc_Var_per_loc(data_plane_vec, loc);
    //     var_moy += Vec_value[loc];
    // }
    // var_moy /= data_plane_vec.nbr_of_locus();
    // // double var_V = var(Vec_value, var_moy);

    // /*******************************************/
    // double nb_allele_moy = 0;
    // for (int loc = 0; loc < Vec_value.size(); ++loc)
    // {
    //     Vec_value[loc] = data_plane_vec.nbr_allele_per_loc(loc);
    //     nb_allele_moy += Vec_value[loc];
    // }

    // nb_allele_moy /= data_plane_vec.nbr_of_locus();
    // // double nb_allele_V = var(Vec_value, nb_allele_moy);

    // /*******************************************/

    // double MGW_moy = 0;
    // for (int loc = 0; loc < Vec_value.size(); ++loc)
    // {
    //     Vec_value[loc] = calc_MGW_per_loc(data_plane_vec, loc);
    //     MGW_moy += Vec_value[loc];
    // }

    // MGW_moy /= data_plane_vec.nbr_of_locus();
    // // double MGW_V = var(Vec_value, MGW_moy);

    // /*******************************************/

    // double Hobs_moy = 0;
    // for (int loc = 0; loc < Vec_value.size(); ++loc)
    // {
    //     Vec_value[loc] = calc_Hobs_per_loc(data_plane_vec, loc);
    //     Hobs_moy += Vec_value[loc];
    // }

    // Hobs_moy /= data_plane_vec.nbr_of_locus();
    // // double Hobs_V = var(Vec_value, Hobs_moy);

    /*******************************************/

    //Fis or Fst
    double F_moy = 0;
    for (int loc = 0; loc < Vec_value.size(); ++loc)
    {
        //    //WARNING : Array <Fis, Fst> => <.at(0), .at(1)>
        auto fract_Fst = Fstat_by_loc_with_indic(data_plane_vec, loc).at(0);
        std::cout << "loc " << fract_Fst.at(0) <<"/"<< fract_Fst.at(1) << "\n";
        Vec_value[loc] = fract_Fst.at(0) / fract_Fst.at(1);
        F_moy += Vec_value[loc];
    }

    F_moy /= data_plane_vec.nbr_of_locus();
    // double F_V = var(Vec_value, F_moy);

    // /*******************************************/

    // auto Qr = calc_qr_all_loc(data_plane_vec);

    // /*******************************************/

    // auto Ar = ar_by_pair(data_plane_vec);
    // for (auto value : Ar)
    // {
    //     std::cout << value.at(0) << " " << value.at(1) << "\n";
    // }
    // auto Ar_regr = linear_regres_X_Y(Ar);

    // /*******************************************/

    // auto er = er_by_pair(data_plane_vec);
    // auto er_regr = linear_regres_X_Y(er);

    // /*******************************************/

    // std::cout << "Hobs_moy : " << Hobs_moy << "\n";
    // std::cout << "Hexp_moy : " << Hexp_moy << "\n";
    // std::cout << "var_moy : " << var_moy << "\n";
    // std::cout << "MGW_moy : " << MGW_moy << "\n";
    // std::cout << "F_moy : " << F_moy << "\n";
    // std::cout << "nb_allele_moy : " << nb_allele_moy << "\n";
    // std::cout << "Qr[0] : " << Qr[0] << "\n";
    // std::cout << "Qr[1] : " << Qr[1] << "\n";
    // std::cout << "Qr[2] : " << Qr[2] << "\n";
    // std::cout << "Qr[3] : " << Qr[3] << "\n";
    // std::cout << "Qr[4] : " << Qr[4] << "\n";
    // std::cout << "Qr[5] : " << Qr[5] << "\n";
    // std::cout << "Qr[6] : " << Qr[6] << "\n";
    // std::cout << "Qr[7] : " << Qr[7] << "\n";
    // std::cout << "Qr[8] : " << Qr[8] << "\n";
    // std::cout << "Qr[9] : " << Qr[9] << "\n";
    // std::cout << "Qr10 : " << Qr[10] << "\n";
    // std::cout << "Qr[11] : " << Qr[11] << "\n";
    // std::cout << "Qr[12] : " << Qr[12] << "\n";
    // std::cout << "Ar_slope : " << Ar_regr.at(0) << "\n";
    // std::cout << "Ar_intercept : " << Ar_regr.at(1) << "\n";
    // std::cout << "Er_slope : " << er_regr.at(0) << "\n";
    // std::cout << "Er_intercept : " << er_regr.at(1) << "\n";

    // std::ofstream ar_GRA("./Result/ar.GRA");
    // if (ar_GRA.is_open())
    // {
    //     for (auto ar : Ar)
    //     {
    //         ar_GRA << std::setprecision(7) << ar.at(0) << " " << ar.at(1) << std::endl;
    //     }
    // }
    // else
    // {
    //     throw std::invalid_argument("Unable to open ar.GRA");
    // }

    // std::ofstream lezard(file_out_name);
    // if (lezard.is_open())
    // {
    //     lezard << "***********Genet Estim*******"
    //            << "\n";
    //     for (auto val : Ar)
    //     {
    //         lezard << std::setprecision(7) << val.at(0) << "\n";
    //     }
    //     lezard << "Hobs_moy\tHexp_moy\tVar_moy\tMGW_moy\tF_moy\tnb_allele_moy\tQr0\tQr1\tQr2\tQr3\tQr4\tQr5\tQr6\tQr7\tQr8\tQr9\tQr10\tQr11\tQr12\tar_slope\tar_intercept\ter_slope\ter_intercept\n";
    //     lezard << std::setprecision(7) << Hobs_moy << "\n";
    //     lezard << Hexp_moy << "\n";
    //     lezard << var_moy << "\n";
    //     lezard << MGW_moy << "\n";
    //     lezard << F_moy << "\n";
    //     lezard << nb_allele_moy << "\n";
    //     lezard << Qr[0] << "\n";
    //     lezard << Qr[1] << "\n";
    //     lezard << Qr[2] << "\n";
    //     lezard << Qr[3] << "\n";
    //     lezard << Qr[4] << "\n";
    //     lezard << Qr[5] << "\n";
    //     lezard << Qr[6] << "\n";
    //     lezard << Qr[7] << "\n";
    //     lezard << Qr[8] << "\n";
    //     lezard << Qr[9] << "\n";
    //     lezard << Qr[10] << "\n";
    //     lezard << Qr[11] << "\n";
    //     lezard << Qr[12] << "\n";
    //     lezard << Ar_regr.at(0) << "\n";
    //     lezard << Ar_regr.at(1) << "\n";
    //     lezard << er_regr.at(0) << "\n";
    //     lezard << er_regr.at(1) << "\n";
    //     lezard.close();
    // }
    // else
    // {
    //     throw std::invalid_argument("Unable to open lezard.txt");
    // }
    // // Hexp_moy nb_allele_moy       Var_moy       MGW_moy      Hobs_moy         F_moy           Qr0           Qr1           Qr2             Qr3           Qr4           Qr5           Qr6           Qr7           Qr8           Qr9          Qr10          Qr11      ar_slope  ar_intercept        er_slope  er_intercept      Hexp_var nb_allele_var       Var_var       MGW_var      Hobs_var         F_var
}