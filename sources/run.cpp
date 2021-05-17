#include "run.hpp"

result_c run(selector_input_c const &selector, data_plane_vec_c const &data_plane_vec)
{
    result_c result(selector.Nbr_class);
    // bool stat = false;
    // if (selector.Hobs)
    // {
    //     stat = true;
    //     std::vector<double> Vec_value(data_plane_vec.base_nbr_locus());
    //     double Hobs_mean = 0;
    //     for (std::size_t loc = 0; loc < Vec_value.size(); ++loc)
    //     {
    //         Vec_value[loc] = calc_Hobs_per_loc(data_plane_vec, loc);
    //         Hobs_mean += Vec_value[loc];
    //     }

    //     result.Hobs_mean = Hobs_mean / data_plane_vec.base_nbr_locus();
    //     result.Hobs_var = var(Vec_value, result.Hobs_mean);
    // }

    // /*******************************************/

    // if (selector.Hexp)
    // {
    //     stat = true;
    //     std::vector<double> Vec_value(data_plane_vec.base_nbr_locus());
    //     double Hexp_mean = 0;
    //     for (std::size_t loc = 0; loc < Vec_value.size(); ++loc)
    //     {
    //         Vec_value[loc] = calc_Hnei_per_loc(data_plane_vec, loc);
    //         Hexp_mean += Vec_value[loc];
    //     }

    //     result.Hexp_mean = Hexp_mean / data_plane_vec.base_nbr_locus();
    //     result.Hexp_var = var(Vec_value, result.Hexp_mean);
    // }

    // /*******************************************/

    // if (selector.Var)
    // {
    //     stat = true;
    //     std::vector<double> Vec_value(data_plane_vec.base_nbr_locus());
    //     double var_mean = 0;
    //     for (std::size_t loc = 0; loc < Vec_value.size(); ++loc)
    //     {
    //         Vec_value[loc] = calc_Var_per_loc(data_plane_vec, loc);
    //         var_mean += Vec_value[loc];
    //     }
    //     result.Var_mean = var_mean / data_plane_vec.base_nbr_locus();
    //     result.Var_var = var(Vec_value, result.Var_mean);
    // }

    // /*******************************************/

    // if (selector.Nb_allele)
    // {
    //     stat = true;
    //     std::vector<double> Vec_value(data_plane_vec.base_nbr_locus());
    //     double nb_allele_mean = 0;
    //     for (std::size_t loc = 0; loc < Vec_value.size(); ++loc)
    //     {
    //         Vec_value[loc] = data_plane_vec.nbr_allele_per_loc(loc);
    //         nb_allele_mean += Vec_value[loc];
    //     }

    //     result.Nb_allele_mean = nb_allele_mean / data_plane_vec.base_nbr_locus();
    //     result.Nb_allele_var = var(Vec_value, result.Nb_allele_mean);
    // }

    // /*******************************************/

    // if (selector.MGW)
    // {
    //     stat = true;
    //     std::vector<double> Vec_value(data_plane_vec.base_nbr_locus());
    //     double MGW_mean = 0;
    //     for (std::size_t loc = 0; loc < Vec_value.size(); ++loc)
    //     {
    //         Vec_value[loc] = calc_MGW_per_loc(data_plane_vec, loc);
    //         MGW_mean += Vec_value[loc];
    //     }

    //     result.MGW_mean = MGW_mean / data_plane_vec.base_nbr_locus();
    //     result.MGW_var = var(Vec_value, result.MGW_mean);
    // }

    // /*******************************************/

    // if (selector.SFS)
    // {
    //     if (selector.min_gene_for_SFS < 0)
    //     {
    //         throw std::logic_error("( Can't calculate sfs if min_gene not set. I exit. )");
    //     }
    //     auto result = calc_SFS(data_plane_vec, selector.min_gene_for_SFS);
    //     std::map<int, double> output_map;
    //     auto result_itr = result.begin();
    //     //all state in result, need to take the n-1 state
    //     for (auto state_order = 1; state_order < result.size(); ++state_order)
    //     {
    //         for (auto count_freq : result_itr->second)
    //         {
    //             auto pair = output_map.emplace(count_freq.first, count_freq.second);
    //             if (!pair.second)
    //             {
    //                 pair.first->second += count_freq.second;
    //             }
    //         }
    //         ++result_itr;
    //     }

    //     output_sfs_stat_files(output_map);
    // }

    // /*******************************************/

    // if (selector.F_stat)
    // {
    //     stat = true;
    //     std::vector<double> Vec_Fis(data_plane_vec.base_nbr_locus());
    //     std::vector<double> Vec_Fst(data_plane_vec.base_nbr_locus());
    //     //Fis or Fst
    //     double Fis_mean = 0;
    //     double Fst_mean = 0;
    //     for (std::size_t loc = 0; loc < Vec_Fis.size(); ++loc)
    //     {
    //         //    //WARNING : Array <Fis, Fst> => <.at(0), .at(1)>
    //         std::array<double, 2> fract_Fis;
    //         std::array<double, 2> fract_Fst;
    //         if (selector.Missing_data)
    //         {
    //             auto temp = Fstat_by_loc_with_indic(data_plane_vec, loc);
    //             fract_Fis = temp.at(0);
    //             fract_Fst = temp.at(1);
    //         }
    //         else
    //         {
    //             auto temp = Fstat_by_loc_with_probid(data_plane_vec, loc);
    //             fract_Fis = temp.at(0);
    //             fract_Fst = temp.at(1);
    //         }
    //         Vec_Fis[loc] = fract_Fis.at(0) / fract_Fis.at(1);
    //         Fis_mean += Vec_Fis[loc];
    //         Vec_Fst[loc] = fract_Fst.at(0) / fract_Fst.at(1);
    //         Fst_mean += Vec_Fst[loc];
    //     }
    //     result.Fis_mean = Fis_mean / data_plane_vec.Nbr_of_locus();
    //     result.Fis_var = var(Vec_Fis, result.Fis_mean);
    //     result.Fst_mean = Fst_mean / data_plane_vec.Nbr_of_locus();
    //     result.Fst_var = var(Vec_Fst, result.Fst_mean);
    // }

    // /*******************************************/

    // if (selector.Qr)
    // {
    //     stat = true;
    //     result.Qr = calc_qr_all_loc(data_plane_vec);
    // }

    // /*******************************************/

    // if (selector.Ar)
    // {
    //     stat = true;
    //     auto Ar = ar_by_pair(data_plane_vec);
    //     result.Ar_reg = linear_regres_X_Y(Ar);
    // }

    // /*******************************************/

    // if (selector.Er)
    // {
    //     stat = true;
    //     auto er = er_by_pair(data_plane_vec);
    //     result.Er_reg = linear_regres_X_Y(er);
    // }

    // /*******************************************/

    // if (selector.Eta)
    // {
    //     if (data_plane_vec.get_Ploidy() == 2)
    //     {
    //         output_eta_stat_files(calc_eta(data_plane_vec));
    //     }
    //     else
    //     {
    //         output_eta_stat_files(calc_eta_q1_version(data_plane_vec));
    //     }
    // }

    // if (stat)
    // {
    //     output_stat_files(selector, result);
    // }

    return result;
}