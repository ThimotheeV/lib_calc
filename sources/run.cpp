#include "run.hpp"

result_c run(selector_input_c const &selector, data_plane_vec_c const &data_plane_vec)
{
    result_c result(selector.Nbr_geo_dist_class);
    bool stat = false;
    if (selector.Hobs)
    {
        std::cout << "######Hobs calculation######" << std::endl;
        stat = true;
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
    }

    /*******************************************/

    if (selector.Hexp)
    {
        std::cout << "######Hexp calculation######" << std::endl;
        stat = true;
        std::vector<double> Vec_value(data_plane_vec.nbr_locus());
        double Hexp_mean = 0;
        int loc_abs = 0;
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
    }

    /*******************************************/

    if (selector.Var)
    {
        std::cout << "######Var calculation######" << std::endl;
        stat = true;
        std::vector<double> Vec_value(data_plane_vec.nbr_locus());
        double var_mean = 0;
        int loc_abs = 0;
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
    }

    /*******************************************/

    if (selector.Nb_allele)
    {
        std::cout << "######Nb_allele calculation######" << std::endl;
        stat = true;
        std::vector<double> Vec_value(data_plane_vec.nbr_locus());
        double nb_allele_mean = 0;
        int loc_abs = 0;
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
    }

    /*******************************************/

    if (selector.MGW)
    {
        std::cout << "######MGW calculation######" << std::endl;
        stat = true;
        std::vector<double> Vec_value(data_plane_vec.nbr_locus());
        double MGW_mean = 0;
        int loc_abs = 0;
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
    }

    /*******************************************/

    if (selector.SFS)
    {
        std::cout << "######SFS calculation######" << std::endl;
        if (selector.min_gene_for_SFS < 0)
        {
            throw std::logic_error("( Can't calculate sfs if min_gene not set. I exit. )");
        }
        auto result = calc_SFS(data_plane_vec, selector.min_gene_for_SFS);
        std::map<int, double> output_map;
        auto result_itr = result.begin();
        //all state in result, need to take the n-1 state
        for (auto state_order = 1; state_order < result.size(); ++state_order)
        {
            for (auto count_freq : result_itr->second)
            {
                auto pair = output_map.emplace(count_freq.first, count_freq.second);
                if (!pair.second)
                {
                    pair.first->second += count_freq.second;
                }
            }
            ++result_itr;
        }

        output_sfs_stat_files(output_map);
    }

    /*******************************************/

    if (selector.F_stat)
    {
        std::cout << "######F_stat calculation######" << std::endl;
        stat = true;
        std::vector<double> Vec_Fis(data_plane_vec.nbr_locus());
        std::vector<double> Vec_Fst(data_plane_vec.nbr_locus());
        //Fis or Fst
        double Fis_mean = 0;
        double Fst_mean = 0;
        int loc_abs = 0;
        for (std::size_t chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
        {
            for (std::size_t loc = 0; loc < data_plane_vec.nbr_locus(chr); ++loc)
            {
                //WARNING : Array <Fis, Fst> => <.at(0), .at(1)>
                std::array<double, 2> fract_Fis;
                std::array<double, 2> fract_Fst;
                if (selector.Missing_data)
                {
                    auto temp = Fstat_per_chr_by_loc_with_indic(data_plane_vec, chr, loc);
                    fract_Fis = temp.at(0);
                    fract_Fst = temp.at(1);
                }
                else
                {
                    auto temp = Fstat_per_chr_by_loc_with_probid(data_plane_vec, chr, loc);
                    fract_Fis = temp.at(0);
                    fract_Fst = temp.at(1);
                }

                Vec_Fis[loc_abs] = fract_Fis.at(0) / fract_Fis.at(1);
                Fis_mean += Vec_Fis[loc_abs];
                Vec_Fst[loc_abs] = fract_Fst.at(0) / fract_Fst.at(1);
                Fst_mean += Vec_Fst[loc_abs];
                ++loc_abs;
            }
        }
        result.Fis_mean = Fis_mean / data_plane_vec.nbr_locus();
        result.Fis_var = var(Vec_Fis, result.Fis_mean);
        result.Fst_mean = Fst_mean / data_plane_vec.nbr_locus();
        result.Fst_var = var(Vec_Fst, result.Fst_mean);
    }

    /*******************************************/

    if (selector.Qr)
    {
        std::cout << "######Qr calculation######" << std::endl;
        stat = true;
        result.Qr = calc_qr_all_loc(data_plane_vec);
    }

    /*******************************************/

    if (selector.Ar)
    {
        std::cout << "######Ar calculation######" << std::endl;
        stat = true;
        auto Ar = ar_by_pair(data_plane_vec);
        result.Ar_reg = linear_regres_X_Y(Ar);
    }

    /*******************************************/

    if (selector.Er)
    {
        std::cout << "######Er calculation######" << std::endl;
        stat = true;
        auto er = er_by_pair(data_plane_vec);
        result.Er_reg = linear_regres_X_Y(er);
    }

    /*******************************************/

    if (selector.Eta)
    {
        std::cout << "######Eta calculation######" << std::endl;
        std::vector<std::array<double, 5>> eta;
        if (data_plane_vec.get_Ploidy() == 2)
        {
            eta = calc_eta(data_plane_vec);
        }
        else
        {
            eta = calc_eta_q1_version(data_plane_vec);
        }

        output_eta_stat_files(eta);
        output_exp_regr_eta_stat_files(eta);
    }

    if (stat)
    {
        std::cout << "-----Writing result-----" << std::endl;
        output_stat_files(selector, result);
    }

    return result;
}