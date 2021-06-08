#include "run.hpp"

result_c run(selector_input_c const &selector, data_plane_vec_c const &data_plane_vec)
{
    result_c result(selector.Geo_dist_class_nbr);
    bool stat = false;
    if (selector.Hobs)
    {
        std::cout << "\n######Hobs calculation######" << std::endl;
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
        std::cout << "\n######Hexp calculation######" << std::endl;
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
        std::cout << "\n######Var calculation######" << std::endl;
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
        std::cout << "\n######Nb_allele calculation######" << std::endl;
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
        std::cout << "\n######MGW calculation######" << std::endl;
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
        std::cout << "\n######SFS calculation######" << std::endl;
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
        std::string str(selector.Generic_data_filename + "_SFS.txt");
        output_sfs_stat_files(output_map, str);
    }

    /*******************************************/

    if (selector.F_stat)
    {
        std::cout << "\n######F_stat calculation######" << std::endl;
        stat = true;
        //Fis or Fst

        auto temp = Fstat_genepop(data_plane_vec, selector.Missing_data);
        result.Fis = temp.at(0);
        result.Fst = temp.at(1);
    }

    /*******************************************/

    if (selector.Qr)
    {
        std::cout << "\n######Qr calculation######" << std::endl;
        stat = true;
        result.Qr = calc_qr_all_loc(data_plane_vec);
    }

    /*******************************************/

    if (selector.Ar)
    {
        std::cout << "\n######Ar calculation######" << std::endl;
        stat = true;
        auto Ar = ar_by_pair(data_plane_vec);
        result.Ar_reg = linear_regres_X_Y(Ar);
    }

    /*******************************************/

    if (selector.Er)
    {
        std::cout << "\n######Er calculation######" << std::endl;
        stat = true;
        auto er = er_by_pair(data_plane_vec);
        result.Er_reg = linear_regres_X_Y(er);
    }

    /*******************************************/

    if (selector.Eta)
    {
        std::cout << "\n######Eta calculation######" << std::endl;
        std::vector<std::array<double, 5>> eta;

        if (data_plane_vec.nbr_of_indiv() == data_plane_vec.nbr_of_deme())
            if (data_plane_vec.get_Ploidy() == 2)
            {
                eta = calc_eta(data_plane_vec);
            }
            else
            {
                eta = calc_eta_q1_version(data_plane_vec);
            }
        else
        {
            eta = calc_eta_1_indiv_deme_v(data_plane_vec);
        }
        std::string str(selector.Generic_data_filename + "_Eta.txt");
        output_eta_stat_files(eta, str);
        //output_exp_regr_eta_stat_files(eta);
    }

    if (stat)
    {
        std::string str(selector.Generic_data_filename + "_Stats.txt");
        output_stat_files(selector, result, str);
    }

    return result;
}