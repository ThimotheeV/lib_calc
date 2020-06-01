//Genedeme constructor
//Pop, indiv, locus => locus, deme, indiv
//locus, deme, indiv
template <std::size_t ploidy>
data_plane_vec_c::data_plane_vec_c(genepop_input_c<ploidy> const &genedeme_data)
{
    Ploidy = ploidy;
    //size of different part
    Nbr_of_deme = genedeme_data.Genotype.size();
    Nbr_dist_class = genedeme_data.Nbr_dist_class;
    if (genedeme_data.Dist_btw_deme.size() > 0)
    {
        Dist_btw_deme = genedeme_data.Dist_btw_deme;
    }
    else
    {
        //If no genedeme_data.Dist_btw_deme; Dist_btw_deme become a identity "inverse" matrix
        Dist_btw_deme.resize(Nbr_of_deme);
        for (auto deme1 = 0; deme1 < Nbr_of_deme; ++deme1)
        {
            Dist_btw_deme[deme1].resize(Nbr_of_deme);
            for (auto deme2 = 0; deme2 < Nbr_of_deme; ++deme2)
            {
                Dist_btw_deme[deme1][deme2] = (deme1 != deme2);
            }
        }
    }

    if (genedeme_data.Dist_class_btw_deme.size() > 0)
    {
        Dist_class_btw_deme = genedeme_data.Dist_class_btw_deme;
    }
    //If no genedeme_data.Dist_class_btw_deme (mean no genedeme_data.Dist_btw_deme); Dist_class_btw_deme = become a identity "inverse" matrix;
    else
    {
        Dist_btw_deme = genedeme_data.Dist_btw_deme;
    }

    Nbr_of_indiv_per_deme.reserve(Nbr_of_deme);

    for (int i = 0; i < Nbr_of_deme; ++i)
    {
        Nbr_of_indiv_per_deme.push_back(genedeme_data.Genotype[i].size());
    }

    Locus_nbr = genedeme_data.Genotype[0][0].size();
    Cumul_nbr_of_indiv_per_deme.reserve(Nbr_of_deme);

    for (auto nbr_indiv : Nbr_of_indiv_per_deme)
    {
        Cumul_nbr_of_indiv_per_deme.push_back(Nbr_of_indiv_tot);
        Nbr_of_indiv_tot += nbr_indiv;
    }

    set_indiv_feature();

    Plane_vec.reserve(Nbr_of_indiv_tot * Locus_nbr);
    Allele_state_per_loc.resize(Locus_nbr);

    //To calc nc for estimate with missing values
    Nomiss_nbr_of_deme_per_loc.resize(Locus_nbr);
    Nomiss_nbr_of_gene_per_loc.resize(Locus_nbr);
    Nomiss_nbr_of_gene_per_loc_per_deme.resize(Locus_nbr);
    Nomiss_nbr_of_indiv_per_loc.resize(Locus_nbr);
    Nomiss_nbr_of_indiv_per_loc_per_deme.resize(Locus_nbr);
    Nomiss_indiv_bool_per_loc.resize(Nbr_of_indiv_tot, bin_vec(Locus_nbr));

    for (int locus = 0; locus < Locus_nbr; ++locus)
    {
        std::map<int, int> temp_count_allele_state;
        Nomiss_nbr_of_gene_per_loc_per_deme[locus].reserve(Nbr_of_deme);
        Nomiss_nbr_of_indiv_per_loc_per_deme[locus].reserve(Nbr_of_deme);
        auto indiv_general = 0;
        for (int deme = 0; deme < Nbr_of_deme; ++deme)
        {
            int nomiss_nbr_of_indiv_in_deme = Nbr_of_indiv_per_deme[deme];
            int nomiss_nbr_of_gene_in_deme = nomiss_nbr_of_indiv_in_deme * Ploidy;
            for (int indiv = 0; indiv < Nbr_of_indiv_per_deme[deme]; ++indiv)
            {
                bool missing_value = false;
                for (int gene = 0; gene < Ploidy; ++gene)
                {
                    int value = genedeme_data.Genotype[deme][indiv][locus][gene];
                    //To calc nc for estimate with missing values
                    if (value == 0)
                    {
                        --nomiss_nbr_of_gene_in_deme;
                        missing_value = true;
                    }
                    else
                    {
                        //emplace return a pair consisting of an iterator to the inserted element (or element already in place) and a bool denoting whether the insertion took place
                        auto pair = temp_count_allele_state.emplace(value, 1);
                        if (!pair.second)
                        {
                            pair.first->second += 1;
                        }
                    }

                    Plane_vec.push_back(value);
                }
                if (missing_value)
                {
                    Nomiss_indiv_bool_per_loc[indiv_general].insert(locus, 0);
                    --nomiss_nbr_of_indiv_in_deme;
                }
                ++indiv_general;
            }
            //To calc nc for estimate with missing values
            if (nomiss_nbr_of_gene_in_deme != 0)
            {
                Nomiss_nbr_of_deme_per_loc[locus] += 1;
                Nomiss_nbr_of_gene_per_loc[locus] += nomiss_nbr_of_gene_in_deme;
                Nomiss_nbr_of_gene_per_loc_per_deme[locus].push_back(nomiss_nbr_of_gene_in_deme);
                Nomiss_nbr_of_indiv_per_loc[locus] += nomiss_nbr_of_indiv_in_deme;
                Nomiss_nbr_of_indiv_per_loc_per_deme[locus].push_back(nomiss_nbr_of_indiv_in_deme);
            }
        }

        Allele_state_per_loc[locus].reserve(temp_count_allele_state.size());
        for (auto const &pair : temp_count_allele_state)
        {
            Allele_state_per_loc[locus].push_back({pair.first, pair.second});
        }
        //size = number of non void deme
        Nomiss_nbr_of_gene_per_loc_per_deme[locus].shrink_to_fit();
        Nomiss_nbr_of_indiv_per_loc_per_deme[locus].shrink_to_fit();
    }

    //Genetic map

    Dist_btw_locus.resize(Locus_nbr);
    for (auto locus1 = 0; locus1 < Locus_nbr; ++locus1)
    {
        Dist_btw_locus[locus1].resize(Locus_nbr);
        for (auto locus2 = 0; locus2 < Locus_nbr; ++locus2)
        {
            Dist_btw_locus[locus1][locus2] = (locus1 != locus2);
        }
    }
}