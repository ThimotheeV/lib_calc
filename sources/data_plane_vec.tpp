//Genedeme constructor
//Pop, indiv, locus => locus, deme, indiv
//locus, deme, indiv
template <std::size_t ploidy>
data_plane_vec_c::data_plane_vec_c(genepop_input_c<ploidy> const &genepop_data)
{
    Ploidy = ploidy;

    //Chr dist and chr class dist between deme
    Nbr_of_locus = genepop_data.Genotype[0][0].size();
    Nbr_of_chr = genepop_data.Chr_dist_btw_loc.size();
    Chr_dist_class_nbr = genepop_data.Chr_dist_class_nbr;

    Chr_dist_btw_loc = std::vector<double>{};
    if (genepop_data.Chr_dist_btw_loc.size() > 0)
    {
        Nbr_of_loc_per_chr.resize(Nbr_of_chr);
        Cumul_nbr_of_loc_per_chr.resize(Nbr_of_chr);
        Chr_dist_btw_loc.resize(Nbr_of_locus * Nbr_of_locus);

        int loc_tot_count = 0;
        int loc_cumul_per_chr = 0;
        for (auto chr = 0; chr < Nbr_of_chr; ++chr)
        {
            int nbr_of_loc = genepop_data.Chr_dist_btw_loc[chr].size();
            Nbr_of_loc_per_chr[chr] = nbr_of_loc;
            Cumul_nbr_of_loc_per_chr[chr] = loc_cumul_per_chr;
            loc_cumul_per_chr += nbr_of_loc;

            for (auto loc1 = 0; loc1 < nbr_of_loc; ++loc1)
            {
                for (auto loc2 = 0; loc2 < nbr_of_loc; ++loc2)
                {
                    Chr_dist_btw_loc[loc_tot_count] = genepop_data.Chr_dist_btw_loc[chr][loc1][loc2];
                    ++loc_tot_count;
                }
            }
        }
    }
    else // If no map
    {
        //When locus are independant, each chr have one locus on it
        Nbr_of_chr = Nbr_of_locus;
        Nbr_of_loc_per_chr.resize(Nbr_of_chr);
        Cumul_nbr_of_loc_per_chr.resize(Nbr_of_chr);
        int cumul = 0;
        for (int i = 0; i < Nbr_of_chr; ++i)
        {
            Nbr_of_loc_per_chr[i] = 1;
            Cumul_nbr_of_loc_per_chr[i] = cumul;
            ++cumul;
        }
    }

    Chr_dist_class_btw_loc = std::vector<int>{};
    if (genepop_data.Chr_dist_class_btw_loc.size() > 0)
    {
        Chr_dist_class_btw_loc.resize(Nbr_of_locus * Nbr_of_locus);
        int count = 0;
        for (auto chr = 0; chr < Nbr_of_chr; ++chr)
        {
            for (auto loc1 = 0; loc1 < genepop_data.Chr_dist_class_btw_loc[chr].size(); ++loc1)
            {
                for (auto loc2 = 0; loc2 < genepop_data.Chr_dist_class_btw_loc[chr].size(); ++loc2)
                {
                    Chr_dist_class_btw_loc[count] = genepop_data.Chr_dist_class_btw_loc[chr][loc1][loc2];
                    ++count;
                }
            }
        }
    }

    // Calc around loc and chr

    //Geo dist and geo class dist between deme
    Nbr_of_deme = genepop_data.Genotype.size();
    Geo_dist_class_nbr = genepop_data.Geo_dist_class_nbr;
    if (genepop_data.Geo_dist_btw_deme.size() > 0)
    {
        Geo_dist_btw_deme.resize(Nbr_of_deme * Nbr_of_deme);
        int count = 0;
        for (auto deme1 = 0; deme1 < Nbr_of_deme; ++deme1)
        {
            for (auto deme2 = 0; deme2 < Nbr_of_deme; ++deme2)
            {
                Geo_dist_btw_deme[count] = genepop_data.Geo_dist_btw_deme[deme1][deme2];
                ++count;
            }
        }
    }
    else
    {
        //If no genepop_data.Geo_dist_btw_deme; Geo_dist_btw_deme become a identity "inverse" matrix
        Geo_dist_btw_deme.resize(Nbr_of_deme * Nbr_of_deme);
        int count = 0;
        for (auto deme1 = 0; deme1 < Nbr_of_deme; ++deme1)
        {
            for (auto deme2 = 0; deme2 < Nbr_of_deme; ++deme2)
            {
                Geo_dist_btw_deme[count] = (deme1 != deme2);
                ++count;
            }
        }
    }

    Geo_dist_class_btw_deme = std::vector<int>{};
    if (genepop_data.Geo_dist_class_btw_deme.size() > 0)
    {
        Geo_dist_class_btw_deme.resize(Geo_dist_btw_deme.size());
        int count = 0;
        for (auto deme1 = 0; deme1 < Nbr_of_deme; ++deme1)
        {
            for (auto deme2 = 0; deme2 < Nbr_of_deme; ++deme2)
            {
                Geo_dist_class_btw_deme[count] = genepop_data.Geo_dist_class_btw_deme[deme1][deme2];
                ++count;
            }
        }
    }

    Nbr_of_indiv_per_deme.reserve(Nbr_of_deme);

    for (int i = 0; i < Nbr_of_deme; ++i)
    {
        Nbr_of_indiv_per_deme.push_back(genepop_data.Genotype[i].size());
    }

    Cumul_nbr_of_indiv_per_deme.reserve(Nbr_of_deme);

    for (auto nbr_indiv : Nbr_of_indiv_per_deme)
    {
        Cumul_nbr_of_indiv_per_deme.push_back(Nbr_of_indiv_tot);
        Nbr_of_indiv_tot += nbr_indiv;
    }

    set_indiv_feature();

    //Nbr indiv * nbr locus
    Plane_vec.reserve(Nbr_of_indiv_tot * Nbr_of_locus * Ploidy);
    Allele_state_per_chr_per_loc.reserve(Nbr_of_locus);
    //To calc nc for estimate with missing values
    Nomiss_nbr_of_deme_per_chr_per_loc.reserve(Nbr_of_locus);
    Nomiss_nbr_of_gene_per_chr_per_loc.reserve(Nbr_of_locus);
    Nomiss_nbr_of_gene_per_chr_per_loc_per_deme.reserve(Nbr_of_locus * Nbr_of_deme);
    Nomiss_nbr_of_indiv_per_loc.resize(Nbr_of_locus);
    Nomiss_nbr_of_indiv_per_chr_per_loc_per_deme.reserve(Nbr_of_locus * Nbr_of_deme);
    Nomiss_per_indiv_per_loc.resize(Nbr_of_indiv_tot, bin_vec(Nbr_of_locus));

    Polymorph_locus_list_per_chr.resize(Nbr_of_chr);

    int locus_index_in_sample = 0;

    for (int chr = 0; chr < Nbr_of_chr; ++chr)
    {
        Polymorph_locus_list_per_chr[chr].reserve(Nbr_of_loc_per_chr[chr]);

        for (int locus_index_in_chr = 0; locus_index_in_chr < Nbr_of_loc_per_chr[chr]; ++locus_index_in_chr)
        {
            std::map<int, int> temp_count_allele_state;
            Nomiss_nbr_of_deme_per_chr_per_loc.push_back(0);
            Nomiss_nbr_of_gene_per_chr_per_loc.push_back(0);
            auto indiv_general = 0;
            for (int deme = 0; deme < Nbr_of_deme; ++deme)
            {
                int nomiss_nbr_of_indiv_in_deme = Nbr_of_indiv_per_deme[deme];
                int nomiss_nbr_of_gene_in_deme = nomiss_nbr_of_indiv_in_deme * Ploidy;
                for (int indiv_index_in_deme = 0; indiv_index_in_deme < Nbr_of_indiv_per_deme[deme]; ++indiv_index_in_deme)
                {
                    bool missing_value = false;
                    for (int gene_index_in_indiv = 0; gene_index_in_indiv < Ploidy; ++gene_index_in_indiv)
                    {
                        int value = genepop_data.Genotype[deme][indiv_index_in_deme][Cumul_nbr_of_loc_per_chr[chr] + locus_index_in_chr][gene_index_in_indiv];
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
                        Nomiss_per_indiv_per_loc[indiv_general].insert(locus_index_in_sample, 0);
                        --nomiss_nbr_of_indiv_in_deme;
                    }
                    ++indiv_general;
                }
                //To calc nc for estimate with missing values
                Nomiss_nbr_of_indiv_per_loc[locus_index_in_sample] += nomiss_nbr_of_indiv_in_deme;
                Nomiss_nbr_of_deme_per_chr_per_loc.back() += 1;
                Nomiss_nbr_of_gene_per_chr_per_loc.back() += nomiss_nbr_of_gene_in_deme;
                Nomiss_nbr_of_gene_per_chr_per_loc_per_deme.push_back(nomiss_nbr_of_gene_in_deme);
                Nomiss_nbr_of_indiv_per_chr_per_loc_per_deme.push_back(nomiss_nbr_of_indiv_in_deme);
            }
            if (temp_count_allele_state.size() > 1)
            {
                Polymorph_locus_list_per_chr[chr].push_back(locus_index_in_chr);
            }

            Allele_state_per_chr_per_loc.push_back(std::move(temp_count_allele_state));

            ++locus_index_in_sample;
        }
        Polymorph_locus_list_per_chr[chr].shrink_to_fit();
    }
}