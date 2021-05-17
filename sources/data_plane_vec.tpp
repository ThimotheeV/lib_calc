//Genedeme constructor
//Pop, indiv, locus => locus, deme, indiv
//locus, deme, indiv
template <std::size_t ploidy>
data_plane_vec_c::data_plane_vec_c(genepop_input_c<ploidy> const &genedeme_data)
{
    Ploidy = ploidy;
    //Dist matrix locus
    Dist_btw_loc = genedeme_data.Dist_btw_loc;
    Dist_class_btw_loc = genedeme_data.Dist_class_btw_loc;
    Nbr_chr_dist_class = genedeme_data.Nbr_chr_dist_class;

    Nbr_of_locus = genedeme_data.Genotype[0][0].size();
    if (Dist_btw_loc.size() == 0)
    {
        //When locus are independant, each chr have one locus on it
        Chr_nbr = Nbr_of_locus;
        Nbr_of_loc_per_chr.resize(Chr_nbr);
        Cumul_nbr_of_loc_per_chr.resize(Chr_nbr);
        int cumul = 0;
        for (int i = 0; i < Chr_nbr; ++i)
        {
            Nbr_of_loc_per_chr[i] = 1;
            Cumul_nbr_of_loc_per_chr[i] = cumul;
            ++cumul;
        }
    }
    else
    {
        Chr_nbr = Dist_btw_loc.size();
        Nbr_of_loc_per_chr.resize(Chr_nbr);
        Cumul_nbr_of_loc_per_chr.resize(Chr_nbr);
        int cumul = 0;
        for (int i = 0; i < Chr_nbr; ++i)
        {
            int nbr_of_loc = Dist_btw_loc[i].size();
            Nbr_of_loc_per_chr[i] = nbr_of_loc;
            Cumul_nbr_of_loc_per_chr[i] = cumul;
            cumul += nbr_of_loc;
        }
    }

    //size of different part
    Nbr_of_deme = genedeme_data.Genotype.size();
    Dist_class_nbr = genedeme_data.Nbr_dist_class;
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

    Cumul_nbr_of_indiv_per_deme.reserve(Nbr_of_deme);

    for (auto nbr_indiv : Nbr_of_indiv_per_deme)
    {
        Cumul_nbr_of_indiv_per_deme.push_back(Nbr_of_indiv_tot);
        Nbr_of_indiv_tot += nbr_indiv;
    }

    set_indiv_feature();

    //Nbr indiv * nbr locus
    Plane_vec.reserve(Nbr_of_indiv_tot * Nbr_of_locus * Ploidy);
    Allele_state_per_chr_per_loc.resize(Chr_nbr);
    //To calc nc for estimate with missing values
    Nomiss_nbr_of_deme_per_chr_per_loc.resize(Chr_nbr);
    Nomiss_nbr_of_gene_per_chr_per_loc.resize(Chr_nbr);
    Nomiss_nbr_of_gene_per_chr_per_loc_per_deme.resize(Chr_nbr);
    Nomiss_nbr_of_indiv_per_loc.resize(Nbr_of_locus);
    Nomiss_nbr_of_indiv_per_chr_per_loc_per_deme.resize(Chr_nbr);
    Nomiss_indiv_bool_per_loc.resize(Nbr_of_indiv_tot, bin_vec(Nbr_of_locus));

    //WARNING : Be sure than min state < 2 ^ 16
    Allele_state_bound = std::array<int, 2>{1 << 16, 0};
    Polymorph_locus_per_chr.resize(Chr_nbr);

    int locus_general = 0;

    for (int chr = 0; chr < Chr_nbr; ++chr)
    {
        Allele_state_per_chr_per_loc[chr].resize(Nbr_of_loc_per_chr[chr]);
        Nomiss_nbr_of_deme_per_chr_per_loc[chr].resize(Nbr_of_loc_per_chr[chr]);
        Nomiss_nbr_of_gene_per_chr_per_loc[chr].resize(Nbr_of_loc_per_chr[chr]);
        Nomiss_nbr_of_gene_per_chr_per_loc_per_deme[chr].resize(Nbr_of_loc_per_chr[chr]);
        Nomiss_nbr_of_indiv_per_chr_per_loc_per_deme[chr].resize(Nbr_of_loc_per_chr[chr]);
        Polymorph_locus_per_chr[chr].reserve(Nbr_of_loc_per_chr[chr]);

        for (int locus = 0; locus < Nbr_of_loc_per_chr[chr]; ++locus)
        {
            std::map<int, int> temp_count_allele_state;
            Nomiss_nbr_of_gene_per_chr_per_loc_per_deme[chr][locus].reserve(Nbr_of_deme);
            Nomiss_nbr_of_indiv_per_chr_per_loc_per_deme[chr][locus].reserve(Nbr_of_deme);
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
                        int value = genedeme_data.Genotype[deme][indiv][Cumul_nbr_of_loc_per_chr[chr] + locus][gene];
                        //To calc nc for estimate with missing values
                        if (value == 0)
                        {
                            --nomiss_nbr_of_gene_in_deme;
                            missing_value = true;
                        }
                        else
                        {
                            if (value < Allele_state_bound.at(0))
                            {
                                Allele_state_bound.at(0) = value;
                            }
                            if (value > Allele_state_bound.at(1))
                            {
                                Allele_state_bound.at(1) = value;
                            }
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
                        Nomiss_indiv_bool_per_loc[indiv_general].insert(locus_general, 0);
                        --nomiss_nbr_of_indiv_in_deme;
                    }
                    ++indiv_general;
                }
                //To calc nc for estimate with missing values
                if (nomiss_nbr_of_gene_in_deme != 0)
                {
                    Nomiss_nbr_of_indiv_per_loc[locus_general] += nomiss_nbr_of_indiv_in_deme;
                    Nomiss_nbr_of_deme_per_chr_per_loc[chr][locus] += 1;
                    Nomiss_nbr_of_gene_per_chr_per_loc[chr][locus] += nomiss_nbr_of_gene_in_deme;
                    Nomiss_nbr_of_gene_per_chr_per_loc_per_deme[chr][locus].push_back(nomiss_nbr_of_gene_in_deme);
                    Nomiss_nbr_of_indiv_per_chr_per_loc_per_deme[chr][locus].push_back(nomiss_nbr_of_indiv_in_deme);
                }
            }
            if (temp_count_allele_state.size() > 1)
            {
                Polymorph_locus_per_chr[chr].push_back(locus);
            }

            Allele_state_per_chr_per_loc[chr][locus] = std::move(temp_count_allele_state);
            //size = number of non void deme
            Nomiss_nbr_of_gene_per_chr_per_loc_per_deme[chr][locus].shrink_to_fit();
            Nomiss_nbr_of_indiv_per_chr_per_loc_per_deme[chr][locus].shrink_to_fit();

            ++locus_general;
        }
        Polymorph_locus_per_chr[chr].shrink_to_fit();
        //Genetic map
    }
}