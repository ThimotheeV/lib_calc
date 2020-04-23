//Genepop constructor
//Pop, indiv, locus => locus, pop, indiv
//locus, pop, indiv
template <std::size_t ploidy>
data_plane_vec_c::data_plane_vec_c(genepop_input_c<ploidy> const &genepop_data)
{
    Ploidy = ploidy;
    //size of different part
    Nbr_of_pop = genepop_data.Genotype.size();
    Nbr_dist_class = genepop_data.Nbr_dist_class;
    if (genepop_data.Dist_btw_pop.size() > 0)
    {
        Dist_btw_pop = genepop_data.Dist_btw_pop;
    }
    else
    {
        //If no genepop_data.Dist_btw_pop; Dist_btw_pop become a identity "inverse" matrix
        Dist_btw_pop.resize(Nbr_of_pop);
        for (auto pop1 = 0; pop1 < Nbr_of_pop; ++pop1)
        {
            Dist_btw_pop[pop1].resize(Nbr_of_pop);
            for (auto pop2 = 0; pop2 < Nbr_of_pop; ++pop2)
            {
                Dist_btw_pop[pop1][pop2] = (pop1 != pop2);
            }
        }
    }

    if (genepop_data.Dist_class_btw_pop.size() > 0)
    {
        Dist_class_btw_pop = genepop_data.Dist_class_btw_pop;
    }
    //If no genepop_data.Dist_class_btw_pop (mean no genepop_data.Dist_btw_pop); Dist_class_btw_pop = become a identity "inverse" matrix;
    else
    {
        Dist_btw_pop = genepop_data.Dist_btw_pop;
    }

    Nbr_of_indiv_per_pop.reserve(Nbr_of_pop);

    for (int i = 0; i < Nbr_of_pop; ++i)
    {
        Nbr_of_indiv_per_pop.push_back(genepop_data.Genotype[i].size());
    }

    Locus_nbr = genepop_data.Genotype[0][0].size();
    Cumul_nbr_of_indiv_per_pop.reserve(Nbr_of_pop);

    for (auto nbr_indiv : Nbr_of_indiv_per_pop)
    {
        Cumul_nbr_of_indiv_per_pop.push_back(Nbr_of_indiv_tot);
        Nbr_of_indiv_tot += nbr_indiv;
    }

    set_indiv_feature();

    Plane_vec.reserve(Nbr_of_indiv_tot * Locus_nbr);
    Allele_state_per_loc.resize(Locus_nbr);

    //To calc nc for estimate with missing values
    Nomiss_nbr_of_pop_per_loc.resize(Locus_nbr);
    Nomiss_nbr_of_gene_per_loc.resize(Locus_nbr);
    Nomiss_nbr_of_gene_per_loc_per_pop.resize(Locus_nbr);
    Nomiss_nbr_of_indiv_per_loc.resize(Locus_nbr);
    Nomiss_nbr_of_indiv_per_loc_per_pop.resize(Locus_nbr);
    Nomiss_indiv_bool_per_loc.resize(Nbr_of_indiv_tot, bin_vec(Locus_nbr));

    for (int locus = 0; locus < Locus_nbr; ++locus)
    {
        std::map<int, int> temp_count_allele_state;
        Nomiss_nbr_of_gene_per_loc_per_pop[locus].reserve(Nbr_of_pop);
        Nomiss_nbr_of_indiv_per_loc_per_pop[locus].reserve(Nbr_of_pop);
        auto indiv_general = 0;
        for (int pop = 0; pop < Nbr_of_pop; ++pop)
        {
            int nomiss_nbr_of_indiv_in_pop = Nbr_of_indiv_per_pop[pop];
            int nomiss_nbr_of_gene_in_pop = nomiss_nbr_of_indiv_in_pop * Ploidy;
            for (int indiv = 0; indiv < Nbr_of_indiv_per_pop[pop]; ++indiv)
            {
                bool missing_value = false;
                for (int gene = 0; gene < Ploidy; ++gene)
                {
                    int value = genepop_data.Genotype[pop][indiv][locus][gene];
                    //To calc nc for estimate with missing values
                    if (value == 0)
                    {
                        --nomiss_nbr_of_gene_in_pop;
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
                    --nomiss_nbr_of_indiv_in_pop;
                }
                ++indiv_general;
            }
            //To calc nc for estimate with missing values
            if (nomiss_nbr_of_gene_in_pop != 0)
            {
                Nomiss_nbr_of_pop_per_loc[locus] += 1;
                Nomiss_nbr_of_gene_per_loc[locus] += nomiss_nbr_of_gene_in_pop;
                Nomiss_nbr_of_gene_per_loc_per_pop[locus].push_back(nomiss_nbr_of_gene_in_pop);
                Nomiss_nbr_of_indiv_per_loc[locus] += nomiss_nbr_of_indiv_in_pop;
                Nomiss_nbr_of_indiv_per_loc_per_pop[locus].push_back(nomiss_nbr_of_indiv_in_pop);
            }
        }

        Allele_state_per_loc[locus].reserve(temp_count_allele_state.size());
        for (auto const &pair : temp_count_allele_state)
        {
            Allele_state_per_loc[locus].push_back({pair.first, pair.second});
        }
        //size = number of non void pop
        Nomiss_nbr_of_gene_per_loc_per_pop[locus].shrink_to_fit();
        Nomiss_nbr_of_indiv_per_loc_per_pop[locus].shrink_to_fit();
    }
}