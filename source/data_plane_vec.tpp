//Genepop constructor
//Pop, indiv, locus => locus, pop, indiv
//locus, pop, indiv
template <std::size_t ploidy>
data_plane_vec_c::data_plane_vec_c(genepop_input_c<ploidy> const &genepop_data)
{
    Ploidy = ploidy;
    //size of different part
    Pop_nbr = genepop_data.Genotype.size();
    if (genepop_data.Dist_btw_pop.size() > 0)
    {
        Dist_btw_pop = genepop_data.Dist_btw_pop;
    }
    else
    {
        //If no genepop_data.Dist_btw_pop; Dist_btw_pop become a identity "inverse" matrix
        Dist_btw_pop.resize(Pop_nbr);
        for (auto pop1 = 0; pop1 < Pop_nbr; ++pop1)
        {
            Dist_btw_pop[pop1].resize(Pop_nbr);
            for (auto pop2 = 0; pop2 < Pop_nbr; ++pop2)
            {
                Dist_btw_pop[pop1][pop2] = (pop1 != pop2);
            }
        }
    }

    Indiv_nbr_per_pop.reserve(Pop_nbr);

    for (int i = 0; i < Pop_nbr; ++i)
    {
        Indiv_nbr_per_pop.push_back(genepop_data.Genotype[i].size());
    }

    Locus_nbr = genepop_data.Genotype[0][0].size();
    Cumul_indiv_nbr_per_pop.reserve(Pop_nbr);

    for (auto nbr_indiv : Indiv_nbr_per_pop)
    {
        Cumul_indiv_nbr_per_pop.push_back(Indiv_nbr_tot);
        Indiv_nbr_tot += nbr_indiv;
    }

    set_indiv_feature();

    Plane_vec.reserve(Indiv_nbr_tot * Locus_nbr);
    Allele_state_per_loc.reserve(Locus_nbr);

    //To calc nc for estimate with missing values
    Nomiss_pop_nbr_per_loc.resize(Locus_nbr);
    Nomiss_gene_nbr_per_loc.resize(Locus_nbr);
    Nomiss_gene_nbr_per_loc_per_pop.resize(Locus_nbr);
    Nomiss_indiv_nbr_per_loc.resize(Locus_nbr);
    Nomiss_indiv_nbr_per_loc_per_pop.resize(Locus_nbr);
    Nomiss_indiv_bool_per_loc.resize(Indiv_nbr_tot, bin_vec(Locus_nbr));

    for (int locus = 0; locus < Locus_nbr; ++locus)
    {
        std::map<int, int> temp_count_allele_state;

        Nomiss_gene_nbr_per_loc_per_pop[locus].reserve(Pop_nbr);
        Nomiss_indiv_nbr_per_loc_per_pop[locus].reserve(Pop_nbr);
        auto indiv_general = 0;
        for (int pop = 0; pop < Pop_nbr; ++pop)
        {
            int nomiss_indiv_nbr_in_pop = Indiv_nbr_per_pop[pop];
            int nomiss_gene_nbr_in_pop = nomiss_indiv_nbr_in_pop * Ploidy;
            for (int indiv = 0; indiv < Indiv_nbr_per_pop[pop]; ++indiv)
            {
                bool missing_value = false;
                for (int gene = 0; gene < Ploidy; ++gene)
                {
                    int value = genepop_data.Genotype[pop][indiv][locus][gene];
                    //To calc nc for estimate with missing values
                    if (value == 0)
                    {
                        --nomiss_gene_nbr_in_pop;
                        missing_value = true;
                    }
                    else
                    {
                        temp_count_allele_state.try_emplace(value, 1);
                    }
                    
                    Plane_vec.push_back(value);
                }
                if (missing_value)
                {
                    Nomiss_indiv_bool_per_loc[indiv_general].insert(locus, 0);
                    --nomiss_indiv_nbr_in_pop;
                }
                ++indiv_general;
            }
            //To calc nc for estimate with missing values
            if (nomiss_gene_nbr_in_pop != 0)
            {
                Nomiss_pop_nbr_per_loc[locus] += 1;
                Nomiss_gene_nbr_per_loc[locus] += nomiss_gene_nbr_in_pop;
                Nomiss_gene_nbr_per_loc_per_pop[locus].push_back(nomiss_gene_nbr_in_pop);
                Nomiss_indiv_nbr_per_loc[locus] += nomiss_indiv_nbr_in_pop;
                Nomiss_indiv_nbr_per_loc_per_pop[locus].push_back(nomiss_indiv_nbr_in_pop);
            }
        }
        
        Allele_state_per_loc.push_back({temp_count_allele_state.begin()->first, (--temp_count_allele_state.end())->first, static_cast<int>(temp_count_allele_state.size())});
        //size = number of non void pop
        Nomiss_gene_nbr_per_loc_per_pop[locus].shrink_to_fit();
        Nomiss_indiv_nbr_per_loc_per_pop[locus].shrink_to_fit();
    }
}