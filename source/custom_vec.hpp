#pragma once

#include <tuple>
#include "input.hpp"

struct feature_c
{
    int Pop{9999};
};

class data_plane_vec_c
{
private:
    std::vector<int> Plane_vec;
    std::vector<feature_c> Indiv_feat;
    std::vector<std::vector<int>> Dist_btw_pop;
    int Pop_nbr{0};
    //Cumul sum of indiv for each pop
    int Sum_indiv_per_locus{0};
    std::vector<int> Indiv_nbr_per_pop;
    std::vector<int> Cumul_indiv_per_pop;
    int Locus_nbr{0};
    int Ploidy{1};

public:
    data_plane_vec_c(){};
    //Genepop constructor
    //Pop, indiv, locus => locus, pop, indiv
    //locus, pop, indiv
    template <std::size_t ploidy>
    data_plane_vec_c(genepop_input_c<ploidy> const &genepop_data)
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
            //If no genepop_data.Dist_btw_pop; Dist_btw_pop become a identity matrix  
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
        Cumul_indiv_per_pop.reserve(Pop_nbr);

        for (auto nbr_indiv : Indiv_nbr_per_pop)
        {
            Cumul_indiv_per_pop.push_back(Sum_indiv_per_locus);
            Sum_indiv_per_locus += nbr_indiv;
        }

        //TODO : Separate in a other function
        //Each indiv have attribut
        Indiv_feat.resize(Sum_indiv_per_locus);
        auto indiv_feat_itr = Indiv_feat.begin();
        for (int pop = 0; pop < Pop_nbr; ++pop)
        {
            for (int indiv = 0; indiv < Indiv_nbr_per_pop[pop]; ++indiv)
            {
                indiv_feat_itr->Pop = pop;
                ++indiv_feat_itr;
            }
        }

        Plane_vec.reserve(Sum_indiv_per_locus * Locus_nbr);

        for (int locus = 0; locus < Locus_nbr; ++locus)
        {
            for (int pop = 0; pop < Pop_nbr; ++pop)
            {
                for (int indiv = 0; indiv < Indiv_nbr_per_pop[pop]; ++indiv)
                {
                    for (int gene = 0; gene < Ploidy; ++gene)
                    {
                        Plane_vec.push_back(genepop_data.Genotype[pop][indiv][locus].at(gene));
                    }
                }
            }
        }
    }

    //Getter
    int size() const;
    int pop_nbr() const;                      //number of pop
    int indiv_nbr_per_pop(int pop_nbr) const; //pop size
    int locus_nbr() const;                    //number of locus in a specifique indiv
    int indiv_nbr() const;                    //
    feature_c const &get_feature(int indiv);
    int get_indiv(int gene) const;
    std::vector<int> const &cumul_indiv_nbr_per_pop();
    std::vector<int> const &get_plane_vec();

    int operator[](int i) const;

    std::vector<int>::const_iterator begin() const;
    std::vector<int>::const_iterator end() const;

    int const &operator()(int locus, int pop, int indiv, int gene) const;
    int index_begin_locus(int locus) const;
    int index_end_locus(int locus) const;

    // in the same locus => gene 1 & gene 2 : same pop, indiv ?
    bool same_indiv(int gene_index1, int gene_index2) const;
    bool pop_at_dist(int gene_index1, int gene_index2, int dist) const;
};