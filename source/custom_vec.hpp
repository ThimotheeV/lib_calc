#pragma once

#include <tuple>
#include "input.hpp"

class data_plane_vec_c
{
private:
    std::vector<int> Plane_vec;
    int Pop_nbr{0};
    //Cumul sum of indiv for each pop
    int Sum_indiv_per_locus{0};
    std::vector<int> Indiv_nbr_per_pop;
    std::vector<int> Cumul_indiv_per_pop;
    int Locus_nbr{0};
    int Ploidy{1};

public:
    data_plane_vec_c() {};
    //Genepop constructor
    //Pop, indiv, locus => locus, pop, indiv
    //locus, pop, indiv
    template <std::size_t ploidy>
    data_plane_vec_c(genepop_input_c<ploidy> const &genepop_data)
    {
        Ploidy = ploidy;
        //size of different part
        Pop_nbr = genepop_data.Genotype.size();
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
    std::vector<int> const &cumul_indiv_nbr_per_pop();
    std::vector<int> const &get_plane_vec();

    int operator[](int i) const;

    std::vector<int>::const_iterator begin() const;
    std::vector<int>::const_iterator end() const;

    int const &operator()(int locus, int pop, int indiv, int gene) const;
    int index_begin_locus(int locus) const;
    int index_end_locus(int locus) const;

    // in the same locus => gene 1 & gene 2 : same pop, indiv ?
    bool same_indiv(int gene1, int gene2) const;
    bool same_pop(int locus, int gene1, int gene2) const;
};