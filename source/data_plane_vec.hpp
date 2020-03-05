#pragma once

#include <tuple>
#include <map>
#include <cmath>

#include "common_tools.hpp"
#include "input.hpp"

struct feature_c
{
    int Pop{-1};
};

class data_plane_vec_c
{
private:
    std::vector<int> Plane_vec;
    std::vector<feature_c> Indiv_feat;
    std::vector<std::vector<float>> Dist_btw_pop;
    //Locus, pop
    int Pop_nbr{0};
    //Cumul sum of indiv for each pop
    int Indiv_nbr_tot{0};
    std::vector<int> Indiv_nbr_per_pop;
    std::vector<int> Cumul_indiv_nbr_per_pop;
    int Locus_nbr{0};
    //Complete indiv sampled
    std::vector<int> Nomiss_gene_nbr_per_loc;
    std::vector<int> Nomiss_indiv_nbr_per_loc;
    std::vector<int> Nomiss_pop_nbr_per_loc;
    //matrix(indiv, loc)
    std::vector<bin_vec> Nomiss_indiv_bool_per_loc;
    std::vector<std::vector<int>> Nomiss_gene_nbr_per_loc_per_pop;
    std::vector<std::vector<int>> Nomiss_indiv_nbr_per_loc_per_pop;
    //Allele state resum by loc (min, max, nbr of state)
    std::vector<std::array<int, 3>> Allele_state_per_loc;
    int Ploidy{-1};

public:
    data_plane_vec_c(){};
    //Genepop constructor
    //Pop, indiv, locus => locus, pop, indiv
    //locus, pop, indiv
    template <std::size_t ploidy>
    data_plane_vec_c(genepop_input_c<ploidy> const &genepop_data);

    //Setter
    void set_indiv_feature();

    //Getter
    int get_Ploidy() const;

    int size() const;
    int pop_nbr() const;
    int locus_nbr() const;                    //number of locus in a specifique indiv
    int indiv_nbr() const;                    //number of pop
    int indiv_nbr_per_pop(int pop_nbr) const; //pop size
    std::vector<int> const &cumul_indiv_nbr_per_pop() const;

    int get_indiv(int gene) const;
    feature_c const &get_feature(int indiv);

    int nomiss_gene_nbr_per_loc(int locus) const;
    int nomiss_indiv_nbr_per_loc(int locus) const;
    std::vector<int> const &nomiss_gene_nbr_per_loc_per_pop(int locus) const;
    int nomiss_gene_nbr_per_loc_per_pop(int locus, int pop) const;
    std::vector<int> const &nomiss_indiv_nbr_per_loc_per_pop(int locus) const;
    int nomiss_indiv_nbr_per_loc_per_pop(int locus, int pop) const;
    int nomiss_pop_nbr_per_loc(int locus) const;

    std::array<int, 3> const &allele_state_per_loc(int locus) const;

    std::vector<int> const &get_plane_vec();
    int operator[](int i) const;

    std::vector<int>::const_iterator begin() const;
    std::vector<int>::const_iterator end() const;

    int const &operator()(int locus, int pop, int indiv, int gene) const;
    //Acces at indiv for one specific locus (pop independent)
    int const &operator()(int locus, int indiv, int gene) const;
    int index_begin_locus(int locus) const;
    int index_end_locus(int locus) const;

    // in the same locus => gene 1 & gene 2 : same pop, indiv ?
    bool same_indiv(int dpv_gene_index1, int dpv_gene_index2) const;
    bool same_pop(int dpv_gene_index1, int dpv_gene_index2) const;
    bin_vec const &nomiss_data_indiv(int indiv) const;
    bool nomiss_data_indiv_per_loc(int indiv, int locus) const;
    float dist_btw_pop(int dpv_gene_index1, int dpv_gene_index2) const;
};

#include "data_plane_vec.tpp"