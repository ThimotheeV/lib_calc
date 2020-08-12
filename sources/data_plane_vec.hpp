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
    std::vector<std::vector<double>> Dist_btw_deme;
    std::vector<std::vector<int>> Dist_class_btw_deme;
    int Nbr_dist_class = 10;

    std::vector<std::vector<double>> Dist_btw_locus;
    //Locus, deme
    int Nbr_of_deme{0};
    //Cumul sum of indiv for each deme
    int Nbr_of_indiv_tot{0};
    std::vector<int> Nbr_of_indiv_per_deme;
    std::vector<int> Cumul_nbr_of_indiv_per_deme;
    int Locus_nbr{0};
    //Complete indiv sampled
    std::vector<int> Nomiss_nbr_of_gene_per_loc;
    std::vector<int> Nomiss_nbr_of_indiv_per_loc;
    std::vector<int> Nomiss_nbr_of_deme_per_loc;
    //matrix(indiv, loc)
    std::vector<bin_vec> Nomiss_indiv_bool_per_loc;
    //[locus][deme]
    std::vector<std::vector<int>> Nomiss_nbr_of_gene_per_loc_per_deme;
    std::vector<std::vector<int>> Nomiss_nbr_of_indiv_per_loc_per_deme;
    //Allele state resum by loc (min, max, nbr of state)
    std::vector<std::vector<std::array<int, 2>>> Allele_state_per_loc;
    int Ploidy{-1};

public:
    data_plane_vec_c(){};
    //Genedeme constructor
    //Pop, indiv, locus => locus, deme, indiv
    //locus, deme, indiv
    template <std::size_t ploidy>
    data_plane_vec_c(genepop_input_c<ploidy> const &genedeme_data);

    //Setter
    void set_indiv_feature();

    //Getter
    int get_Ploidy() const;

    int size() const;
    int nbr_of_deme() const;
    int nbr_of_locus() const; //number of locus in a specifique indiv
    int nbr_of_gene() const;
    int nbr_of_indiv() const;
    int nbr_of_indiv_per_deme(int nbr_of_deme) const; //deme size
    std::vector<int> const &cumul_nbr_of_indiv_per_deme() const;

    int get_indiv(int gene) const;
    feature_c const &get_feature(int indiv);

    int nomiss_nbr_of_gene_per_loc(int locus) const;
    int nomiss_nbr_of_indiv_per_loc(int locus) const;
    std::vector<int> const &nomiss_nbr_of_gene_per_loc_per_deme(int locus) const;
    int nomiss_nbr_of_gene_per_loc_per_deme(int locus, int deme) const;
    std::vector<int> const &nomiss_nbr_of_indiv_per_loc_per_deme(int locus) const;
    int nomiss_nbr_of_indiv_per_loc_per_deme(int locus, int deme) const;
    int nomiss_nbr_of_deme_per_loc(int locus) const;

    int nbr_allele_per_loc(int locus) const;
    std::vector<std::array<int, 2>> const &allele_state_per_loc(int locus) const;

    std::vector<int> const &get_plane_vec();
    int operator[](int i) const;

    std::vector<int>::const_iterator begin() const;
    std::vector<int>::const_iterator end() const;

    int const &operator()(int locus, int deme, int indiv, int gene) const;
    //Acces at indiv for one specific locus (deme independent)
    int const &operator()(int locus, int indiv, int gene) const;
    int index_begin_locus(int locus) const;
    int index_end_locus(int locus) const;

    // in the same locus => gene 1 & gene 2 : same deme, indiv ?
    bool same_indiv(int dpv_gene_index1, int dpv_gene_index2) const;
    bool same_deme(int dpv_gene_index1, int dpv_gene_index2) const;
    bin_vec const &nomiss_data_indiv(int indiv) const;
    bool nomiss_data_indiv_per_loc(int indiv, int locus) const;
    double dist_btw_deme(int dpv_gene_index1, int dpv_gene_index2) const;
    double dist_btw_deme_with_deme(int deme_index1, int deme_index2) const;
    int nbr_of_dist_class() const;
    int dist_class_btw_deme(int dpv_gene_index1, int dpv_gene_index2) const;

    double dist_btw_locus(int locus_index1, int locus_index2) const;
};

#include "data_plane_vec.tpp"