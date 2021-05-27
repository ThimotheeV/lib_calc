#pragma once

#include <list>
#include <tuple>
#include <map>
#include <cmath>

#include "common_tools.hpp"
#include "input.hpp"

struct feature_c
{
    int Deme{-1};
};

class data_plane_vec_c
{
    //Needed to have wrapper in GSpace
protected:
    std::vector<int> Plane_vec;
    std::vector<feature_c> Indiv_feat;
    std::vector<std::vector<double>> Dist_btw_deme;
    std::vector<std::vector<int>> Dist_class_btw_deme;
    int Dist_class_nbr = 0;

    //Chr<matrix(loc_i, loc_j)>
    std::vector<std::vector<std::vector<double>>> Dist_btw_loc;
    std::vector<std::vector<std::vector<int>>> Dist_class_btw_loc;
    int Nbr_chr_dist_class = 0;

    //Locus, deme
    int Nbr_of_deme{0};
    //Cumul sum of indiv for each deme
    int Nbr_of_indiv_tot{0};
    std::vector<int> Nbr_of_indiv_per_deme;
    std::vector<int> Cumul_nbr_of_indiv_per_deme;
    int Nbr_of_chr{0};
    //Each chr can have a heterogene number of locus, each loc have the same number of deme and each deme can have a heterogene number of indiv
    int Nbr_of_locus{0};
    std::vector<int> Nbr_of_loc_per_chr;
    std::vector<int> Cumul_nbr_of_loc_per_chr;
    //Complete indiv sampled
    std::vector<std::vector<int>> Nomiss_nbr_of_gene_per_chr_per_loc;
    //Work with absolute locus index
    std::vector<int> Nomiss_nbr_of_indiv_per_loc;
    std::vector<std::vector<int>> Nomiss_nbr_of_deme_per_chr_per_loc;
    //for each indiv => if locus have (ploidy) data present (absolute locus index)
    std::vector<bin_vec> Nomiss_per_indiv_per_loc;
    //[chr][locus][deme]
    std::vector<std::vector<std::vector<int>>> Nomiss_nbr_of_gene_per_chr_per_loc_per_deme;
    std::vector<std::vector<std::vector<int>>> Nomiss_nbr_of_indiv_per_chr_per_loc_per_deme;
    //Allele state resum by loc (state, nbr of allele in this state)
    std::vector<std::vector<std::map<int, int>>> Allele_state_per_chr_per_loc;
    //Chr<Indent of polymorphic locus>
    std::vector<std::vector<int>> Polymorph_locus_per_chr;
    //Total allele bound(min, max)
    std::array<int, 2> Allele_state_bound;
    int Ploidy{-1};

public:
    data_plane_vec_c(){};
    //Genedeme constructor
    //Pop, indiv, locus => chr, locus, deme, indiv
    //locus, deme, indiv
    template <std::size_t ploidy>
    data_plane_vec_c(genepop_input_c<ploidy> const &genepop_data);

    //Setter
    void set_indiv_feature();

    //Getter
    int get_Ploidy() const;

    int size() const;
    int nbr_of_deme() const;
    int nbr_locus() const;        //number of locus in a specifique indiv
    int nbr_locus(int chr) const; //number of locus in a specifique indiv
    int nbr_of_gene() const;
    int nbr_of_chr() const;
    int nbr_of_indiv() const;
    int nbr_of_indiv_per_deme(int nbr_of_deme) const; //deme size
    std::vector<int> const &cumul_nbr_of_indiv_per_deme() const;

    int get_indiv(int gene) const;
    feature_c const &get_feature(int indiv);

    int nomiss_nbr_of_gene(int chr, int locus) const;
    //in haploid data, indiv = gene ; in diploid data, indiv = 2 genes
    //absolute locus index
    int nomiss_nbr_of_indiv(int locus) const;
    // chr relative locus index
    int nomiss_nbr_of_indiv(int chr, int locus) const;
    std::vector<int> const &nomiss_nbr_of_gene_per_deme(int chr, int locus) const;
    int nomiss_nbr_of_gene(int chr, int locus, int deme) const;
    std::vector<int> const &nomiss_nbr_of_indiv_per_deme(int chr, int locus) const;
    int nomiss_nbr_of_indiv(int chr, int locus, int deme) const;
    int nomiss_nbr_of_deme(int chr, int locus) const;

    int nbr_allele(int chr, int locus) const;
    std::map<int, int> const &allele_state(int chr, int locus) const;
    int state_min() const;
    int state_max() const;
    std::vector<int> const &polymorph_locus(int chr) const;

    std::vector<int> const &get_plane_vec();
    int operator[](int i) const;

    std::vector<int>::const_iterator begin() const;
    std::vector<int>::const_iterator end() const;

    int const &operator()(int chr, int locus, int deme, int indiv, int gene) const;
    //Acces at indiv for one specific locus (deme independent)
    int const &operator()(int chr, int locus, int indiv, int gene) const;
    int const &operator()(int locus, int indiv, int gene) const;
    int index_begin_locus(int chr, int locus) const;
    int index_end_locus(int chr, int locus) const;

    // in the same locus => gene 1 & gene 2 : same deme, indiv ?
    bool same_indiv(int dpv_gene_index1, int dpv_gene_index2) const;
    bool same_deme(int dpv_gene_index1, int dpv_gene_index2) const;
    bin_vec const &nomiss_data(int indiv) const;
    double geo_dist_btw_gene(int dpv_gene_index1, int dpv_gene_index2) const;
    double geo_dist_btw_deme(int deme_index1, int deme_index2) const;
    int nbr_geo_dist_class() const;
    int geo_dist_class_btw_gene(int dpv_gene_index1, int dpv_gene_index2) const;
    int geo_dist_class_btw_deme(int deme1, int deme2) const;

    int nbr_chr_dist_class() const;
    double chr_dist_btw_locus(int chr, int locus_index1, int locus_index2) const;
    int chr_dist_class_btw_locus(int chr, int locus_index1, int locus_index2) const;
};

#include "data_plane_vec.tpp"