#include "data_plane_vec.hpp"

void data_plane_vec_c::set_indiv_feature()
{
    //Each indiv have attribut
    Indiv_feat.resize(Nbr_of_indiv_tot);
    auto indiv_feat_itr = Indiv_feat.begin();
    for (int deme = 0; deme < Nbr_of_deme; ++deme)
    {
        for (int indiv = 0; indiv < Nbr_of_indiv_per_deme[deme]; ++indiv)
        {
            indiv_feat_itr->Deme = deme;
            ++indiv_feat_itr;
        }
    }
}

int data_plane_vec_c::get_Ploidy() const
{
    return Ploidy;
}

int data_plane_vec_c::size() const
{
    return Plane_vec.size();
}
int data_plane_vec_c::base_nbr_locus_per_indiv() const
{
    return Locus_nbr;
}
int data_plane_vec_c::nbr_of_deme() const
{
    return Nbr_of_deme;
}
int data_plane_vec_c::nbr_of_locus_tot() const
{
    return Nbr_of_indiv_tot * Ploidy * Locus_nbr;
}
int data_plane_vec_c::nbr_of_indiv() const
{
    return Nbr_of_indiv_tot;
}
int data_plane_vec_c::nbr_of_indiv_per_deme(int nbr_of_deme) const
{
    return Nbr_of_indiv_per_deme[nbr_of_deme];
}
std::vector<int> const &data_plane_vec_c::cumul_nbr_of_indiv_per_deme() const
{
    return Cumul_nbr_of_indiv_per_deme;
}

int data_plane_vec_c::get_indiv(int gene_index) const
{
    int place_in_locus = gene_index % (Nbr_of_indiv_tot * Ploidy);

    return place_in_locus / Ploidy;
}
feature_c const &data_plane_vec_c::get_feature(int indiv)
{
    return Indiv_feat[indiv];
}

int data_plane_vec_c::nomiss_nbr_of_gene_per_loc(int locus) const
{
    return Nomiss_nbr_of_gene_per_loc[locus];
}
int data_plane_vec_c::nomiss_nbr_of_indiv_per_loc(int locus) const
{
    return Nomiss_nbr_of_indiv_per_loc[locus];
}
std::vector<int> const &data_plane_vec_c::nomiss_nbr_of_gene_per_loc_per_deme(int locus) const
{
    return Nomiss_nbr_of_gene_per_loc_per_deme[locus];
}
int data_plane_vec_c::nomiss_nbr_of_gene_per_loc_per_deme(int locus, int deme) const
{
    return Nomiss_nbr_of_gene_per_loc_per_deme[locus][deme];
}
std::vector<int> const &data_plane_vec_c::nomiss_nbr_of_indiv_per_loc_per_deme(int locus) const
{
    return Nomiss_nbr_of_indiv_per_loc_per_deme[locus];
}
int data_plane_vec_c::nomiss_nbr_of_indiv_per_loc_per_deme(int locus, int deme) const
{
    return Nomiss_nbr_of_indiv_per_loc_per_deme[locus][deme];
}
int data_plane_vec_c::nomiss_nbr_of_deme_per_loc(int locus) const
{
    return Nomiss_nbr_of_deme_per_loc[locus];
}

int data_plane_vec_c::nbr_allele_per_loc(int locus) const
{
    return Allele_state_per_loc[locus].size();
}

//map(state, nbr of allele in this state)
std::map<int, int> const &data_plane_vec_c::allele_state_per_loc(int locus) const
{
    return Allele_state_per_loc[locus];
}

int data_plane_vec_c::state_min() const
{
    return Allele_state_bound.at(0);
}

int data_plane_vec_c::state_max() const
{
    return Allele_state_bound.at(1);
}

std::vector<int> const &data_plane_vec_c::get_plane_vec()
{
    return Plane_vec;
}

int data_plane_vec_c::operator[](int i) const
{
    return Plane_vec[i];
}

std::vector<int>::const_iterator data_plane_vec_c::begin() const
{
    //cbegin => const_iter
    return Plane_vec.cbegin();
}
std::vector<int>::const_iterator data_plane_vec_c::end() const
{
    return Plane_vec.cend();
}

int const &data_plane_vec_c::operator()(int locus, int deme, int indiv, int gene) const
{
    if (gene > Ploidy - 1)
    {
        ++gene;
        throw std::logic_error("Can't show the gene " + std::to_string(gene) + " when max gene by indiv is " + std::to_string(Ploidy));
    }
    return Plane_vec[(Nbr_of_indiv_tot * locus + Cumul_nbr_of_indiv_per_deme[deme] + indiv) * Ploidy + gene]; //(Ploidy - 1) if operator use in haploid return the same for gene 0 and gene 1
}

int const &data_plane_vec_c::operator()(int locus, int indiv, int gene) const
{
    if (gene > Ploidy - 1)
    {
        ++gene;
        throw std::logic_error("Can't show the gene " + std::to_string(gene) + " when max gene by indiv is " + std::to_string(Ploidy));
    }

    if (indiv > Nbr_of_indiv_tot - 1)
    {
        throw std::logic_error("Only " + std::to_string(Nbr_of_indiv_tot) + " in locus " + std::to_string(locus));
    }

    return Plane_vec[(Nbr_of_indiv_tot * locus + indiv) * Ploidy + gene]; //(Ploidy - 1) if operator use in haploid return the same for gene 0 and gene 1
}

int data_plane_vec_c::index_begin_locus(int locus) const
{
    return Nbr_of_indiv_tot * locus * Ploidy;
}

int data_plane_vec_c::index_end_locus(int locus) const
{
    return Nbr_of_indiv_tot * (locus + 1) * Ploidy;
}

//In same locus
bool data_plane_vec_c::same_loc_in_indiv(int dpv_gene_index1, int dpv_gene_index2) const
{
    if (dpv_gene_index1 == dpv_gene_index2)
    {
        return true;
    }
    if (dpv_gene_index1 > dpv_gene_index2)
    {
        auto temp = dpv_gene_index1;
        dpv_gene_index1 = dpv_gene_index2;
        dpv_gene_index2 = temp;
    }

    bool result{false};
    if (Ploidy == 2)
    {
        result = (dpv_gene_index1 % 2 == 0) && (dpv_gene_index2 - dpv_gene_index1 == 1);
    }
    return result;
}

bool data_plane_vec_c::same_deme(int dpv_gene_index1, int dpv_gene_index2) const
{
    if (dpv_gene_index1 == dpv_gene_index2)
    {
        return true;
    }

    return Indiv_feat[get_indiv(dpv_gene_index1)].Deme == Indiv_feat[get_indiv(dpv_gene_index2)].Deme;
}

bin_vec const &data_plane_vec_c::nomiss_data_indiv(int indiv) const
{
    return Nomiss_indiv_bool_per_loc[indiv];
}

bool data_plane_vec_c::nomiss_data_indiv_per_loc(int indiv, int locus) const
{
    return Nomiss_indiv_bool_per_loc[indiv].at(locus);
}

//Passer par un tableau d'attribut des indivs
double data_plane_vec_c::dist_btw_deme(int dpv_gene_index1, int dpv_gene_index2) const
{
    if (dpv_gene_index1 == dpv_gene_index2)
    {
        return 0;
    }
    auto const deme_gen1 = Indiv_feat[get_indiv(dpv_gene_index1)].Deme;
    auto const deme_gen2 = Indiv_feat[get_indiv(dpv_gene_index2)].Deme;

    return Dist_btw_deme[deme_gen1][deme_gen2];
}

double data_plane_vec_c::dist_btw_deme_with_deme(int deme_index1, int deme_index2) const
{
    if (deme_index1 == deme_index2)
    {
        return 0;
    }

    return Dist_btw_deme[deme_index1][deme_index2];
}

int data_plane_vec_c::nbr_of_dist_class() const
{
    return Dist_class_nbr;
}

int data_plane_vec_c::dist_class_btw_deme(int dpv_gene_index1, int dpv_gene_index2) const
{
    if (dpv_gene_index1 == dpv_gene_index2)
    {
        return 0;
    }
    auto const deme_gen1 = Indiv_feat[get_indiv(dpv_gene_index1)].Deme;
    auto const deme_gen2 = Indiv_feat[get_indiv(dpv_gene_index2)].Deme;

    return Dist_class_btw_deme[deme_gen1][deme_gen2];
}

double data_plane_vec_c::dist_btw_locus(int locus_index1, int locus_index2) const
{
    if (locus_index1 == locus_index2)
    {
        return 0;
    }

    return Dist_btw_locus[locus_index1][locus_index2];
}