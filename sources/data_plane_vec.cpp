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
int data_plane_vec_c::nbr_locus() const
{
    return Nbr_of_locus;
}
int data_plane_vec_c::nbr_locus(int chr) const
{
    return Nbr_of_loc_per_chr[chr];
}
int data_plane_vec_c::nbr_of_deme() const
{
    return Nbr_of_deme;
}
int data_plane_vec_c::nbr_of_gene() const
{
    return Nbr_of_indiv_tot * Nbr_of_locus * Ploidy;
}
int data_plane_vec_c::nbr_of_chr() const
{
    return Chr_nbr;
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

int data_plane_vec_c::nomiss_nbr_of_gene(int chr, int locus) const
{
    return Nomiss_nbr_of_gene_per_chr_per_loc[chr][locus];
}
int data_plane_vec_c::nomiss_nbr_of_indiv(int locus) const
{
    return Nomiss_nbr_of_indiv_per_loc[locus];
}
int data_plane_vec_c::nomiss_nbr_of_indiv(int chr, int locus) const
{
    return Nomiss_nbr_of_indiv_per_loc[Cumul_nbr_of_loc_per_chr[chr] + locus];
}
std::vector<int> const &data_plane_vec_c::nomiss_nbr_of_gene_per_deme(int chr, int locus) const
{
    return Nomiss_nbr_of_gene_per_chr_per_loc_per_deme[chr][locus];
}
int data_plane_vec_c::nomiss_nbr_of_gene(int chr, int locus, int deme) const
{
    return Nomiss_nbr_of_gene_per_chr_per_loc_per_deme[chr][locus][deme];
}
std::vector<int> const &data_plane_vec_c::nomiss_nbr_of_indiv_per_deme(int chr, int locus) const
{
    return Nomiss_nbr_of_indiv_per_chr_per_loc_per_deme[chr][locus];
}
int data_plane_vec_c::nomiss_nbr_of_indiv(int chr, int locus, int deme) const
{
    return Nomiss_nbr_of_indiv_per_chr_per_loc_per_deme[chr][locus][deme];
}
int data_plane_vec_c::nomiss_nbr_of_deme(int chr, int locus) const
{
    return Nomiss_nbr_of_deme_per_chr_per_loc[chr][locus];
}

int data_plane_vec_c::nbr_allele(int chr, int locus) const
{
    return Allele_state_per_chr_per_loc[chr][locus].size();
}

//map(state, nbr of allele in this state)
std::map<int, int> const &data_plane_vec_c::allele_state(int chr, int locus) const
{
    return Allele_state_per_chr_per_loc[chr][locus];
}

int data_plane_vec_c::state_min() const
{
    return Allele_state_bound.at(0);
}

int data_plane_vec_c::state_max() const
{
    return Allele_state_bound.at(1);
}

std::vector<int> const &data_plane_vec_c::polymorph_locus(int chr) const
{
    return Polymorph_locus_per_chr[chr];
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

int const &data_plane_vec_c::operator()(int chr, int locus, int deme, int indiv, int gene) const
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
    //Cumul_nbr_of_loc_per_chr[chr] * Nbr_indiv * Ploidy + locus * Nbr_indiv * Ploidy + Cumul_nbr_of_indiv_per_deme[deme] * Ploidy + indiv(deme relative index) * Ploidy + gene
    return Plane_vec[((Cumul_nbr_of_loc_per_chr[chr] + locus) * Nbr_of_indiv_tot + Cumul_nbr_of_indiv_per_deme[deme] + indiv) * Ploidy + gene]; //(Ploidy - 1) if operator use in haploid return the same for gene 0 and gene 1
}

int const &data_plane_vec_c::operator()(int chr, int locus, int indiv, int gene) const
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
    //Cumul_nbr_of_loc_per_chr[chr] * Nbr_indiv * Ploidy + locus * Nbr_indiv * Ploidy + indiv(absolute index) * Ploidy + gene
    return Plane_vec[((Cumul_nbr_of_loc_per_chr[chr] + locus) * Nbr_of_indiv_tot + indiv) * Ploidy + gene]; //(Ploidy - 1) if operator use in haploid return the same for gene 0 and gene 1
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
    //Cumul_nbr_of_loc_per_chr[chr] * Nbr_indiv * Ploidy + locus * Nbr_indiv * Ploidy + indiv(absolute index) * Ploidy + gene
    return Plane_vec[(locus * Nbr_of_indiv_tot + indiv) * Ploidy + gene]; //(Ploidy - 1) if operator use in haploid return the same for gene 0 and gene 1
}

int data_plane_vec_c::index_begin_locus(int chr, int locus) const
{
    return (Cumul_nbr_of_loc_per_chr[chr] + locus) * Nbr_of_indiv_tot * Ploidy;
}

int data_plane_vec_c::index_end_locus(int chr, int locus) const
{
    return (Cumul_nbr_of_loc_per_chr[chr] + locus + 1) * Nbr_of_indiv_tot * Ploidy;
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

bin_vec const &data_plane_vec_c::nomiss_data_indiv_per_loc(int indiv) const
{
    return Nomiss_indiv_bool_per_loc[indiv];
}

bool data_plane_vec_c::nomiss_data_indiv(int locus, int indiv) const
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

int data_plane_vec_c::dist_class_btw_deme_with_deme(int deme1, int deme2) const
{
    if (deme1 == deme2)
    {
        return 0;
    }

    return Dist_class_btw_deme[deme1][deme2];
}

int data_plane_vec_c::nbr_of_chr_dist_class() const
{
    return Nbr_chr_dist_class;
}

double data_plane_vec_c::dist_btw_locus(int chr, int locus_index1, int locus_index2) const
{
    if (locus_index1 == locus_index2)
    {
        return 0;
    }

    return Dist_btw_loc[chr][locus_index1][locus_index2];
}

int data_plane_vec_c::dist_class_btw_locus(int chr, int locus_index1, int locus_index2) const
{
    if (locus_index1 == locus_index2)
    {
        return 0;
    }

    return Dist_class_btw_loc[chr][locus_index1][locus_index2];
}