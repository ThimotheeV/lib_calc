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
int data_plane_vec_c::nbr_of_gene_per_loc() const
{
    return Nbr_of_indiv_tot * Ploidy;
}
int data_plane_vec_c::nbr_of_chr() const
{
    return Nbr_of_chr;
}
int data_plane_vec_c::nbr_of_indiv() const
{
    return Nbr_of_indiv_tot;
}
int data_plane_vec_c::nbr_of_indiv(int nbr_of_deme) const
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
feature_c const &data_plane_vec_c::get_feature(int indiv_index_in_sample)
{
    return Indiv_feat[indiv_index_in_sample];
}

int data_plane_vec_c::nomiss_nbr_of_gene(int chr, int locus_index_in_chr) const
{
    return Nomiss_nbr_of_gene_per_chr_per_loc[Cumul_nbr_of_loc_per_chr[chr] + locus_index_in_chr];
}
int data_plane_vec_c::nomiss_nbr_of_gene(int chr, int locus_index_in_chr, int deme) const
{
    return Nomiss_nbr_of_gene_per_chr_per_loc_per_deme[(Cumul_nbr_of_loc_per_chr[chr] + locus_index_in_chr) * Nbr_of_deme + deme];
}
int data_plane_vec_c::nomiss_nbr_of_indiv(int locus) const
{
    return Nomiss_nbr_of_indiv_per_loc[locus];
}
int data_plane_vec_c::nomiss_nbr_of_indiv(int chr, int locus_index_in_chr) const
{
    return Nomiss_nbr_of_indiv_per_loc[Cumul_nbr_of_loc_per_chr[chr] + locus_index_in_chr];
}
int data_plane_vec_c::nomiss_nbr_of_indiv(int chr, int locus_index_in_chr, int deme) const
{
    return Nomiss_nbr_of_indiv_per_chr_per_loc_per_deme[(Cumul_nbr_of_loc_per_chr[chr] + locus_index_in_chr) * Nbr_of_deme + deme];
}
int data_plane_vec_c::nomiss_nbr_of_deme(int chr, int locus_index_in_chr) const
{
    return Nomiss_nbr_of_deme_per_chr_per_loc[Cumul_nbr_of_loc_per_chr[chr] + locus_index_in_chr];
}

int data_plane_vec_c::nbr_allele(int chr, int locus_index_in_chr) const
{
    return Allele_state_per_chr_per_loc[Cumul_nbr_of_loc_per_chr[chr] + locus_index_in_chr].size();
}

//map(state, nbr of allele in this state)
std::map<int, int> const &data_plane_vec_c::allele_state(int chr, int locus_index_in_chr) const
{
    return Allele_state_per_chr_per_loc[Cumul_nbr_of_loc_per_chr[chr] + locus_index_in_chr];
}

std::vector<int> const &data_plane_vec_c::polymorph_locus_list(int chr) const
{
    return Polymorph_locus_list_per_chr[chr];
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

int const &data_plane_vec_c::operator()(int chr, int locus_index_in_chr, int deme, int indiv_index_in_deme, int gene_index_in_indiv) const
{
    if (gene_index_in_indiv > Ploidy - 1)
    {
        ++gene_index_in_indiv;
        throw std::logic_error("Can't show the gene " + std::to_string(gene_index_in_indiv) + " when max gene by indiv is " + std::to_string(Ploidy));
    }

    if (indiv_index_in_deme > Nbr_of_indiv_tot - 1)
    {
        throw std::logic_error("Only " + std::to_string(Nbr_of_indiv_tot) + " in locus " + std::to_string(locus_index_in_chr));
    }
    //Cumul_nbr_of_loc_per_chr[chr] * Nbr_indiv * Ploidy + locus * Nbr_indiv * Ploidy + Cumul_nbr_of_indiv_per_deme[deme] * Ploidy + indiv(deme relative index) * Ploidy + gene
    return Plane_vec[((Cumul_nbr_of_loc_per_chr[chr] + locus_index_in_chr) * Nbr_of_indiv_tot + Cumul_nbr_of_indiv_per_deme[deme] + indiv_index_in_deme) * Ploidy + gene_index_in_indiv]; //(Ploidy - 1) if operator use in haploid return the same for gene 0 and gene 1
}

int const &data_plane_vec_c::operator()(int chr, int locus_index_in_chr, int indiv_index_in_sample, int gene_index_in_indiv) const
{
    if (gene_index_in_indiv > Ploidy - 1)
    {
        ++gene_index_in_indiv;
        throw std::logic_error("Can't show the gene " + std::to_string(gene_index_in_indiv) + " when max gene by indiv is " + std::to_string(Ploidy));
    }

    if (indiv_index_in_sample > Nbr_of_indiv_tot - 1)
    {
        throw std::logic_error("Only " + std::to_string(Nbr_of_indiv_tot) + " in locus " + std::to_string(locus_index_in_chr));
    }
    //Cumul_nbr_of_loc_per_chr[chr] * Nbr_indiv * Ploidy + locus * Nbr_indiv * Ploidy + indiv(absolute index) * Ploidy + gene
    return Plane_vec[((Cumul_nbr_of_loc_per_chr[chr] + locus_index_in_chr) * Nbr_of_indiv_tot + indiv_index_in_sample) * Ploidy + gene_index_in_indiv]; //(Ploidy - 1) if operator use in haploid return the same for gene 0 and gene 1
}

int const &data_plane_vec_c::operator()(int locus_index_in_sample, int indiv_index_in_sample, int gene_index_in_indiv) const
{
    if (gene_index_in_indiv > Ploidy - 1)
    {
        ++gene_index_in_indiv;
        throw std::logic_error("Can't show the gene " + std::to_string(gene_index_in_indiv) + " when max gene by indiv is " + std::to_string(Ploidy));
    }

    if (indiv_index_in_sample > Nbr_of_indiv_tot - 1)
    {
        throw std::logic_error("Only " + std::to_string(Nbr_of_indiv_tot) + " in locus " + std::to_string(locus_index_in_sample));
    }
    //Cumul_nbr_of_loc_per_chr[chr] * Nbr_indiv * Ploidy + locus * Nbr_indiv * Ploidy + indiv(absolute index) * Ploidy + gene
    return Plane_vec[(locus_index_in_sample * Nbr_of_indiv_tot + indiv_index_in_sample) * Ploidy + gene_index_in_indiv]; //(Ploidy - 1) if operator use in haploid return the same for gene 0 and gene 1
}

int data_plane_vec_c::index_begin_locus(int chr, int locus_index_in_chr) const
{
    return (Cumul_nbr_of_loc_per_chr[chr] + locus_index_in_chr) * Nbr_of_indiv_tot * Ploidy;
}

int data_plane_vec_c::index_end_locus(int chr, int locus_index_in_chr) const
{
    return (Cumul_nbr_of_loc_per_chr[chr] + locus_index_in_chr + 1) * Nbr_of_indiv_tot * Ploidy;
}

//In same locus
bool data_plane_vec_c::same_indiv(int gene1_index_in_sample, int gene2_index_in_sample) const
{
    if (gene1_index_in_sample == gene2_index_in_sample)
    {
        return true;
    }
    if (gene1_index_in_sample > gene2_index_in_sample)
    {
        auto temp = gene1_index_in_sample;
        gene1_index_in_sample = gene2_index_in_sample;
        gene2_index_in_sample = temp;
    }

    bool result{false};
    if (Ploidy == 2)
    {
        //Need to be adjacent
        result = (gene1_index_in_sample % 2 == 0) && (gene2_index_in_sample - gene1_index_in_sample == 1);
    }
    return result;
}

bool data_plane_vec_c::same_deme(int gene1_index_in_sample, int gene2_index_in_sample) const
{
    if (gene1_index_in_sample == gene2_index_in_sample)
    {
        return true;
    }

    return Indiv_feat[get_indiv(gene1_index_in_sample)].Deme == Indiv_feat[get_indiv(gene2_index_in_sample)].Deme;
}

bin_vec const &data_plane_vec_c::nomiss_data(int indiv_index_in_sample) const
{
    return Nomiss_per_indiv_per_loc[indiv_index_in_sample];
}

//Passer par un tableau d'attribut des indivs
double data_plane_vec_c::geo_dist_btw_gene(int gene1_index_in_sample, int gene2_index_in_sample) const
{
    if (gene1_index_in_sample == gene2_index_in_sample)
    {
        return 0;
    }
    auto const deme_gen1 = Indiv_feat[get_indiv(gene1_index_in_sample)].Deme;
    auto const deme_gen2 = Indiv_feat[get_indiv(gene2_index_in_sample)].Deme;

    return Geo_dist_btw_deme[deme_gen1 * Nbr_of_deme + deme_gen2];
}

double data_plane_vec_c::geo_dist_btw_deme(int deme1, int deme2) const
{
    if (deme1 == deme2)
    {
        return 0;
    }

    return Geo_dist_btw_deme[deme1 * Nbr_of_deme + deme2];
}

int data_plane_vec_c::nbr_geo_dist_class() const
{
    return Geo_dist_class_nbr;
}

int data_plane_vec_c::geo_dist_class_btw_gene(int gene1_index_in_sample, int gene2_index_in_sample) const
{
    if (gene1_index_in_sample == gene2_index_in_sample)
    {
        return 0;
    }
    auto const deme_gen1 = Indiv_feat[get_indiv(gene1_index_in_sample)].Deme;
    auto const deme_gen2 = Indiv_feat[get_indiv(gene2_index_in_sample)].Deme;

    return Geo_dist_class_btw_deme[deme_gen1 * Nbr_of_deme + deme_gen2];
}

int data_plane_vec_c::geo_dist_class_btw_deme(int deme1, int deme2) const
{
    if (deme1 == deme2)
    {
        return 0;
    }

    return Geo_dist_class_btw_deme[deme1 * Nbr_of_deme + deme2];
}

int data_plane_vec_c::nbr_chr_dist_class() const
{
    return Chr_dist_class_nbr;
}

double data_plane_vec_c::chr_dist_btw_locus(int chr, int locus1_index_in_chr, int locus2_index_in_chr) const
{
    if (locus1_index_in_chr == locus2_index_in_chr)
    {
        return 0;
    }

    return Chr_dist_btw_loc[Cumul_nbr_of_loc_per_chr[chr] + locus1_index_in_chr * Nbr_of_loc_per_chr[chr] + locus2_index_in_chr];
}

int data_plane_vec_c::chr_dist_class_btw_locus(int chr, int locus1_index_in_chr, int locus2_index_in_chr) const
{
    if (locus1_index_in_chr == locus2_index_in_chr)
    {
        return 0;
    }

    return Chr_dist_class_btw_loc[Cumul_nbr_of_loc_per_chr[chr] + locus1_index_in_chr * Nbr_of_loc_per_chr[chr] + locus2_index_in_chr];
}