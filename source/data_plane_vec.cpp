#include "data_plane_vec.hpp"

void data_plane_vec_c::set_indiv_feature()
{
    //Each indiv have attribut
    Indiv_feat.resize(Indiv_nbr_tot);
    auto indiv_feat_itr = Indiv_feat.begin();
    for (int pop = 0; pop < Pop_nbr; ++pop)
    {
        for (int indiv = 0; indiv < Indiv_nbr_per_pop[pop]; ++indiv)
        {
            indiv_feat_itr->Pop = pop;
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
int data_plane_vec_c::locus_nbr() const
{
    return Locus_nbr;
}
int data_plane_vec_c::pop_nbr() const
{
    return Pop_nbr;
}
int data_plane_vec_c::indiv_nbr() const
{
    return Indiv_nbr_tot;
}
int data_plane_vec_c::indiv_nbr_per_pop(int pop_nbr) const
{
    return Indiv_nbr_per_pop[pop_nbr];
}
std::vector<int> const &data_plane_vec_c::cumul_indiv_nbr_per_pop() const
{
    return Cumul_indiv_nbr_per_pop;
}

int data_plane_vec_c::get_indiv(int gene_index) const
{
    int place_in_locus = gene_index % (Indiv_nbr_tot * Ploidy);

    return place_in_locus / Ploidy;
}
feature_c const &data_plane_vec_c::get_feature(int indiv)
{
    return Indiv_feat[indiv];
}

int data_plane_vec_c::nomiss_gene_nbr_per_loc(int locus) const
{
    return Nomiss_gene_nbr_per_loc[locus];
}
int data_plane_vec_c::nomiss_indiv_nbr_per_loc(int locus) const
{
    return Nomiss_indiv_nbr_per_loc[locus];
}
std::vector<int> const &data_plane_vec_c::nomiss_gene_nbr_per_loc_per_pop(int locus) const
{
    return Nomiss_gene_nbr_per_loc_per_pop[locus];
}
int data_plane_vec_c::nomiss_gene_nbr_per_loc_per_pop(int locus, int pop) const
{
    return Nomiss_gene_nbr_per_loc_per_pop[locus][pop];
}
std::vector<int> const &data_plane_vec_c::nomiss_indiv_nbr_per_loc_per_pop(int locus) const
{
    return Nomiss_indiv_nbr_per_loc_per_pop[locus];
}
int data_plane_vec_c::nomiss_indiv_nbr_per_loc_per_pop(int locus, int pop) const
{
    return Nomiss_indiv_nbr_per_loc_per_pop[locus][pop];
}
int data_plane_vec_c::nomiss_pop_nbr_per_loc(int locus) const
{
    return Nomiss_pop_nbr_per_loc[locus];
}

std::array<int, 3> const &data_plane_vec_c::allele_state_per_loc(int locus) const
{
    return Allele_state_per_loc[locus];
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

int const &data_plane_vec_c::operator()(int locus, int pop, int indiv, int gene) const
{
    if (gene > Ploidy - 1)
    {
        ++gene;
        throw std::logic_error("Can't show the gene " + std::to_string(gene) + " when max gene by indiv is " + std::to_string(Ploidy));
    }
    return Plane_vec[(Indiv_nbr_tot * locus + Cumul_indiv_nbr_per_pop[pop] + indiv) * Ploidy + gene]; //(Ploidy - 1) if operator use in haploid return the same for gene 0 and gene 1
}

int const &data_plane_vec_c::operator()(int locus, int indiv, int gene) const
{
    if (gene > Ploidy - 1)
    {
        ++gene;
        throw std::logic_error("Can't show the gene " + std::to_string(gene) + " when max gene by indiv is " + std::to_string(Ploidy));
    }

    if (indiv > Indiv_nbr_tot - 1)
    {
        throw std::logic_error("Only " + std::to_string(Indiv_nbr_tot) + " in locus " + std::to_string(locus));
    }

    return Plane_vec[(Indiv_nbr_tot * locus + indiv) * Ploidy + gene]; //(Ploidy - 1) if operator use in haploid return the same for gene 0 and gene 1
}

int data_plane_vec_c::index_begin_locus(int locus) const
{
    return Indiv_nbr_tot * locus * Ploidy;
}

int data_plane_vec_c::index_end_locus(int locus) const
{
    return Indiv_nbr_tot * (locus + 1) * Ploidy;
}

bool data_plane_vec_c::same_indiv(int dpv_gene_index1, int dpv_gene_index2) const
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

bool data_plane_vec_c::same_pop(int dpv_gene_index1, int dpv_gene_index2) const
{
    if (dpv_gene_index1 == dpv_gene_index2)
    {
        return true;
    }

    return Indiv_feat[get_indiv(dpv_gene_index1)].Pop == Indiv_feat[get_indiv(dpv_gene_index2)].Pop;
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
float data_plane_vec_c::dist_btw_pop(int dpv_gene_index1, int dpv_gene_index2) const
{
    if (dpv_gene_index1 == dpv_gene_index2)
    {
        return true;
    }
    auto const pop_gen1 = Indiv_feat[get_indiv(dpv_gene_index1)].Pop;
    auto const pop_gen2 = Indiv_feat[get_indiv(dpv_gene_index2)].Pop;

    return Dist_btw_pop[pop_gen1][pop_gen2];
}