#include "data_plane_vec.hpp"

void data_plane_vec_c::set_indiv_feature()
{
    //Each indiv have attribut
    Indiv_feat.resize(Nbr_of_indiv_tot);
    auto indiv_feat_itr = Indiv_feat.begin();
    for (int pop = 0; pop < Nbr_of_pop; ++pop)
    {
        for (int indiv = 0; indiv < Nbr_of_indiv_per_pop[pop]; ++indiv)
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
int data_plane_vec_c::nbr_of_locus() const
{
    return Locus_nbr;
}
int data_plane_vec_c::nbr_of_pop() const
{
    return Nbr_of_pop;
}
int data_plane_vec_c::nbr_of_gene() const
{
    return Nbr_of_indiv_tot * Ploidy;
}
int data_plane_vec_c::nbr_of_indiv() const
{
    return Nbr_of_indiv_tot;
}
int data_plane_vec_c::nbr_of_indiv_per_pop(int nbr_of_pop) const
{
    return Nbr_of_indiv_per_pop[nbr_of_pop];
}
std::vector<int> const &data_plane_vec_c::cumul_nbr_of_indiv_per_pop() const
{
    return Cumul_nbr_of_indiv_per_pop;
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
std::vector<int> const &data_plane_vec_c::nomiss_nbr_of_gene_per_loc_per_pop(int locus) const
{
    return Nomiss_nbr_of_gene_per_loc_per_pop[locus];
}
int data_plane_vec_c::nomiss_nbr_of_gene_per_loc_per_pop(int locus, int pop) const
{
    return Nomiss_nbr_of_gene_per_loc_per_pop[locus][pop];
}
std::vector<int> const &data_plane_vec_c::nomiss_nbr_of_indiv_per_loc_per_pop(int locus) const
{
    return Nomiss_nbr_of_indiv_per_loc_per_pop[locus];
}
int data_plane_vec_c::nomiss_nbr_of_indiv_per_loc_per_pop(int locus, int pop) const
{
    return Nomiss_nbr_of_indiv_per_loc_per_pop[locus][pop];
}
int data_plane_vec_c::nomiss_nbr_of_pop_per_loc(int locus) const
{
    return Nomiss_nbr_of_pop_per_loc[locus];
}

int data_plane_vec_c::nbr_allele_per_loc(int locus) const
{
    return Allele_state_per_loc[locus].size();
}

std::vector<std::array<int, 2>> const &data_plane_vec_c::allele_state_per_loc(int locus) const
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
    return Plane_vec[(Nbr_of_indiv_tot * locus + Cumul_nbr_of_indiv_per_pop[pop] + indiv) * Ploidy + gene]; //(Ploidy - 1) if operator use in haploid return the same for gene 0 and gene 1
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
double data_plane_vec_c::dist_btw_pop(int dpv_gene_index1, int dpv_gene_index2) const
{
    if (dpv_gene_index1 == dpv_gene_index2)
    {
        return true;
    }
    auto const pop_gen1 = Indiv_feat[get_indiv(dpv_gene_index1)].Pop;
    auto const pop_gen2 = Indiv_feat[get_indiv(dpv_gene_index2)].Pop;

    return Dist_btw_pop[pop_gen1][pop_gen2];
}
int data_plane_vec_c::nbr_of_dist_class() const
{
    return Nbr_dist_class;
}

int data_plane_vec_c::dist_class_btw_pop(int dpv_gene_index1, int dpv_gene_index2) const
{
    if (dpv_gene_index1 == dpv_gene_index2)
    {
        return true;
    }
    auto const pop_gen1 = Indiv_feat[get_indiv(dpv_gene_index1)].Pop;
    auto const pop_gen2 = Indiv_feat[get_indiv(dpv_gene_index2)].Pop;

    return Dist_class_btw_pop[pop_gen1][pop_gen2];
}