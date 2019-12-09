#include "custom_vec.hpp"

int data_plane_vec_c::size() const
{
    return Plane_vec.size();
}

int data_plane_vec_c::pop_nbr() const
{
    return Pop_nbr;
}
int data_plane_vec_c::indiv_nbr_per_pop(int pop_nbr) const
{
    return Indiv_nbr_per_pop[pop_nbr];
}
int data_plane_vec_c::locus_nbr() const
{
    return Locus_nbr;
}

int data_plane_vec_c::indiv_nbr() const
{
    return Sum_indiv_per_locus;
}

feature_c const &data_plane_vec_c::get_feature(int indiv)
{
    return Indiv_feat[indiv];
}

int data_plane_vec_c::get_indiv(int gene_index) const
{
    int place_in_locus = gene_index % (Sum_indiv_per_locus * Ploidy);

    return place_in_locus / Ploidy;
}

std::vector<int> const &data_plane_vec_c::cumul_indiv_nbr_per_pop()
{
    return Cumul_indiv_per_pop;
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
    return Plane_vec[(Sum_indiv_per_locus * locus + Cumul_indiv_per_pop[pop] + indiv) * Ploidy + gene]; //(Ploidy - 1) if operator use in haploid return the same for gene 0 and gene 1
}

int data_plane_vec_c::index_begin_locus(int locus) const
{
    return Sum_indiv_per_locus * locus * Ploidy;
}

int data_plane_vec_c::index_end_locus(int locus) const
{
    return Sum_indiv_per_locus * (locus + 1) * Ploidy;
}

bool data_plane_vec_c::same_indiv(int gene_index1, int gene_index2) const
{
    if (gene_index1 == gene_index2)
    {
        return true;
    }
    if (gene_index1 > gene_index2)
    {
        auto temp = gene_index1;
        gene_index1 = gene_index2;
        gene_index2 = temp;
    }

    bool result{false};
    if (Ploidy == 2)
    {
        result = (gene_index1 % 2 == 0) && (gene_index2 - gene_index1 == 1);
    }
    return result;
}

//Passer par un tableau d'attribut des indivs
bool data_plane_vec_c::pop_at_dist(int gene_index1, int gene_index2, int dist) const
{
    if (gene_index1 == gene_index2)
    {
        return true;
    }
    auto const pop_gen1 = Indiv_feat[get_indiv(gene_index1)].Pop;
    auto const pop_gen2 = Indiv_feat[get_indiv(gene_index2)].Pop;

    return Dist_btw_pop[pop_gen1][pop_gen2] == dist;
}