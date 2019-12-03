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
    return Plane_vec.cbegin();
}
std::vector<int>::const_iterator data_plane_vec_c::end() const
{
    return Plane_vec.cend();
}

int const &data_plane_vec_c::operator()(int locus, int pop, int indiv, int gene) const
{
    if(gene > Ploidy - 1)
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

bool data_plane_vec_c::same_indiv(int gene1, int gene2) const
{
    if (gene1 == gene2)
    {
        return true;
    }
    if (gene1 > gene2)
    {
        auto temp = gene1;
        gene1 = gene2;
        gene2 = temp;
    }

    bool result{false};
    if (Ploidy == 2)
    {
        result = (gene1 % 2 == 0) && (gene2 - gene1 == 1);
    }
    return result;
}

bool data_plane_vec_c::same_pop(int locus, int gene1, int gene2) const
{
    if (gene1 == gene2)
    {
        return true;
    }

    if (gene1 > gene2)
    {
        auto temp = gene1;
        gene1 = gene2;
        gene2 = temp;
    }

    int cumul_pop_size = Indiv_nbr_per_pop[0] * Ploidy + index_begin_locus(locus);
    int end_locus = index_end_locus(locus);
    int nbr_pop = 0;

    //Gene 1 come from nbr_pop - 1 ?
    while (cumul_pop_size <= end_locus)
    {
        if (gene1 < cumul_pop_size)
        {
            if (gene2 < cumul_pop_size)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        //cumul for each locus
        ++nbr_pop;
        nbr_pop %= Pop_nbr;
        cumul_pop_size += Indiv_nbr_per_pop[nbr_pop] * Ploidy;
    }

    return false;
}