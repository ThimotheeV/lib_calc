#include <tuple>

#include "calc_stat.hpp"

//calc_Q0 => calc_Q0_intra_ind
double calc_Q0(data_plane_vec_c const &genotype)
{
    auto Q0 = std::array<int, 2>{0, 0};
    for (int i = 0; i < genotype.size(); i += 2)
    {
        if (genotype[i] == genotype[i + 1])
        {
            ++Q0.at(0);
        }
        ++Q0.at(1);
    }
    return static_cast<double>(Q0.at(0)) / Q0.at(1);
}

//pop/locus/indiv
double calc_Q1(data_plane_vec_c const &genotype)
{
    auto Q1 = std::array<int, 2>{0, 0};
    for (int locus = 0; locus < genotype.locus_nbr(); ++locus)
    {
        for (int pop = 0; pop < genotype.pop_nbr(); ++pop)
        {
            for (int indiv = 0; indiv < genotype.indiv_nbr_per_pop(pop) - 1; ++indiv)
            {
                for (int next_indiv = indiv + 1; next_indiv < genotype.indiv_nbr_per_pop(pop); ++next_indiv)
                {

                    if (genotype(locus, pop, indiv, 0) == genotype(locus, pop, indiv, 1))
                    {
                        ++Q1.at(0);
                    }
                    if (genotype(locus, pop, indiv, 0) == genotype(locus, pop, next_indiv, 0))
                    {
                        ++Q1.at(0);
                    }
                    if (genotype(locus, pop, indiv, 0) == genotype(locus, pop, next_indiv, 1))
                    {
                        ++Q1.at(0);
                    }
                    if (genotype(locus, pop, indiv, 1) == genotype(locus, pop, next_indiv, 0))
                    {
                        ++Q1.at(0);
                    }
                    if (genotype(locus, pop, indiv, 1) == genotype(locus, pop, next_indiv, 1))
                    {
                        ++Q1.at(0);
                    }
                    Q1.at(1) += 5;
                }
            }
            //Last locus
            if (genotype(locus, pop, genotype.indiv_nbr_per_pop(pop) - 1, 0) == genotype(locus, pop, genotype.indiv_nbr_per_pop(pop) - 1, 1))
            {
                ++Q1.at(0);
            }
            ++Q1.at(1);
        }
    }
    return static_cast<double>(Q1.at(0)) / Q1.at(1);
}

//pop/locus/indiv
double calc_Q2(data_plane_vec_c const &genotype)
{
    auto Q2 = std::array<int, 2>{0, 0};
    for (int locus = 0; locus < genotype.locus_nbr(); ++locus)
    {
        for (int pop = 0; pop < genotype.pop_nbr() - 1; ++pop)
        {
            for (int other_pop = pop + 1; other_pop < genotype.pop_nbr(); ++other_pop)
            {
                for (int indiv = 0; indiv < genotype.indiv_nbr_per_pop(pop); ++indiv)
                {
                    for (int indiv_other_pop = 0; indiv_other_pop < genotype.indiv_nbr_per_pop(other_pop); ++indiv_other_pop)
                    {
                        if (genotype(locus, pop, indiv, 0) == genotype(locus, other_pop, indiv_other_pop, 0))
                        {
                            ++Q2.at(0);
                        }
                        if (genotype(locus, pop, indiv, 0) == genotype(locus, other_pop, indiv_other_pop, 1))
                        {
                            ++Q2.at(0);
                        }
                        if (genotype(locus, pop, indiv, 1) == genotype(locus, other_pop, indiv_other_pop, 0))
                        {
                            ++Q2.at(0);
                        }
                        if (genotype(locus, pop, indiv, 1) == genotype(locus, other_pop, indiv_other_pop, 1))
                        {
                            ++Q2.at(0);
                        }
                        Q2.at(1) += 4;
                    }
                }
            }
        }
    }
    return static_cast<double>(Q2.at(0)) / Q2.at(1);
}

double calc_Fstat(double Qx, double Qy)
{
    return (Qx - Qy) / (1 - Qy);
}

result_qstat calc_qstat(data_plane_vec_c const &genotype)
{
    auto Q0_intra_ind = std::array<int, 2>{0, 0};
    auto Q1_intra_pop = std::array<int, 2>{0, 0};
    auto Q2_inter_pop = std::array<int, 2>{0, 0};

    for (int locus = 0; locus < genotype.locus_nbr(); ++locus)
    {
        int locus_end = genotype.index_end_locus(locus);
        for (int gene1 = genotype.index_begin_locus(locus); gene1 < locus_end - 1; ++gene1)
        {
            for (int gene2 = gene1 + 1; gene2 < locus_end; ++gene2)
            {
                if (genotype.same_indiv(gene1, gene2))
                {
                    if (genotype[gene1] == genotype[gene2])
                    {
                        ++Q0_intra_ind.at(0);
                    }
                    ++Q0_intra_ind.at(1);
                }

                if (genotype.pop_at_dist(gene1, gene2, 0))
                {
                    if (genotype[gene1] == genotype[gene2])
                    {
                        ++Q1_intra_pop.at(0);
                    }
                    ++Q1_intra_pop.at(1);
                }
                else
                {
                    if (genotype[gene1] == genotype[gene2])
                    {
                        ++Q2_inter_pop.at(0);
                    }
                    ++Q2_inter_pop.at(1);
                }
            }
        }
    }

    result_qstat result;
    result.Q0_intra_ind = static_cast<double>(Q0_intra_ind.at(0)) / Q0_intra_ind.at(1);
    result.Q1_intra_pop = static_cast<double>(Q1_intra_pop.at(0)) / Q1_intra_pop.at(1);
    result.Q2_inter_pop = static_cast<double>(Q2_inter_pop.at(0)) / Q2_inter_pop.at(1);
    return result;
}

std::vector<double> calc_qr(data_plane_vec_c const &genotype, int dist_max)
{
    std::vector<std::array<int, 2>> result_fract(dist_max + 1);
    for (int locus = 0; locus < genotype.locus_nbr(); ++locus)
    {
        int locus_end = genotype.index_end_locus(locus);
        for (int gene1 = genotype.index_begin_locus(locus); gene1 < locus_end - 1; ++gene1)
        {
            for (int gene2 = gene1 + 1; gene2 < locus_end; ++gene2)
            {
                for (int dist = 0; dist <= dist_max; ++dist)
                {
                    if (genotype.pop_at_dist(gene1, gene2, dist))
                    {
                        if (genotype[gene1] == genotype[gene2])
                        {
                            ++result_fract[dist].at(0);
                        }
                        ++result_fract[dist].at(1);
                    }
                }
            }
        }
    }

    std::vector<double> result(result_fract.size());
    auto result_itr = result.begin();
    for (auto const &frac : result_fract)
    {
        *result_itr = static_cast<double>(frac.at(0)) / frac.at(1);
        ++result_itr;
    }
    return result;
}

template <typename value>
double mean(std::vector<value> X)
{
    double result;
    return result;
}

template <typename value>
double var(std::vector<value> X)
{
    double result;
    return result;
}

template <typename value>
double cov_X_Y(std::vector<value> X, std::vector<value> Y)
{
    double result;
    return result;
}