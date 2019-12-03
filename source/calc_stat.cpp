#include <tuple>

#include "calc_stat.hpp"

std::array<int, 2> calc_Q0(data_plane_vec_c const &genotype)
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
    return Q0;
}

//pop/locus/indiv
std::array<int, 2> calc_Q1(data_plane_vec_c const &genotype)
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
    return Q1;
}

//pop/locus/indiv
std::array<int, 2> calc_Q2(data_plane_vec_c const &genotype)
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
    return Q2;
}

std::array<int, 2> calc_Fst(std::array<int, 2> Q1, std::array<int, 2> Q2)
{
    std::array<int, 2> fst;
    //Q1 - Q2 / 1 - Q2 => {a/b - c/d} / {1 - c/d} => {ad-cb}/{bd-cb}
    fst.at(0) = Q1.at(0) * Q2.at(1) - Q2.at(0) * Q1.at(1);
    fst.at(1) = Q1.at(1) * Q2.at(1) - Q2.at(0) * Q1.at(1);
    return fst;
}

result_qstat calc_qstat(data_plane_vec_c const &genotype)
{
    result_qstat result;

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
                        ++result.Q0.at(0);
                    }
                    ++result.Q0.at(1);
                }

                if (genotype.same_pop(locus, gene1, gene2))
                {
                    if (genotype[gene1] == genotype[gene2])
                    {
                        ++result.Q1.at(0);
                    }
                    ++result.Q1.at(1);
                }
                else
                {
                    if (genotype[gene1] == genotype[gene2])
                    {
                        ++result.Q2.at(0);
                    }
                    ++result.Q2.at(1);
                }
            }
        }
    }
    return result;
}