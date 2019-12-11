#include <tuple>
#include <cmath>

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

result_qstat calc_qstat_loc_by_loc(data_plane_vec_c const &genotype, int locus)
{
    auto Q0_intra_ind = std::array<int, 2>{0, 0};
    auto Q1_intra_pop = std::array<int, 2>{0, 0};
    auto Q2_inter_pop = std::array<int, 2>{0, 0};

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

    result_qstat result;
    result.Q0_intra_ind = static_cast<double>(Q0_intra_ind.at(0)) / Q0_intra_ind.at(1);
    result.Q1_intra_pop = static_cast<double>(Q1_intra_pop.at(0)) / Q1_intra_pop.at(1);
    result.Q2_inter_pop = static_cast<double>(Q2_inter_pop.at(0)) / Q2_inter_pop.at(1);
    return result;
}

//TODO : A optimiser
result_qstat calc_qstat_all_loc(data_plane_vec_c const &genotype)
{
    result_qstat temp;
    result_qstat result;

    for (int locus = 0; locus < genotype.locus_nbr(); ++locus)
    {
        temp = calc_qstat_loc_by_loc(genotype, locus);
        result.Q0_intra_ind += temp.Q0_intra_ind;
        result.Q1_intra_pop += temp.Q1_intra_pop;
        result.Q2_inter_pop += temp.Q2_inter_pop;
    }

    result.Q0_intra_ind /= genotype.locus_nbr();
    result.Q1_intra_pop /= genotype.locus_nbr();
    result.Q2_inter_pop /= genotype.locus_nbr();

    return result;
}

std::vector<double> calc_qr_all_loc(data_plane_vec_c const &genotype, int dist_max)
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

std::vector<double> Fst_qstat_all_loc(std::vector<double> const &calc_qr_all_loc)
{
    auto Q0 = calc_qr_all_loc[0];
    std::vector<double> result(calc_qr_all_loc.size());
    auto result_itr = result.begin();

    for (auto const &value : calc_qr_all_loc)
    {
        *result_itr = (Q0 - value) / (1 - value);
        ++result_itr;
    }
    return result;
}

template <typename value>
double mean(std::vector<value> X_vec)
{
    double result{0};
    for (value const &val : X_vec)
    {
        result += val;
    }

    return result / X_vec.size();
}

template <typename value>
double var(std::vector<value> X_vec, double meanX)
{
    double result{0};
    for (value const &val : X_vec)
    {
        result += std::pow(val - meanX, 2);
    }
    return result / X_vec.size();
}

template <typename value1, typename value2>
double cov_X_Y(std::vector<value1> X_vec, double meanX, std::vector<value2> Y_vec, double meanY)
{
    if (X_vec.size() != Y_vec.size())
    {
        throw std::logic_error("Vector don't have the same size");
    }

    double result{0};
    auto Y_itr = Y_vec.begin();

    for (value1 const &val : X_vec)
    {
        result += (val - meanX) * (*Y_itr - meanY);
        ++Y_itr;
    }

    return result / X_vec.size();
}

template <typename value1, typename value2>
std::array<double, 2> linear_regres_X_Y(std::vector<value1> X_vec, std::vector<value2> Y_vec)
{
    double meanX = mean(X_vec);
    double meanY = mean(Y_vec);
    double a = cov_X_Y(X_vec, meanX, Y_vec, meanY) / var(X_vec, meanX);
    double b = meanY - meanX * a;

    return {a, b};
}

//Explicit instantiation
template double mean(std::vector<int> X);
template double mean(std::vector<double> X);

template double var(std::vector<int> X, double meanX);
template double var(std::vector<double> X, double meanX);

template double cov_X_Y(std::vector<int> X, double meanX, std::vector<int> Y, double meanY);
template double cov_X_Y(std::vector<int> X, double meanX, std::vector<double> Y, double meanY);
template double cov_X_Y(std::vector<double> X, double meanX, std::vector<int> Y, double meanY);
template double cov_X_Y(std::vector<double> X, double meanX, std::vector<double> Y, double meanY);

template std::array<double, 2> linear_regres_X_Y(std::vector<int> X, std::vector<int> Y);
template std::array<double, 2> linear_regres_X_Y(std::vector<int> X, std::vector<double> Y);
template std::array<double, 2> linear_regres_X_Y(std::vector<double> X, std::vector<int> Y);
template std::array<double, 2> linear_regres_X_Y(std::vector<double> X, std::vector<double> Y);