#include <tuple>
#include <map>
#include <cmath>

#include "calc_stat.hpp"

//calc_Q_intra_indiv => calc_Qwi_frac
double calc_Q_intra_indiv(data_plane_vec_c const &data_plane_vec)
{
    auto Q0 = std::array<int, 2>{0, 0};
    for (int i = 0; i < data_plane_vec.size(); i += 2)
    {
        if ((data_plane_vec[i] != 0) && (data_plane_vec[i + 1] != 0))
        {
            if (data_plane_vec[i] == data_plane_vec[i + 1])
            {
                ++Q0.at(0);
            }
            ++Q0.at(1);
        }
    }
    return static_cast<double>(Q0.at(0)) / Q0.at(1);
}

//pop/locus/indiv
double calc_Q_inter_indiv_intra_pop(data_plane_vec_c const &data_plane_vec)
{
    auto Q1 = std::array<int, 2>{0, 0};
    for (int locus = 0; locus < data_plane_vec.nbr_of_locus(); ++locus)
    {
        for (int pop = 0; pop < data_plane_vec.nbr_of_pop(); ++pop)
        {
            //Handle last indiv
            for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv_per_pop(pop) - 1; ++indiv)
            {
                int indiv1_gene1 = data_plane_vec(locus, pop, indiv, 0);
                int indiv1_gene2 = data_plane_vec(locus, pop, indiv, 1);

                for (int next_indiv = indiv + 1; next_indiv < data_plane_vec.nbr_of_indiv_per_pop(pop); ++next_indiv)
                {
                    int indiv2_gene1 = data_plane_vec(locus, pop, next_indiv, 0);
                    int indiv2_gene2 = data_plane_vec(locus, pop, next_indiv, 1);
                    if ((indiv1_gene1 != 0))
                    {
                        if ((indiv2_gene1 != 0))
                        {
                            if (indiv1_gene1 == indiv2_gene1)
                            {
                                ++Q1.at(0);
                            }
                            ++Q1.at(1);
                        }
                        if ((indiv2_gene2 != 0))
                        {
                            if (indiv1_gene1 == indiv2_gene2)
                            {
                                ++Q1.at(0);
                            }
                            ++Q1.at(1);
                        }
                    }
                    if ((indiv1_gene2 != 0))
                    {
                        if ((indiv2_gene1 != 0))
                        {
                            if (indiv1_gene2 == indiv2_gene1)
                            {
                                ++Q1.at(0);
                            }
                            ++Q1.at(1);
                        }
                        if ((indiv2_gene2 != 0))
                        {
                            if (indiv1_gene2 == indiv2_gene2)
                            {
                                ++Q1.at(0);
                            }
                            ++Q1.at(1);
                        }
                    }
                }
            }
        }
    }
    return static_cast<double>(Q1.at(0)) / Q1.at(1);
}

//pop/locus/indiv
double calc_Q_inter_pop(data_plane_vec_c const &data_plane_vec)
{
    auto Q2 = std::array<int, 2>{0, 0};
    for (int locus = 0; locus < data_plane_vec.nbr_of_locus(); ++locus)
    {
        for (int pop = 0; pop < data_plane_vec.nbr_of_pop() - 1; ++pop)
        {
            for (int other_pop = pop + 1; other_pop < data_plane_vec.nbr_of_pop(); ++other_pop)
            {
                for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv_per_pop(pop); ++indiv)
                {
                    int indiv1_gene1 = data_plane_vec(locus, pop, indiv, 0);
                    int indiv1_gene2 = data_plane_vec(locus, pop, indiv, 1);

                    for (int indiv_other_pop = 0; indiv_other_pop < data_plane_vec.nbr_of_indiv_per_pop(other_pop); ++indiv_other_pop)
                    {
                        int indiv2_gene1 = data_plane_vec(locus, other_pop, indiv_other_pop, 0);
                        int indiv2_gene2 = data_plane_vec(locus, other_pop, indiv_other_pop, 1);
                        if ((indiv1_gene1 != 0))
                        {
                            if ((indiv2_gene1 != 0))
                            {
                                if (indiv1_gene1 == indiv2_gene1)
                                {
                                    ++Q2.at(0);
                                }
                                ++Q2.at(1);
                            }
                            if ((indiv2_gene2 != 0))
                            {
                                if (indiv1_gene1 == indiv2_gene2)
                                {
                                    ++Q2.at(0);
                                }
                                ++Q2.at(1);
                            }
                        }
                        if ((indiv1_gene2 != 0))
                        {
                            if ((indiv2_gene1 != 0))
                            {
                                if (indiv1_gene2 == indiv2_gene1)
                                {
                                    ++Q2.at(0);
                                }
                                ++Q2.at(1);
                            }
                            if ((indiv2_gene2 != 0))
                            {
                                if (indiv1_gene2 == indiv2_gene2)
                                {
                                    ++Q2.at(0);
                                }
                                ++Q2.at(1);
                            }
                        }
                    }
                }
            }
        }
    }
    return static_cast<double>(Q2.at(0)) / Q2.at(1);
}

double calc_Hnei_per_loc(data_plane_vec_c const &data_plane_vec, int locus)
{
    double result = 0;
    for (auto freq_al : data_plane_vec.allele_state_per_loc(locus))
    {
        result += pow(freq_al.at(1), 2);
    }
    double nomiss_nbr_of_gene = data_plane_vec.nomiss_nbr_of_gene_per_loc(locus);
    result = (1 - (result / pow(nomiss_nbr_of_gene, 2))) * nomiss_nbr_of_gene / (nomiss_nbr_of_gene - 1);
    return result;
}

double calc_Hnei(data_plane_vec_c const &data_plane_vec)
{
    double result = 0;
    int nbr_of_locus = data_plane_vec.nbr_of_locus();
    for (auto locus = 0; locus < nbr_of_locus; ++locus)
    {
        result += calc_Hnei_per_loc(data_plane_vec, locus);
    }

    return result / nbr_of_locus;
}

//WARNING : For microsat only
double calc_Var_per_loc(data_plane_vec_c const &data_plane_vec, int locus)
{
    double moy = 0;
    int nbr_of_count = 0;
    for (auto freq_al : data_plane_vec.allele_state_per_loc(locus))
    {
        //freq_al.at(0) = value of microsat ; freq_al.at(1) = number of microsat
        moy += freq_al.at(0) * freq_al.at(1);
        nbr_of_count += freq_al.at(1);
    }
    moy /= nbr_of_count;
    double result = 0;
    for (auto freq_al : data_plane_vec.allele_state_per_loc(locus))
    {
        result += freq_al.at(1) * pow(freq_al.at(0) - moy, 2);
    }
    //Correction n/n-1
    double n = data_plane_vec.nomiss_nbr_of_gene_per_loc(locus);
    return (result / nbr_of_count) * (n / (n - 1));
}
//calc_Q_intra_indiv => calc_Qwi_frac
double calc_Hobs_per_loc(data_plane_vec_c const &data_plane_vec, int locus)
{
    auto Q0 = std::array<int, 2>{0, 0};

    for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv(); ++indiv)
    {
        int indiv1_gene1 = data_plane_vec(locus, indiv, 0);
        int indiv1_gene2 = data_plane_vec(locus, indiv, 1);
        if ((indiv1_gene1 != 0) && (indiv1_gene2 != 0))
        {
            if (indiv1_gene1 != indiv1_gene2)
            {
                ++Q0.at(0);
            }
            ++Q0.at(1);
        }
    }
    return static_cast<double>(Q0.at(0)) / Q0.at(1);
}

double calc_MGW_per_loc(data_plane_vec_c const &data_plane_vec, int locus)
{
    auto const &allele_state_per_loc = data_plane_vec.allele_state_per_loc(locus);
    //nbr allele / rage allele
    return static_cast<double>(allele_state_per_loc.size()) / (allele_state_per_loc[allele_state_per_loc.size() - 1].at(0) - allele_state_per_loc[0].at(0) + 1);
}

double calc_MGW(data_plane_vec_c const &data_plane_vec)
{
    double result = 0;
    int nbr_of_locus = data_plane_vec.nbr_of_locus();
    for (auto locus = 0; locus < nbr_of_locus; ++locus)
    {
        result += calc_MGW_per_loc(data_plane_vec, locus);
    }

    return result / nbr_of_locus;
}

#include <iostream>
// std::vector<qr_num, qr_denum>
std::vector<std::array<int, 2>> calc_qr_loc_by_loc(data_plane_vec_c const &data_plane_vec, int locus)
{
    std::map<double, std::array<int, 2>> result_fract;
    int locus_end = data_plane_vec.index_end_locus(locus);
    for (int gene1 = data_plane_vec.index_begin_locus(locus); gene1 < locus_end - 1; ++gene1)
    {
        for (int gene2 = gene1 + 1; gene2 < locus_end; ++gene2)
        {
            auto dist_class = data_plane_vec.dist_class_btw_pop(gene1, gene2);
            auto &frac = result_fract[dist_class];
            if (data_plane_vec[gene1] == data_plane_vec[gene2])
            {
                ++frac.at(0);
            }
            ++frac.at(1);
        }
    }

    std::vector<std::array<int, 2>> result(result_fract.size());
    auto result_itr = result.begin();
    for (auto const &frac : result_fract)
    {
        *result_itr = {frac.second.at(0), frac.second.at(1)};
        ++result_itr;
    }
    return result;
}

// std::vector<qr>
std::vector<double> calc_qr_all_loc(data_plane_vec_c const &data_plane_vec)
{
    std::vector<std::array<int, 2>> result_frac(data_plane_vec.nbr_of_dist_class());

    for (int locus = 0; locus < data_plane_vec.nbr_of_locus(); ++locus)
    {
        auto temp = calc_qr_loc_by_loc(data_plane_vec, locus);

        for (int i = 0; i < temp.size(); ++i)
        {
            result_frac[i].at(0) += temp[i].at(0);
            result_frac[i].at(1) += temp[i].at(1);
        }
    }

    std::vector<double> result(data_plane_vec.nbr_of_dist_class());
    auto result_frac_itr = result_frac.begin();
    for (auto &res : result)
    {
        res = static_cast<double>(result_frac_itr->at(0)) / result_frac_itr->at(1);
        ++result_frac_itr;
    }

    return result;
}

std::vector<std::array<double, 2>> ar_by_pair(data_plane_vec_c const &data_plane_vec)
{
    int Ploidy = data_plane_vec.get_Ploidy();
    if (Ploidy != 2)
    {
        throw std::logic_error("Can't calculate ar if ploidy was different than 2.");
    }

    int nbr_of_pair = combination(2, data_plane_vec.nbr_of_indiv());

    //numerator and denominator
    std::vector<std::array<double, 2>> result(nbr_of_pair, {0, 0});
    //In prob_id case Qw with a missing value = 0; Denom by locus will be the same for all pair
    std::vector<double> Qw_by_locus(data_plane_vec.nbr_of_locus(), 0);
    //
    double sum_Qw_all_loc = 0;

    //Multilocus estimates are defined as the sum of locus-specific numerators divided by the sum of locus-specific denominators (Weir & Cockerham, 1984)
    for (int locus = 0; locus < data_plane_vec.nbr_of_locus(); ++locus)
    {
        auto result_itr = result.begin();
        for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv() - 1; ++indiv)
        {
            int indiv1_gene1 = data_plane_vec(locus, indiv, 0);
            int indiv1_gene2 = data_plane_vec(locus, indiv, 1);
            if ((indiv1_gene1 != 0) && (indiv1_gene2 != 0))
            {
                double semi_Qw = (indiv1_gene1 == indiv1_gene2);

                for (int next_indiv = indiv + 1; next_indiv < data_plane_vec.nbr_of_indiv(); ++next_indiv)
                {
                    int indiv2_gene1 = data_plane_vec(locus, next_indiv, 0);
                    int indiv2_gene2 = data_plane_vec(locus, next_indiv, 1);
                    if ((indiv2_gene1 != 0) && (indiv2_gene2 != 0))
                    {
                        double Qw = (semi_Qw + (indiv2_gene1 == indiv2_gene2)) / 2.0;
                        double Qr = ((indiv1_gene1 == indiv2_gene1) + (indiv1_gene1 == indiv2_gene2) + (indiv1_gene2 == indiv2_gene1) + (indiv1_gene2 == indiv2_gene2)) / 4.0;

                        if (locus == 0)
                        {
                            result_itr->at(0) = data_plane_vec.dist_btw_pop(indiv * Ploidy, next_indiv * Ploidy);
                        }
                        result_itr->at(1) += Qw - 2 * Qr + 1;
                    }
                    ++result_itr;
                }
                //Use to estimate Qw for all pair
                Qw_by_locus[locus] += semi_Qw;
            }
            else
            {
                //ptr arithmetic (if indiv have missing value it's unnecessary to compute all his pair)
                result_itr = result_itr + (data_plane_vec.nbr_of_indiv() - indiv - 1);
            }
        }
        //HAndle last indiv
        Qw_by_locus[locus] += (data_plane_vec(locus, data_plane_vec.nbr_of_indiv() - 1, 0) == data_plane_vec(locus, data_plane_vec.nbr_of_indiv() - 1, 1));
        Qw_by_locus[locus] /= data_plane_vec.nomiss_nbr_of_indiv_per_loc(locus);
        sum_Qw_all_loc += Qw_by_locus[locus];
    }

    //compute denom because in missing value denom will be diff between pair (When the data_plane_vec is unknown for some locus in one individual from a pair, this locus is not included in the sum over loci Rousset, 2000)
    //compute by remove Qw from sum_Qw_all_loc
    std::vector<std::array<double, 2>> sum_Qw_use_by_pair(nbr_of_pair, {sum_Qw_all_loc, static_cast<double>(data_plane_vec.nbr_of_locus())});

    auto sum_Qw_use_by_pair_itr = sum_Qw_use_by_pair.begin();
    for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv() - 1; ++indiv)
    {
        for (int next_indiv = indiv + 1; next_indiv < data_plane_vec.nbr_of_indiv(); ++next_indiv)
        {
            std::vector<bool> nomiss_data = bin_vec::and_(data_plane_vec.nomiss_data_indiv(indiv), data_plane_vec.nomiss_data_indiv(next_indiv));
            for (int locus = 0; locus < data_plane_vec.nbr_of_locus(); ++locus)
            {
                //if indiv1 OR indiv2 have missing data at loc, remove it from Qw_sum
                if (!nomiss_data[locus])
                {
                    sum_Qw_use_by_pair_itr->at(0) -= Qw_by_locus[locus];
                    //to concerve how many locus have been use in this sum
                    --sum_Qw_use_by_pair_itr->at(1);
                }
            }
            ++sum_Qw_use_by_pair_itr;
        }
    }

    sum_Qw_use_by_pair_itr = sum_Qw_use_by_pair.begin();
    for (auto &value : result)
    {
        value.at(1) = value.at(1) / (2 * (sum_Qw_use_by_pair_itr->at(1) - sum_Qw_use_by_pair_itr->at(0))) - 0.5;
        ++sum_Qw_use_by_pair_itr;
    }

    return result;
}

std::vector<std::array<double, 2>> er_by_pair(data_plane_vec_c const &data_plane_vec)
{
    int Ploidy = data_plane_vec.get_Ploidy();
    if (Ploidy != 2)
    {
        throw std::logic_error("Can't calculate ar if ploidy was different than 2.");
    }

    int nbr_indiv = data_plane_vec.nbr_of_indiv();
    int nbr_of_pair = combination(2, nbr_indiv);

    //numerator and denominator
    std::vector<std::array<double, 2>> save_value(nbr_of_pair, {0, 0});
    //To calculate e_r like Genepop (Rousset, ...) with Loiselle F statistic equivalent d_k^2 terme
    std::vector<double> sum_Qij_by_locus(data_plane_vec.nbr_of_locus(), 0);
    //In prob_id case Qw with a missing value = 0; Denom by locus will be the same for all pair
    std::vector<double> Qw_by_locus(data_plane_vec.nbr_of_locus(), 0);
    double sum_Qw_all_loc = 0;

    //All Qi need to be concerve for later correction
    std::vector<std::vector<double>> Qi_by_indiv_by_locus(nbr_indiv, std::vector<double>(data_plane_vec.nbr_of_locus(), 0));
    std::vector<double> Qi_sum_all_loc(nbr_indiv, 0);

    //Multilocus estimates are defined as the sum of locus-specific numerators divided by the sum of locus-specific denominators (Weir & Cockerham, 1984)
    for (int locus = 0; locus < data_plane_vec.nbr_of_locus(); ++locus)
    {
        auto nomiss_indiv = data_plane_vec.nomiss_nbr_of_indiv_per_loc(locus);
        auto save_value_itr = save_value.begin();
        //Handle all indiv
        for (int indiv = 0; indiv < nbr_indiv; ++indiv)
        {
            int indiv1_gene1 = data_plane_vec(locus, indiv, 0);
            int indiv1_gene2 = data_plane_vec(locus, indiv, 1);
            if ((indiv1_gene1 != 0) && (indiv1_gene2 != 0))
            {
                double semi_Qw = (indiv1_gene1 == indiv1_gene2);
                for (int next_indiv = indiv + 1; next_indiv < nbr_indiv; ++next_indiv)
                {
                    int indiv2_gene1 = data_plane_vec(locus, next_indiv, 0);
                    int indiv2_gene2 = data_plane_vec(locus, next_indiv, 1);
                    if ((indiv2_gene1 != 0) && (indiv2_gene2 != 0))
                    {
                        double Qij = ((indiv1_gene1 == indiv2_gene1) + (indiv1_gene1 == indiv2_gene2) + (indiv1_gene2 == indiv2_gene1) + (indiv1_gene2 == indiv2_gene2)) / 4.0;

                        if (locus == 0)
                        {
                            save_value_itr->at(0) = data_plane_vec.dist_btw_pop(indiv * Ploidy, next_indiv * Ploidy);
                        }
                        save_value_itr->at(1) += Qij;
                        sum_Qij_by_locus[locus] += Qij;
                        //Qi = mean Qij between i and j = all no missing indiv
                        Qij /= nomiss_indiv;
                        //Indiv and next_indiv have to be update to calculate Qindiv and Qnext_indiv mean by loc
                        Qi_by_indiv_by_locus[indiv][locus] += Qij;
                        Qi_by_indiv_by_locus[next_indiv][locus] += Qij;
                        //compute Qi. and Qw because in missing value will be diff between pair (When the data_plane_vec is unknown for some locus in one individual from a pair, this locus is not included in the sum over loci Rousset, 2000)
                        //compute by remove Qi.[locus] from Qi_sum_by_indiv and Qw[locus] from sum_Qw_all_loc
                        Qi_sum_all_loc[indiv] += Qij;
                        Qi_sum_all_loc[next_indiv] += Qij;
                    }
                    ++save_value_itr;
                }
                //2 for auto comparaison
                double Qii = (2 + 2 * (indiv1_gene1 == indiv1_gene2)) / 4.0;
                //Handle auto comparison
                Qii /= nomiss_indiv;
                Qi_by_indiv_by_locus[indiv][locus] += Qii;
                Qi_sum_all_loc[indiv] += Qii;
                //Use to estimate Qi and Qw for all pair
                Qw_by_locus[locus] += semi_Qw;
            }
            else
            {
                //ptr arithmetic (if indiv have missing value it's unnecessary to compute all his pair)
                save_value_itr = save_value_itr + (nbr_indiv - indiv - 1);
            }
        }
        Qw_by_locus[locus] /= data_plane_vec.nomiss_nbr_of_indiv_per_loc(locus);
        sum_Qw_all_loc += Qw_by_locus[locus];
    }

    //calc of terme for match Genepop based on Loisel F statistic

    std::vector<double> Qi_use_by_pair(nbr_of_pair, 0);
    //nbr of pair use to calc sum_Qij (bijection with the number of indiv)
    std::vector<double> l_term_use_by_pair(nbr_of_pair, 0);
    //1 - sum_Qw_all_loc to match denom in cal er
    std::vector<double> sum_Qw_use_by_pair(nbr_of_pair, data_plane_vec.nbr_of_locus() - sum_Qw_all_loc);

    auto Qi_use_by_pair_itr = Qi_use_by_pair.begin();
    auto l_term_use_by_pair_itr = l_term_use_by_pair.begin();
    auto sum_Qw_use_by_pair_itr = sum_Qw_use_by_pair.begin();
    for (int indiv = 0; indiv < nbr_indiv - 1; ++indiv)
    {
        for (int next_indiv = indiv + 1; next_indiv < nbr_indiv; ++next_indiv)
        {
            //file Qi_use_by_pair with sum of Qi and Qj (if i or j have missing value, Qi or Qj = 0, just need to remove the other one)
            *Qi_use_by_pair_itr = Qi_sum_all_loc[indiv] + Qi_sum_all_loc[next_indiv];
            std::vector<bool> nomiss_data = bin_vec::and_(data_plane_vec.nomiss_data_indiv(indiv), data_plane_vec.nomiss_data_indiv(next_indiv));
            for (int locus = 0; locus < data_plane_vec.nbr_of_locus(); ++locus)
            {
                //if indiv1 OR indiv2 have missing data at loc, remove it from Qw_sum
                if (!nomiss_data[locus])
                {
                    *Qi_use_by_pair_itr -= (Qi_by_indiv_by_locus[indiv][locus] + Qi_by_indiv_by_locus[next_indiv][locus]);

                    *sum_Qw_use_by_pair_itr -= (1 - Qw_by_locus[locus]);
                }
                else
                {
                    //TODO : Essayer de faire les calculs par décroissance
                    //Loiselle F statistic equivalent d_k^2 terme
                    int nomiss_nbr_indiv = data_plane_vec.nomiss_nbr_of_indiv_per_loc(locus);
                    int nomiss_nbr_pair_indiv = combination(2, nomiss_nbr_indiv);
                    *l_term_use_by_pair_itr += (sum_Qij_by_locus[locus] + nomiss_nbr_indiv * ((1 + Qw_by_locus[locus]) / 2.0)) / (nomiss_nbr_indiv + nomiss_nbr_pair_indiv);
                }
            }
            ++Qi_use_by_pair_itr;
            ++l_term_use_by_pair_itr;
            ++sum_Qw_use_by_pair_itr;
        }
    }

    //Calc e_r
    std::vector<std::array<double, 2>> result(nbr_of_pair, {0, 0});
    auto result_itr = result.begin();
    auto save_value_itr = save_value.begin();
    Qi_use_by_pair_itr = Qi_use_by_pair.begin();
    l_term_use_by_pair_itr = l_term_use_by_pair.begin();
    sum_Qw_use_by_pair_itr = sum_Qw_use_by_pair.begin();

    for (int indiv = 0; indiv < nbr_indiv - 1; ++indiv)
    {
        for (int next_indiv = indiv + 1; next_indiv < nbr_indiv; ++next_indiv)
        {
            result_itr->at(0) = save_value_itr->at(0);
            //in paper is -e
            result_itr->at(1) = (-save_value_itr->at(1) + *Qi_use_by_pair_itr - *l_term_use_by_pair_itr) / *sum_Qw_use_by_pair_itr;
            ++result_itr;
            ++save_value_itr;
            ++l_term_use_by_pair_itr;
            ++Qi_use_by_pair_itr;
            ++sum_Qw_use_by_pair_itr;
        }
    }

    return result;
}

//esti Fis, Fst
//WARNING : Not know pertinence if missing value
std::array<std::array<double, 2>, 2> Fstat_by_loc_with_probid(data_plane_vec_c const &data_plane_vec, int locus)
{
    std::array<std::array<double, 2>, 2> result{{{0}, {0}}};

    int Ploidy = data_plane_vec.get_Ploidy();
    if (Ploidy != 2)
    {
        throw std::logic_error("in Fstat_by_loc_with_probid : Can't calculate Fstat if ploidy was different than 2.");
    }
    //
    auto Qwi_frac = std::array<int, 2>{0, 0};
    //+
    auto Qwp_frac_v = std::vector<std::array<int, 2>>(data_plane_vec.nbr_of_pop(), {0, 0});
    auto Qbp_frac = std::array<int, 2>{0, 0};

    int locus_end = data_plane_vec.index_end_locus(locus);
    auto Qwp_frac_v_itr = Qwp_frac_v.begin();
    auto cumul_pop_size_itr = data_plane_vec.cumul_nbr_of_indiv_per_pop().begin() + 1;
    int indiv_relatif_nbr = -1;

    for (int gene1 = data_plane_vec.index_begin_locus(locus); gene1 < locus_end - 1; ++gene1)
    {
        ++indiv_relatif_nbr;
        if (indiv_relatif_nbr == (*cumul_pop_size_itr) * data_plane_vec.get_Ploidy())
        {
            ++cumul_pop_size_itr;
            ++Qwp_frac_v_itr;
        }
        if (data_plane_vec[gene1] != 0)
        {
            for (int gene2 = gene1 + 1; gene2 < locus_end; ++gene2)
            {
                if (data_plane_vec[gene2] != 0)
                {
                    if (data_plane_vec.same_indiv(gene1, gene2))
                    {
                        if (data_plane_vec[gene1] == data_plane_vec[gene2])
                        {
                            ++Qwi_frac.at(0);
                        }
                        ++Qwi_frac.at(1);
                    }
                    else
                    {
                        if (data_plane_vec.same_pop(gene1, gene2))
                        {
                            if (data_plane_vec[gene1] == data_plane_vec[gene2])
                            {
                                ++(Qwp_frac_v_itr->at(0));
                            }
                            ++(Qwp_frac_v_itr->at(1));
                        }
                        else
                        {
                            if (data_plane_vec[gene1] == data_plane_vec[gene2])
                            {
                                ++Qbp_frac.at(0);
                            }
                            ++Qbp_frac.at(1);
                        }
                    }
                }
            }
        }
    }

    double Qwi = static_cast<double>(Qwi_frac.at(0)) / Qwi_frac.at(1);
    double num_Qwp = 0;
    double denum_Qwp = 0;
    for (auto num_denum_Qwp : Qwp_frac_v)
    {
        num_Qwp += num_denum_Qwp.at(0);
        denum_Qwp += num_denum_Qwp.at(1);
    }
    double Qwp_g = num_Qwp / denum_Qwp;
    double Qbp = static_cast<double>(Qbp_frac.at(0)) / Qbp_frac.at(1);

    double S1 = 0;
    double S2 = 0;
    double ns = data_plane_vec.nomiss_nbr_of_pop_per_loc(locus);

    for (auto pop_size : data_plane_vec.nomiss_nbr_of_indiv_per_loc_per_pop(locus))
    {
        S1 += pop_size;
        S2 += pow(pop_size, 2);
    }
    double nc = (S1 - S2 / S1) / (ns - 1);

    double MSG = 1 - Qwi;
    double MSI = 1 + Qwi - 2 * Qwp_g;
    double MSP = 2 * nc * (Qwp_g - Qbp) + MSI;

    //Fis
    result.at(0).at(0) = (MSI - MSG);
    result.at(0).at(1) = (MSI + MSG);
    //Fst
    result.at(1).at(0) = (MSP - MSI);
    result.at(1).at(1) = (MSP + (nc - 1) * MSI + nc * MSG);

    return result;
}

//esti Fis, Fst
std::array<std::array<double, 2>, 2> Fstat_by_loc_with_indic(data_plane_vec_c const &data_plane_vec, int locus)
{
    std::vector<std::array<int, 2>> const &allele_state = data_plane_vec.allele_state_per_loc(locus);

    int Ploidy = data_plane_vec.get_Ploidy();
    if (Ploidy != 2)
    {
        throw std::logic_error("in Fstat_by_loc_with_indic : Can't calculate Fstat if ploidy was different than 2.");
    }

    double S1 = 0;
    double S2 = 0;
    double ns = data_plane_vec.nomiss_nbr_of_pop_per_loc(locus);

    for (auto pop_size : data_plane_vec.nomiss_nbr_of_indiv_per_loc_per_pop(locus))
    {
        S1 += pop_size;
        S2 += pow(pop_size, 2);
    }
    //Indicatrix matrix mean for each pop
    std::vector<std::vector<double>> pop_state_mean(data_plane_vec.nbr_of_pop(), std::vector<double>(allele_state.size(), 0));

    //Calc mean for sample
    std::vector<double> state_mean_samp(allele_state.size(), 0);
    for (int pop = 0; pop < data_plane_vec.nbr_of_pop(); ++pop)
    {
        //for handle missing value need to remove all indiv with missing value
        int gene_in_nomiss_nbr_of_indiv_per_pop = data_plane_vec.nomiss_nbr_of_indiv_per_loc_per_pop(locus, pop) * 2;
        int tot_gene_in_nomiss_nbr_of_indiv = (S1 * 2);
        for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv_per_pop(pop); ++indiv)
        {
            int indiv1_gene1 = data_plane_vec(locus, pop, indiv, 0);
            int indiv1_gene2 = data_plane_vec(locus, pop, indiv, 1);

            if ((indiv1_gene1 != 0) && (indiv1_gene2 != 0))
            {
                int count = 0;
                for (auto const &state : allele_state)
                {
                    double temp = (static_cast<double>(indiv1_gene1 == state.at(0)) + static_cast<double>(indiv1_gene2 == state.at(0)));
                    pop_state_mean[pop][count] += (temp / gene_in_nomiss_nbr_of_indiv_per_pop);
                    state_mean_samp[count] += (temp / tot_gene_in_nomiss_nbr_of_indiv);
                    ++count;
                }
            }
        }
    }

    double SSg = 0, SSi = 0, SSp = 0;

    for (int pop = 0; pop < data_plane_vec.nbr_of_pop(); ++pop)
    {
        for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv_per_pop(pop); ++indiv)
        {
            int indiv1_gene1 = data_plane_vec(locus, pop, indiv, 0);
            int indiv1_gene2 = data_plane_vec(locus, pop, indiv, 1);

            if ((indiv1_gene1 != 0) && (indiv1_gene2 != 0))
            {
                int count = 0;
                for (auto const &state : allele_state)
                {
                    int Xij1 = (indiv1_gene1 == state.at(0));
                    int Xij2 = (indiv1_gene2 == state.at(0));
                    double Xij = (Xij1 + Xij2) / 2.0;
                    SSg += pow(Xij1 - Xij, 2) + pow(Xij2 - Xij, 2);
                    SSi += (pow(Xij - pop_state_mean[pop][count], 2)) * 2;
                    //Two genes
                    SSp += (pow(pop_state_mean[pop][count] - state_mean_samp[count], 2) * 2);
                    ++count;
                }
            }
        }
    }

    double nc = (S1 - S2 / S1) / (ns - 1);

    double MSG = SSg / S1;
    double MSI = SSi / (S1 - ns);
    double MSP = SSp / (ns - 1);

    //Fis, Fst
    std::array<std::array<double, 2>, 2> result{{{0}, {0}}};
    //Fis
    result.at(0).at(0) = nc * (MSI - MSG);
    result.at(0).at(1) = nc * (MSI + MSG);
    //Fst
    result.at(1).at(0) = (MSP - MSI);
    result.at(1).at(1) = (MSP + (nc - 1) * MSI + nc * MSG);

    return result;
}

std::array<double, 2> Fstat_genepop(data_plane_vec_c const &data_plane_vec)
{
    std::array<double, 2> result{0, 0};

    double fis_num = 0, fis_denum = 0, fst_num = 0, fst_denum = 0;
    for (int locus = 0; locus < data_plane_vec.nbr_of_locus(); ++locus)
    {
        auto temp = Fstat_by_loc_with_indic(data_plane_vec, locus);
        fis_num += temp.at(0).at(0);
        fis_denum += temp.at(0).at(1);
        fst_num += temp.at(1).at(0);
        fst_denum += temp.at(1).at(1);
    }

    result.at(0) = fis_num / fis_denum;
    result.at(1) = fst_num / fst_denum;
    return result;
}
