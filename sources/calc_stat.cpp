#include <tuple>
#include <map>
#include <cmath>
#include <random>

#include "calc_stat.hpp"

//calc_Q_intra_indiv => calc_Qwi_frac
double calc_Q_intra_indiv_per_chr_per_locus(data_plane_vec_c const &data_plane_vec, int chr, int locus)
{
    auto Q0 = std::array<int, 2>{0, 0};
    for (int indiv = 0; indiv < data_plane_vec.nomiss_nbr_of_indiv(chr, locus); ++indiv)
    {
        int indiv1_gene1 = data_plane_vec(chr, locus, indiv, 0);
        int indiv1_gene2 = data_plane_vec(chr, locus, indiv, 1);
        if ((indiv1_gene1 != 0) && (indiv1_gene2 != 0))
        {
            if (indiv1_gene1 == indiv1_gene2)
            {
                ++Q0.at(0);
            }
            ++Q0.at(1);
        }
    }
    return static_cast<double>(Q0.at(0)) / Q0.at(1);
}

double calc_Q_intra_indiv(data_plane_vec_c const &data_plane_vec)
{
    double Q0 = 0;
    for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
    {
        for (int locus = 0; locus < data_plane_vec.nbr_locus(chr); ++locus)
        {
            Q0 += calc_Q_intra_indiv_per_chr_per_locus(data_plane_vec, chr, locus);
        }
    }
    return Q0 / data_plane_vec.nbr_locus();
}

double calc_Q_inter_indiv_per_chr_per_locus_per_deme(data_plane_vec_c const &data_plane_vec, int chr, int locus, int deme)
{
    auto Q1 = std::array<int, 2>{0, 0};

    if (data_plane_vec.get_Ploidy() == 1)
    {
        for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv(deme) - 1; ++indiv)
        {
            int indiv1_gene = data_plane_vec(chr, locus, deme, indiv, 0);

            for (int next_indiv = indiv + 1; next_indiv < data_plane_vec.nbr_of_indiv(deme); ++next_indiv)
            {
                int indiv2_gene = data_plane_vec(chr, locus, deme, next_indiv, 0);
                if ((indiv1_gene != 0) && (indiv2_gene != 0))
                {
                    if (indiv1_gene == indiv2_gene)
                    {
                        ++Q1.at(0);
                    }
                    ++Q1.at(1);
                }
            }
        }
    }
    if (data_plane_vec.get_Ploidy() == 2)
    {
        for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv(deme) - 1; ++indiv)
        {
            int indiv1_gene1 = data_plane_vec(chr, locus, deme, indiv, 0);
            int indiv1_gene2 = data_plane_vec(chr, locus, deme, indiv, 1);

            for (int next_indiv = indiv + 1; next_indiv < data_plane_vec.nbr_of_indiv(deme); ++next_indiv)
            {
                int indiv2_gene1 = data_plane_vec(chr, locus, deme, next_indiv, 0);
                int indiv2_gene2 = data_plane_vec(chr, locus, deme, next_indiv, 1);
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
    return static_cast<double>(Q1.at(0)) / Q1.at(1);
}

double calc_Q_inter_indiv_per_chr_per_locus(data_plane_vec_c const &data_plane_vec, int chr, int locus)
{
    double Q1 = 0;

    for (int deme = 0; deme < data_plane_vec.nbr_of_deme(); ++deme)
    {
        Q1 += calc_Q_inter_indiv_per_chr_per_locus_per_deme(data_plane_vec, chr, locus, deme);
    }
    return Q1 / data_plane_vec.nbr_of_deme();
}

double calc_Q_inter_indiv_intra_deme(data_plane_vec_c const &data_plane_vec)
{
    double Q1 = 0;

    for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
    {
        for (int locus = 0; locus < data_plane_vec.nbr_locus(chr); ++locus)
        {
            Q1 += calc_Q_inter_indiv_per_chr_per_locus(data_plane_vec, chr, locus);
        }
    }
    return Q1 / data_plane_vec.nbr_locus();
}

double calc_Q_inter_deme_per_chr_per_locus(data_plane_vec_c const &data_plane_vec, int chr, int locus)
{
    auto Q2 = std::array<int, 2>{0, 0};

    if (data_plane_vec.get_Ploidy() == 1)
    {
        for (int deme = 0; deme < data_plane_vec.nbr_of_deme() - 1; ++deme)
        {
            for (int other_deme = deme + 1; other_deme < data_plane_vec.nbr_of_deme(); ++other_deme)
            {
                for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv(deme); ++indiv)
                {
                    int indiv1_gene = data_plane_vec(chr, locus, deme, indiv, 0);

                    for (int indiv_other_deme = 0; indiv_other_deme < data_plane_vec.nbr_of_indiv(other_deme); ++indiv_other_deme)
                    {
                        int indiv2_gene = data_plane_vec(chr, locus, other_deme, indiv_other_deme, 0);
                        if ((indiv1_gene != 0) && (indiv2_gene != 0))
                        {
                            if (indiv1_gene == indiv2_gene)
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

    if (data_plane_vec.get_Ploidy() == 2)
    {
        for (int deme = 0; deme < data_plane_vec.nbr_of_deme() - 1; ++deme)
        {
            for (int other_deme = deme + 1; other_deme < data_plane_vec.nbr_of_deme(); ++other_deme)
            {
                for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv(deme); ++indiv)
                {
                    int indiv1_gene1 = data_plane_vec(chr, locus, deme, indiv, 0);
                    int indiv1_gene2 = data_plane_vec(chr, locus, deme, indiv, 1);

                    for (int indiv_other_deme = 0; indiv_other_deme < data_plane_vec.nbr_of_indiv(other_deme); ++indiv_other_deme)
                    {
                        int indiv2_gene1 = data_plane_vec(chr, locus, other_deme, indiv_other_deme, 0);
                        int indiv2_gene2 = data_plane_vec(chr, locus, other_deme, indiv_other_deme, 1);
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

//deme/locus/indiv
double calc_Q_inter_deme(data_plane_vec_c const &data_plane_vec)
{
    double Q2 = 0;
    for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
    {
        for (int locus = 0; locus < data_plane_vec.nbr_locus(chr); ++locus)
        {
            Q2 += calc_Q_inter_deme_per_chr_per_locus(data_plane_vec, chr, locus);
        }
    }
    return Q2 / data_plane_vec.nbr_locus();
}

double calc_Hnei_per_chr_per_loc(data_plane_vec_c const &data_plane_vec, int chr, int locus)
{
    double result = 0;
    for (auto freq_al : data_plane_vec.allele_state(chr, locus))
    {
        result += pow(freq_al.second, 2);
    }
    double nomiss_nbr_of_gene = data_plane_vec.nomiss_nbr_of_gene(chr, locus);
    result = (1 - (result / pow(nomiss_nbr_of_gene, 2))) * nomiss_nbr_of_gene / (nomiss_nbr_of_gene - 1);
    return result;
}

double calc_Hnei(data_plane_vec_c const &data_plane_vec)
{
    double result = 0;
    for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
    {
        for (int locus = 0; locus < data_plane_vec.nbr_locus(chr); ++locus)
        {
            result += calc_Hnei_per_chr_per_loc(data_plane_vec, chr, locus);
        }
    }

    return result / data_plane_vec.nbr_locus();
}

//WARNING : For microsat only
double calc_Var_per_chr_per_loc(data_plane_vec_c const &data_plane_vec, int chr, int locus)
{
    double moy = 0;
    int nbr_of_count = 0;
    for (auto freq_al : data_plane_vec.allele_state(chr, locus))
    {
        //freq_al.first = state ; freq_al.second = number of allele
        moy += freq_al.first * freq_al.second;
        nbr_of_count += freq_al.second;
    }
    moy /= nbr_of_count;
    double result = 0;
    for (auto freq_al : data_plane_vec.allele_state(chr, locus))
    {
        result += freq_al.second * pow(freq_al.first - moy, 2);
    }
    //Correction n/n-1
    double n = data_plane_vec.nomiss_nbr_of_gene(chr, locus);
    return (result / nbr_of_count) * (n / (n - 1));
}

//calc_Q_intra_indiv => calc_Qwi_frac
double calc_Hobs_per_chr_per_loc(data_plane_vec_c const &data_plane_vec, int chr, int locus)
{
    auto Q0 = std::array<int, 2>{0, 0};

    for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv(); ++indiv)
    {
        int indiv1_gene1 = data_plane_vec(chr, locus, indiv, 0);
        int indiv1_gene2 = data_plane_vec(chr, locus, indiv, 1);
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

double calc_Hobs(data_plane_vec_c const &data_plane_vec)
{
    double Hobs = 0;

    for (std::size_t chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
    {
        for (std::size_t loc = 0; loc < data_plane_vec.nbr_locus(chr); ++loc)
        {
            Hobs += calc_Hobs_per_chr_per_loc(data_plane_vec, chr, loc);
        }
    }

    return Hobs / data_plane_vec.nbr_locus();
}

double calc_MGW_per_chr_per_loc(data_plane_vec_c const &data_plane_vec, int chr, int locus)
{
    auto const &allele_state = data_plane_vec.allele_state(chr, locus);
    //nbr allele / allele range (last of the allele state map - first of the allele state map)
    return static_cast<double>(allele_state.size()) / (allele_state.rbegin()->first - allele_state.begin()->first + 1);
}

double calc_MGW(data_plane_vec_c const &data_plane_vec)
{
    double result = 0;
    int nbr_value = data_plane_vec.nbr_locus();
    for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
    {
        for (int locus = 0; locus < data_plane_vec.nbr_locus(chr); ++locus)
        {
            result += calc_MGW_per_chr_per_loc(data_plane_vec, chr, locus);
        }
    }

    return result / nbr_value;
}

#include <iostream>
// std::vector<qr_num, qr_denum>
std::vector<double> calc_qr_per_chr_by_loc(data_plane_vec_c const &data_plane_vec, int chr, int locus)
{
    std::map<double, std::array<int, 2>> result_fract;
    int ploidy = data_plane_vec.get_Ploidy();
    int locus_end = data_plane_vec.index_end_locus(chr, locus);
    for (int gene1 = data_plane_vec.index_begin_locus(chr, locus); gene1 < locus_end - 1; ++gene1)
    {
        // gene1 + ploidy - (gene1 % ploidy) haploid case = gene1 + 1
        //gene1 + ploidy - (gene1 % ploidy) diploid case = gene1 + 1 (impair) or + 2 (pair)
        for (int gene2 = gene1 + ploidy - (gene1 % ploidy); gene2 < locus_end; ++gene2)
        {
            auto dist_class = data_plane_vec.geo_dist_class_btw_gene(gene1, gene2);
            auto &frac = result_fract[dist_class];
            if (data_plane_vec[gene1] == data_plane_vec[gene2])
            {
                ++frac.at(0);
            }
            ++frac.at(1);
        }
    }

    int dist_class_nbr = data_plane_vec.nbr_geo_dist_class();
    std::vector<double> result(dist_class_nbr);
    for (int i = 0; i < dist_class_nbr; ++i)
    {
        result[i] = static_cast<double>(result_fract[i].at(0)) / result_fract[i].at(1);
    }
    return result;
}

// std::vector<qr>
std::vector<double> calc_qr(data_plane_vec_c const &data_plane_vec)
{
    std::vector<double> result_frac(data_plane_vec.nbr_geo_dist_class());

    for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
    {
        for (int locus = 0; locus < data_plane_vec.nbr_locus(chr); ++locus)
        {
            auto temp = calc_qr_per_chr_by_loc(data_plane_vec, chr, locus);

            for (std::size_t i = 0; i < temp.size(); ++i)
            {
                result_frac[i] += temp[i];
            }
        }
    }

    std::vector<double> result(data_plane_vec.nbr_geo_dist_class());
    auto result_frac_itr = result_frac.begin();
    for (auto &res : result)
    {
        res = *result_frac_itr / data_plane_vec.nbr_locus();
        ++result_frac_itr;
    }

    return result;
}

// Fst /(1-Fst) = (Qr - Qid) / (1-Qid) => see Genepop doc
std::vector<std::array<double, 2>> lin_Fst_by_deme_pair(data_plane_vec_c const &data_plane_vec)
{
    int deme_pair_nbr = (data_plane_vec.nbr_of_deme() * (data_plane_vec.nbr_of_deme() - 1)) / 2;

    std::vector<double> Qr(deme_pair_nbr);
    auto Qr_iter = Qr.begin();

    std::vector<double> deme_dist(deme_pair_nbr);
    auto deme_dist_iter = deme_dist.begin();

    double Q1 = calc_Q_inter_indiv_intra_deme(data_plane_vec);
    std::vector<std::array<double, 2>> result(deme_pair_nbr);

    int ploidy = data_plane_vec.get_Ploidy();

    //lin_Fst for all deme pair
    for (int deme1 = 0; deme1 < data_plane_vec.nbr_of_deme(); ++deme1)
    {
        for (int deme2 = deme1 + 1; deme2 < data_plane_vec.nbr_of_deme(); ++deme2)
        {
            *deme_dist_iter = data_plane_vec.geo_dist_btw_deme(deme1, deme2);

            for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
            {
                for (int locus = 0; locus < data_plane_vec.nbr_locus(chr); ++chr)
                {
                    std::array<int, 2> Qr_loc = {0,0};
                    for (int indiv1 = 0; indiv1 < data_plane_vec.nbr_of_indiv(deme1); ++indiv1)
                    {
                        int indiv1_gene1 = data_plane_vec(chr, locus, deme1, indiv1, 0);
                        int indiv1_gene2 = 0;
                        if (ploidy == 2)
                        {
                            indiv1_gene2 = data_plane_vec(chr, locus, deme1, indiv1, 1);
                        }

                        if (indiv1_gene1 != 0 || indiv1_gene2 != 0)
                        {
                            for (int indiv2 = 0; indiv2 < data_plane_vec.nbr_of_indiv(deme2); ++indiv2)
                            {
                                int indiv2_gene1 = data_plane_vec(chr, locus, deme2, indiv2, 0);
                                int indiv2_gene2 = 0;
                                if (ploidy == 2)
                                {
                                    indiv2_gene2 = data_plane_vec(chr, locus, deme2, indiv2, 1);
                                }

                                if (indiv1_gene1 != 0)
                                {
                                    if (indiv2_gene1 != 0)
                                    {
                                        if (indiv1_gene1 == indiv2_gene1)
                                        {
                                            ++(Qr_loc.at(0));
                                        }
                                        ++(Qr_loc.at(1));
                                    }

                                    if (indiv2_gene2 != 0)
                                    {
                                        if (indiv1_gene1 == indiv2_gene2)
                                        {
                                            ++(Qr_loc.at(0));
                                        }
                                        ++(Qr_loc.at(1));
                                    }
                                }
                                if (indiv1_gene2 != 0)
                                {
                                    if (indiv2_gene1 != 0)
                                    {
                                        if (indiv1_gene2 == indiv2_gene1)
                                        {
                                            ++(Qr_loc.at(0));
                                        }
                                        ++(Qr_loc.at(1));
                                    }

                                    if (indiv2_gene2 != 0)
                                    {
                                        if (indiv1_gene2 == indiv2_gene2)
                                        {
                                            ++(Qr_loc.at(0));
                                        }
                                        ++(Qr_loc.at(1));
                                    }
                                }
                            }
                        }
                    }
                    *Qr_iter += (static_cast<double>(Qr_loc.at(0)) / Qr_loc.at(1));
                }
            }
            //Mean of Qr by loc
            *Qr_iter /= data_plane_vec.nbr_locus();
            ++deme_dist_iter;
            ++Qr_iter;
        }
    }

    auto result_iter = result.begin();
    Qr_iter = Qr.begin();

    for (auto &dist : deme_dist)
    {
        result_iter->at(0) = dist;
        result_iter->at(1) = (*Qr_iter - Q1) / (1 - Q1);
        ++Qr_iter;
        ++result_iter;
    }

    return result;
}

std::vector<std::array<double, 2>> ar_by_indiv_pair(data_plane_vec_c const &data_plane_vec)
{
    int Ploidy = data_plane_vec.get_Ploidy();
    if (Ploidy != 2)
    {
        throw std::logic_error("( Can't calculate ar if ploidy was different than 2. I exit. )");
    }

    int nbr_of_pair = combination(2, data_plane_vec.nbr_of_indiv());

    //numerator and denominator
    std::vector<std::array<double, 2>> result(nbr_of_pair, {0, 0});
    //In prob_id case Qw with a missing value = 0; Denom by locus will be the same for all pair
    std::vector<double> Qw_by_locus(data_plane_vec.nbr_locus(), 0);
    //
    double sum_Qw_all_loc = 0;

    //Multilocus estimates are defined as the sum of locus-specific numerators divided by the sum of locus-specific denominators (Weir & Cockerham, 1984)
    for (int locus = 0; locus < data_plane_vec.nbr_locus(); ++locus)
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
                            result_itr->at(0) = data_plane_vec.geo_dist_btw_gene(indiv * Ploidy, next_indiv * Ploidy);
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
        Qw_by_locus[locus] /= data_plane_vec.nomiss_nbr_of_indiv(locus);
        sum_Qw_all_loc += Qw_by_locus[locus];
    }

    //compute denom because in missing value denom will be diff between pair (When the data_plane_vec is unknown for some locus in one individual from a pair, this locus is not included in the sum over loci Rousset, 2000)
    //compute by remove Qw from sum_Qw_all_loc
    std::vector<std::array<double, 2>> sum_Qw_use_by_pair(nbr_of_pair, {sum_Qw_all_loc, static_cast<double>(data_plane_vec.nbr_locus())});

    auto sum_Qw_use_by_pair_itr = sum_Qw_use_by_pair.begin();
    for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv() - 1; ++indiv)
    {
        for (int next_indiv = indiv + 1; next_indiv < data_plane_vec.nbr_of_indiv(); ++next_indiv)
        {
            std::vector<bool> nomiss_data = bin_vec::and_(data_plane_vec.nomiss_data(indiv), data_plane_vec.nomiss_data(next_indiv));
            for (int locus = 0; locus < data_plane_vec.nbr_locus(); ++locus)
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

std::vector<std::array<double, 2>> er_by_indiv_pair(data_plane_vec_c const &data_plane_vec)
{
    int Ploidy = data_plane_vec.get_Ploidy();
    if (Ploidy != 2)
    {
        throw std::logic_error("( Can't calculate ar if ploidy was different than 2. I exit. )");
    }

    int nbr_indiv = data_plane_vec.nbr_of_indiv();
    int nbr_of_pair = combination(2, nbr_indiv);

    //numerator and denominator
    std::vector<std::array<double, 2>> save_value(nbr_of_pair, {0, 0});
    //To calculate e_r like Genedeme (Rousset, ...) with Loiselle F statistic equivalent d_k^2 terme
    std::vector<double> sum_Qij_by_locus(data_plane_vec.nbr_locus(), 0);
    //In prob_id case Qw with a missing value = 0; Denom by locus will be the same for all pair
    std::vector<double> Qw_by_locus(data_plane_vec.nbr_locus(), 0);
    double sum_Qw_all_loc = 0;

    //All Qi need to be stored for later correction
    std::vector<std::vector<double>> Qi_by_indiv_by_locus(nbr_indiv, std::vector<double>(data_plane_vec.nbr_locus(), 0));
    std::vector<double> Qi_sum_all_loc(nbr_indiv, 0);

    //Multilocus estimates are defined as the sum of locus-specific numerators divided by the sum of locus-specific denominators (Weir & Cockerham, 1984)
    for (int locus = 0; locus < data_plane_vec.nbr_locus(); ++locus)
    {
        auto nomiss_indiv = data_plane_vec.nomiss_nbr_of_indiv(locus);
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
                            save_value_itr->at(0) = data_plane_vec.geo_dist_btw_gene(indiv * Ploidy, next_indiv * Ploidy);
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
        Qw_by_locus[locus] /= data_plane_vec.nomiss_nbr_of_indiv(locus);
        sum_Qw_all_loc += Qw_by_locus[locus];
    }

    //calc of terme for match Genedeme based on Loisel F statistic

    std::vector<double> Qi_use_by_pair(nbr_of_pair, 0);
    //nbr of pair use to calc sum_Qij (bijection with the number of indiv)
    std::vector<double> l_term_use_by_pair(nbr_of_pair, 0);
    //1 - sum_Qw_all_loc to match denom in cal er
    std::vector<double> sum_Qw_use_by_pair(nbr_of_pair, data_plane_vec.nbr_locus() - sum_Qw_all_loc);

    auto Qi_use_by_pair_itr = Qi_use_by_pair.begin();
    auto l_term_use_by_pair_itr = l_term_use_by_pair.begin();
    auto sum_Qw_use_by_pair_itr = sum_Qw_use_by_pair.begin();
    for (int indiv = 0; indiv < nbr_indiv - 1; ++indiv)
    {
        for (int next_indiv = indiv + 1; next_indiv < nbr_indiv; ++next_indiv)
        {
            //file Qi_use_by_pair with sum of Qi and Qj (if i or j have missing value, Qi or Qj = 0, just need to remove the other one)
            *Qi_use_by_pair_itr = Qi_sum_all_loc[indiv] + Qi_sum_all_loc[next_indiv];
            std::vector<bool> nomiss_data = bin_vec::and_(data_plane_vec.nomiss_data(indiv), data_plane_vec.nomiss_data(next_indiv));
            for (int locus = 0; locus < data_plane_vec.nbr_locus(); ++locus)
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
                    int nomiss_nbr_indiv = data_plane_vec.nomiss_nbr_of_indiv(locus);
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
//WARNING : If no missing value and indiv_per_deme > 1
std::array<std::array<double, 2>, 2> Fstat_per_chr_by_loc_with_probid(data_plane_vec_c const &data_plane_vec, int chr, int locus)
{
    std::array<std::array<double, 2>, 2> result{{{0}, {0}}};

    int Ploidy = data_plane_vec.get_Ploidy();
    if (Ploidy != 2)
    {
        throw std::logic_error("( In Fstat_per_chr_by_loc_with_probid : Can't calculate Fstat if ploidy was different than 2. I exit. )");
    }
    //
    auto Qwi_frac = std::array<int, 2>{0, 0};
    //+
    auto Qwp_frac_v = std::vector<std::array<int, 2>>(data_plane_vec.nbr_of_deme(), {0, 0});
    auto Qbp_frac = std::array<int, 2>{0, 0};

    int locus_end = data_plane_vec.index_end_locus(chr, locus);
    auto Qwp_frac_v_itr = Qwp_frac_v.begin();
    auto cumul_deme_size_itr = data_plane_vec.cumul_nbr_of_indiv_per_deme().begin() + 1;
    int indiv_relatif_nbr = -1;

    for (int gene1 = data_plane_vec.index_begin_locus(chr, locus); gene1 < locus_end - 1; ++gene1)
    {
        ++indiv_relatif_nbr;
        if (indiv_relatif_nbr == (*cumul_deme_size_itr) * data_plane_vec.get_Ploidy())
        {
            ++cumul_deme_size_itr;
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
                        if (data_plane_vec.same_deme(gene1, gene2))
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
    double ns = data_plane_vec.nomiss_nbr_of_deme(chr, locus);

    for (int deme = 0; deme < data_plane_vec.nbr_of_deme(); ++deme)
    {
        int deme_size = data_plane_vec.nomiss_nbr_of_indiv(chr, locus, deme);
        S1 += deme_size;
        S2 += pow(deme_size, 2);
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
//WARNING : Work if indiv_per_deme > 1 but slower
std::array<std::array<double, 2>, 2> Fstat_per_chr_by_loc_with_indic(data_plane_vec_c const &data_plane_vec, int chr, int locus)
{
    auto const &allele_state = data_plane_vec.allele_state(chr, locus);

    int Ploidy = data_plane_vec.get_Ploidy();
    if (Ploidy != 2)
    {
        throw std::logic_error("( In Fstat_per_chr_by_loc_with_indic : Can't calculate Fstat if ploidy was different than 2. I exit. )");
    }

    double S1 = 0;
    double S2 = 0;
    double ns = data_plane_vec.nomiss_nbr_of_deme(chr, locus);

    for (int deme = 0; deme < data_plane_vec.nbr_of_deme(); ++deme)
    {
        int deme_size = data_plane_vec.nomiss_nbr_of_indiv(chr, locus, deme);
        S1 += deme_size;
        S2 += pow(deme_size, 2);
    }
    //Indicatrix matrix mean for each deme
    std::vector<std::vector<double>> deme_state_mean(data_plane_vec.nbr_of_deme(), std::vector<double>(allele_state.size(), 0));

    //Calc mean for sample
    std::vector<double> state_mean_samp(allele_state.size(), 0);
    for (int deme = 0; deme < data_plane_vec.nbr_of_deme(); ++deme)
    {
        //for handle missing value need to remove all indiv with missing value
        int gene_in_nomiss_nbr_of_indiv_per_deme = data_plane_vec.nomiss_nbr_of_indiv(chr, locus, deme) * 2;
        int tot_gene_in_nomiss_nbr_of_indiv = (S1 * 2);
        for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv(deme); ++indiv)
        {
            int indiv1_gene1 = data_plane_vec(chr, locus, deme, indiv, 0);
            int indiv1_gene2 = data_plane_vec(chr, locus, deme, indiv, 1);

            if ((indiv1_gene1 != 0) && (indiv1_gene2 != 0))
            {
                int count = 0;
                for (auto const &state : allele_state)
                {
                    double temp = (static_cast<double>(indiv1_gene1 == state.first) + static_cast<double>(indiv1_gene2 == state.first));
                    deme_state_mean[deme][count] += (temp / gene_in_nomiss_nbr_of_indiv_per_deme);
                    state_mean_samp[count] += (temp / tot_gene_in_nomiss_nbr_of_indiv);
                    ++count;
                }
            }
        }
    }

    double SSg = 0, SSi = 0, SSp = 0;

    for (int deme = 0; deme < data_plane_vec.nbr_of_deme(); ++deme)
    {
        for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv(deme); ++indiv)
        {
            int indiv1_gene1 = data_plane_vec(chr, locus, deme, indiv, 0);
            int indiv1_gene2 = data_plane_vec(chr, locus, deme, indiv, 1);

            if ((indiv1_gene1 != 0) && (indiv1_gene2 != 0))
            {
                int count = 0;
                for (auto const &state : allele_state)
                {
                    int Xij1 = (indiv1_gene1 == state.first);
                    int Xij2 = (indiv1_gene2 == state.first);
                    double Xij = (Xij1 + Xij2) / 2.0;
                    SSg += pow(Xij1 - Xij, 2) + pow(Xij2 - Xij, 2);
                    SSi += (pow(Xij - deme_state_mean[deme][count], 2)) * 2;
                    //Two genes
                    SSp += (pow(deme_state_mean[deme][count] - state_mean_samp[count], 2) * 2);
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

std::array<double, 2> Fstat_genepop(data_plane_vec_c const &data_plane_vec, bool mising_data)
{
    std::array<double, 2> result{0, 0};
    double fis_num = 0, fis_denum = 0, fst_num = 0, fst_denum = 0;
    for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
    {
        for (int locus = 0; locus < data_plane_vec.nbr_locus(chr); ++locus)
        {
            std::array<std::array<double, 2>, 2> temp;
            if (mising_data)
            {
                temp = Fstat_per_chr_by_loc_with_indic(data_plane_vec, chr, locus);
            }
            else
            {
                temp = Fstat_per_chr_by_loc_with_probid(data_plane_vec, chr, locus);
            }
            fis_num += temp.at(0).at(0);
            fis_denum += temp.at(0).at(1);
            fst_num += temp.at(1).at(0);
            fst_denum += temp.at(1).at(1);
        }
    }

    result.at(0) = fis_num / fis_denum;
    result.at(1) = fst_num / fst_denum;
    return result;
}

//Return <state, <frequence, nbr of locus>>
std::map<int, std::map<int, double>> calc_AFS(data_plane_vec_c const &data_plane_vec, int limit_min_gene_per_locus)
{
    //chr<>
    std::vector<std::vector<int>> locus_usable_for_AFS;
    locus_usable_for_AFS.resize(data_plane_vec.nbr_of_chr());

    int real_min_limit = data_plane_vec.nbr_of_gene_per_loc();
    bool can_calc_AFS = false;

    for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
    {
        locus_usable_for_AFS[chr].reserve(data_plane_vec.nbr_locus(chr));
        for (int locus = 0; locus < data_plane_vec.nbr_locus(chr); ++locus)
        {
            int nbr_gene = data_plane_vec.nomiss_nbr_of_gene(chr, locus);
            if (nbr_gene >= limit_min_gene_per_locus)
            {
                locus_usable_for_AFS[chr].push_back(locus);
                if (nbr_gene < real_min_limit)
                {
                    real_min_limit = nbr_gene;
                }
            }
        }
        locus_usable_for_AFS[chr].shrink_to_fit();
        if (locus_usable_for_AFS[chr].size() > 0)
        {
            can_calc_AFS = true;
        }
    }

    if (!can_calc_AFS)
    {
        throw std::logic_error("Can't calculate AFS, no locus with number of non missing gene > 'min gene'. I exit.");
    }

    std::map<int, std::map<int, double>> result;

    for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
    {
        for (int locus = 0; locus < data_plane_vec.nbr_locus(chr); ++locus)
        {
            for (auto state : data_plane_vec.allele_state(chr, locus))
            {
                result.emplace(state.first, std::map<int, double>{});
            }
        }
    }

    if (real_min_limit == data_plane_vec.nbr_of_gene_per_loc())
    {
        for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
        {
            for (int locus = 0; locus < data_plane_vec.nbr_locus(chr); ++locus)
            {
                //map(state, nbr of allele in this state)
                auto temp_locus = data_plane_vec.allele_state(chr, locus);
                for (auto state : data_plane_vec.allele_state(chr, locus))
                {
                    auto pair = result.at(state.first).emplace(state.second, 1.0);
                    if (!pair.second)
                    {
                        pair.first->second += 1.0;
                    }
                }
            }
        }
    }
    else
    {
        for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
        {
            for (auto locus : locus_usable_for_AFS.at(chr))
            {
                //map(state, nbr of allele in this state)
                int n = data_plane_vec.nomiss_nbr_of_gene(chr, locus);
                for (auto state : data_plane_vec.allele_state(chr, locus))
                {
                    // k = real_min_limit; n = nomiss_nbr_of_gene_per_loc ; nb = state.second ; b = state
                    int nb = state.second;
                    double total_comb = combination(real_min_limit, n);
                    int max_count_stat = min(real_min_limit, nb);
                    int min_count_stat = max(0, real_min_limit + nb - n); // k-na : na = n-nb

                    for (int count = min_count_stat; count <= max_count_stat; ++count)
                    {
                        double freq = combination(count, nb) * combination(real_min_limit - count, n - nb) / total_comb;
                        auto pair = result.at(state.first).emplace(count, freq);
                        if (!pair.second)
                        {
                            pair.first->second += freq;
                        }
                    }
                }
            }
        }
    }

    return result;
}