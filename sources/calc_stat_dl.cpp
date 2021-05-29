#include "calc_stat_dl.hpp"
#include "calc_stat.hpp"

//<locus<chr>>
double calc_phi_ij_xy(std::array<std::array<int, 4>, 2> const &locus_chr)
{
    double minor_phi = (((locus_chr[0][0] == locus_chr[0][2]) && (locus_chr[1][0] == locus_chr[1][2])) + ((locus_chr[0][0] == locus_chr[0][3]) && (locus_chr[1][0] == locus_chr[1][3])) + ((locus_chr[0][1] == locus_chr[0][2]) && (locus_chr[1][1] == locus_chr[1][2])) + ((locus_chr[0][1] == locus_chr[0][3]) && (locus_chr[1][1] == locus_chr[1][3]))) / 16.0;
    double gamma1 = (((locus_chr[0][0] == locus_chr[0][2]) && (locus_chr[1][0] == locus_chr[1][3])) + ((locus_chr[0][0] == locus_chr[0][3]) && (locus_chr[1][0] == locus_chr[1][2])) + ((locus_chr[0][1] == locus_chr[0][2]) && (locus_chr[1][1] == locus_chr[1][3])) + ((locus_chr[0][1] == locus_chr[0][3]) && (locus_chr[1][1] == locus_chr[1][2]))) / 16.0;
    double gamma2 = (((locus_chr[0][0] == locus_chr[0][2]) && (locus_chr[1][1] == locus_chr[1][2])) + ((locus_chr[0][1] == locus_chr[0][2]) && (locus_chr[1][0] == locus_chr[1][2])) + ((locus_chr[0][0] == locus_chr[0][3]) && (locus_chr[1][1] == locus_chr[1][3])) + ((locus_chr[0][1] == locus_chr[0][3]) && (locus_chr[1][0] == locus_chr[1][3]))) / 16.0;
    double delta = (((locus_chr[0][0] == locus_chr[0][2]) && (locus_chr[1][1] == locus_chr[1][3])) + ((locus_chr[0][0] == locus_chr[0][3]) && (locus_chr[1][1] == locus_chr[1][2])) + ((locus_chr[0][1] == locus_chr[0][2]) && (locus_chr[1][0] == locus_chr[1][3])) + ((locus_chr[0][1] == locus_chr[0][3]) && (locus_chr[1][0] == locus_chr[1][2]))) / 16.0;
    return (minor_phi + gamma1 + gamma2 + delta);
}
//Phi for all indiv, calc without missing data
std::vector<double> calc_phi_ij(data_plane_vec_c const &data_plane_vec, int ploidy)
{
    int locus_pair_nbr = (data_plane_vec.nbr_locus() * (data_plane_vec.nbr_locus() - 1)) / 2;

    std::vector<double> result;
    result.reserve(locus_pair_nbr);

    for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
    {
        for (int locus_i = 0; locus_i < data_plane_vec.nbr_locus(); ++locus_i)
        {
            for (int locus_j = locus_i + 1; locus_j < data_plane_vec.nbr_locus(); ++locus_j)
            {
                int indiv_pair_nbr = 0;
                double value = 0;
                if (ploidy == 2)
                {
                    for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv(); ++indiv)
                    {
                        int locus_i_indiv1_gene1 = data_plane_vec(chr, locus_i, indiv, 0);
                        int locus_i_indiv1_gene2 = data_plane_vec(chr, locus_i, indiv, 1);

                        int locus_j_indiv1_gene1 = data_plane_vec(chr, locus_j, indiv, 0);
                        int locus_j_indiv1_gene2 = data_plane_vec(chr, locus_j, indiv, 1);

                        for (int other_indiv = indiv + 1; other_indiv < data_plane_vec.nbr_of_indiv(); ++other_indiv)
                        {
                            int locus_i_indiv2_gene1 = data_plane_vec(chr, locus_i, other_indiv, 0);
                            int locus_i_indiv2_gene2 = data_plane_vec(chr, locus_i, other_indiv, 1);

                            int locus_j_indiv2_gene1 = data_plane_vec(chr, locus_j, other_indiv, 0);
                            int locus_j_indiv2_gene2 = data_plane_vec(chr, locus_j, other_indiv, 1);

                            value += calc_phi_ij_xy({{{locus_i_indiv1_gene1, locus_i_indiv1_gene2, locus_i_indiv2_gene1, locus_i_indiv2_gene2},
                                                      {locus_j_indiv1_gene1, locus_j_indiv1_gene2, locus_j_indiv2_gene1, locus_j_indiv2_gene2}}});
                            ++indiv_pair_nbr;
                        }
                    }
                }
                else
                {
                    for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv(); ++indiv)
                    {
                        int locus_i_indiv1 = data_plane_vec(chr, locus_i, indiv, 0);
                        int locus_j_indiv1 = data_plane_vec(chr, locus_j, indiv, 0);

                        for (int other_indiv = indiv + 1; other_indiv < data_plane_vec.nbr_of_indiv(); ++other_indiv)
                        {
                            int locus_i_indiv2 = data_plane_vec(chr, locus_i, other_indiv, 0);
                            int locus_j_indiv2 = data_plane_vec(chr, locus_j, other_indiv, 0);

                            value += (locus_i_indiv1 == locus_i_indiv2) && (locus_j_indiv1 == locus_j_indiv2);
                            ++indiv_pair_nbr;
                        }
                    }
                }
                result.push_back(value / indiv_pair_nbr);
            }
        }
    }
    return result;
}

//<locus_pair_nbr<phi, prob_intra_loc, denom>>
std::array<double, 3> calc_eta_ij_xy(data_plane_vec_c const &data_plane_vec, int chr, int locus_i, double Q2_loc_i, double Q1_loc_i, int locus_j, double Q2_loc_j, double Q1_loc_j, int deme_x, int deme_y)
{
    double phi = 0;
    double div = 0;
    if (data_plane_vec.get_Ploidy() == 2)
    {
        for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv_per_deme(deme_x); ++indiv)
        {
            int locus_i_indiv1_gene1 = data_plane_vec(chr, locus_i, deme_x, indiv, 0);
            int locus_i_indiv1_gene2 = data_plane_vec(chr, locus_i, deme_x, indiv, 1);

            int locus_j_indiv1_gene1 = data_plane_vec(chr, locus_j, deme_x, indiv, 0);
            int locus_j_indiv1_gene2 = data_plane_vec(chr, locus_j, deme_x, indiv, 1);
            //In AAaaBBbb if A&B missing => phi = 0, etc
            if (((locus_i_indiv1_gene1 != 0) || (locus_i_indiv1_gene2 != 0)) && ((locus_j_indiv1_gene1 != 0) || (locus_j_indiv1_gene2 != 0)))
            {
                for (int indiv_other_deme = 0; indiv_other_deme < data_plane_vec.nbr_of_indiv_per_deme(deme_y); ++indiv_other_deme)
                {
                    int locus_i_indiv2_gene1 = data_plane_vec(chr, locus_i, deme_y, indiv_other_deme, 0);
                    int locus_i_indiv2_gene2 = data_plane_vec(chr, locus_i, deme_y, indiv_other_deme, 1);

                    int locus_j_indiv2_gene1 = data_plane_vec(chr, locus_j, deme_y, indiv_other_deme, 0);
                    int locus_j_indiv2_gene2 = data_plane_vec(chr, locus_j, deme_y, indiv_other_deme, 1);
                    if (((locus_i_indiv2_gene1 != 0) || (locus_i_indiv2_gene2 != 0)) && ((locus_j_indiv2_gene1 != 0) || (locus_j_indiv2_gene2 != 0)))
                    {
                        int miss = 0;
                        double pond;
                        //handle missing value
                        std::array<std::array<int, 4>, 2> locus_chr = {
                            {{locus_i_indiv1_gene1, locus_i_indiv1_gene2, locus_i_indiv2_gene1, locus_i_indiv2_gene2},
                             {locus_j_indiv1_gene1, locus_j_indiv1_gene2, locus_j_indiv2_gene1, locus_j_indiv2_gene2}}};
                        for (auto &value_loc : locus_chr)
                        {
                            for (auto &value : value_loc)
                            { //To avoid 0==0 test in phi calculation
                                if (value == 0)
                                {
                                    value = --miss;
                                }
                            }
                        }
                        //pond depend of missing value number
                        switch (miss)
                        {
                        case 0:
                        {
                            pond = 1;
                            break;
                        }
                        case -1:
                        {
                            pond = 0.5;
                            break;
                        }
                        case -2:
                        {
                            pond = 0.25;
                            break;
                        }
                        case -3:
                        {
                            pond = 0.125;
                            break;
                        }
                        default:
                        {
                            pond = 0.0625;
                            break;
                        }
                        }

                        phi += pond * calc_phi_ij_xy(locus_chr);
                        div += pond;
                    }
                }
            }
        }
    }
    else
    {
        for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv_per_deme(deme_x); ++indiv)
        {
            int locus_i_indiv1 = data_plane_vec(chr, locus_i, deme_x, indiv, 0);
            int locus_j_indiv1 = data_plane_vec(chr, locus_j, deme_x, indiv, 0);
            if ((locus_i_indiv1 != 0) && (locus_j_indiv1 != 0))
            {
                for (int indiv_other_deme = 0; indiv_other_deme < data_plane_vec.nbr_of_indiv_per_deme(deme_y); ++indiv_other_deme)
                {
                    int locus_i_indiv2 = data_plane_vec(chr, locus_i, deme_y, indiv_other_deme, 0);
                    int locus_j_indiv2 = data_plane_vec(chr, locus_j, deme_y, indiv_other_deme, 0);
                    if ((locus_i_indiv2 != 0) && (locus_j_indiv2 != 0))
                    {
                        phi += (locus_i_indiv1 == locus_i_indiv2) && (locus_j_indiv1 == locus_j_indiv2);
                        ++div;
                    }
                }
            }
        }
    }

    return std::array<double, 3>{phi / div, Q1_loc_i * Q1_loc_j, (1 - Q2_loc_i) * (1 - Q2_loc_j)};
}

//<deme_pair_nbr,<dist-deme, dist-locus, value eta, value eta denom>>
std::vector<std::array<double, 5>> calc_eta_ij(data_plane_vec_c const &data_plane_vec, int chr, int locus_i, double Q2_loc_i, int locus_j, double Q2_loc_j)
{
    int deme_pair_nbr = (data_plane_vec.nbr_of_deme() * (data_plane_vec.nbr_of_deme() - 1)) / 2;
    std::vector<std::array<double, 5>> result;
    result.reserve(deme_pair_nbr);

    std::vector<std::array<int, 2>> Q1_locus_i_deme;
    Q1_locus_i_deme.reserve(data_plane_vec.nbr_of_deme());

    std::vector<std::array<int, 2>> Q1_locus_j_deme;
    Q1_locus_j_deme.reserve(data_plane_vec.nbr_of_deme());

    for (int deme = 0; deme < data_plane_vec.nbr_of_deme(); ++deme)
    {
        Q1_locus_i_deme.push_back(calc_Q_inter_indiv_per_chr_per_locus_per_deme(data_plane_vec, chr, locus_i, deme));

        Q1_locus_j_deme.push_back(calc_Q_inter_indiv_per_chr_per_locus_per_deme(data_plane_vec, chr, locus_j, deme));
    }

    double dist_locus;
    if (data_plane_vec.nbr_chr_dist_class() == 0)
    {
        dist_locus = data_plane_vec.chr_dist_btw_locus(chr, locus_i, locus_j);
    }
    else
    {
        dist_locus = data_plane_vec.chr_dist_class_btw_locus(chr, locus_i, locus_j);
    }

    for (int deme_x = 0; deme_x < data_plane_vec.nbr_of_deme(); ++deme_x)
    {
        for (int deme_y = deme_x + 1; deme_y < data_plane_vec.nbr_of_deme(); ++deme_y)
        {
            auto Q1_loc_i_xy = static_cast<double>(Q1_locus_i_deme[deme_x].at(0) + Q1_locus_i_deme[deme_y].at(0)) / (Q1_locus_i_deme[deme_x].at(1) + Q1_locus_i_deme[deme_y].at(1));
            auto Q1_loc_j_xy = static_cast<double>(Q1_locus_j_deme[deme_x].at(0) + Q1_locus_j_deme[deme_y].at(0)) / (Q1_locus_j_deme[deme_x].at(1) + Q1_locus_j_deme[deme_y].at(1));
            auto eta = calc_eta_ij_xy(data_plane_vec, chr, locus_i, Q2_loc_i, Q1_loc_i_xy, locus_j, Q2_loc_j, Q1_loc_j_xy, deme_x, deme_y);
            //dist-deme, dist-locus, value eta
            double dist_deme;
            if (data_plane_vec.nbr_geo_dist_class() > 1)
            {
                dist_deme = data_plane_vec.geo_dist_class_btw_deme(deme_x, deme_y);
            }
            else
            {
                dist_deme = data_plane_vec.geo_dist_btw_deme(deme_x, deme_y);
            }
            result.push_back({{dist_deme, dist_locus, eta.at(0), eta.at(1), eta.at(2)}});
        }
    }

    return result;
}

//Haploide version with > 1 indiv/deme
//<pair of deme * pair of locus,<dist-locus, dist-deme, value eta>>
std::vector<std::array<double, 5>> calc_eta(data_plane_vec_c const &data_plane_vec)
{
    std::vector<std::array<double, 5>> result;
    for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
    {
        std::vector<int> const &poly_loc = data_plane_vec.polymorph_locus_list(chr);
        //TODO : Need to be handle up to avoid general throws
        if (poly_loc.size() < 2)
        {
            throw std::logic_error("Less than 2 polymorphics locus in data file. Can't calculate eta.");
        }

        int locus_pair_nbr = (poly_loc.size() * (poly_loc.size() - 1)) / 2;
        int deme_pair_nbr = (data_plane_vec.nbr_of_deme() * (data_plane_vec.nbr_of_deme() - 1)) / 2;

        result.reserve(result.size() + locus_pair_nbr * deme_pair_nbr);

        std::vector<double> Q2_locus;
        Q2_locus.reserve(poly_loc.size());
        for (int locus : poly_loc)
        {
            auto temp = calc_Q_inter_deme_per_chr_per_locus(data_plane_vec, chr, locus);
            Q2_locus.push_back(static_cast<double>(temp.at(0)) / temp.at(1));
        }

        for (auto locus_i = poly_loc.begin(); locus_i < poly_loc.end(); ++locus_i)
        {
            for (auto locus_j = locus_i + 1; locus_j < poly_loc.end(); ++locus_j)
            {
                auto temp_vec = calc_eta_ij(data_plane_vec, chr, *locus_i, Q2_locus[*locus_i], *locus_j, Q2_locus[*locus_j]);
                result.insert(result.end(), temp_vec.begin(), temp_vec.end());
            }
        }
    }

    return result;
}

//Diploide version with > 1 indiv/deme
std::vector<std::array<double, 5>> calc_eta_q1_version(data_plane_vec_c const &data_plane_vec)
{
    std::vector<std::array<double, 5>> result;
    for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
    {
        std::vector<int> const &poly_loc = data_plane_vec.polymorph_locus_list(chr);
        //TODO : Need to be handle up to avoid general throws
        if (poly_loc.size() < 2)
        {
            throw std::logic_error("Less than 2 polymorphics locus in data file. Can't calculate eta.");
        }

        int locus_pair_nbr = (poly_loc.size() * (poly_loc.size() - 1)) / 2;
        int deme_pair_nbr = (data_plane_vec.nbr_of_deme() * (data_plane_vec.nbr_of_deme() - 1)) / 2;

        result.reserve(result.size() + locus_pair_nbr * deme_pair_nbr);

        for (auto locus_i = poly_loc.begin(); locus_i < poly_loc.end(); ++locus_i)
        {
            for (auto locus_j = locus_i + 1; locus_j < poly_loc.end(); ++locus_j)
            {
                std::vector<std::array<int, 2>> Q1_locus_i_deme;
                Q1_locus_i_deme.reserve(data_plane_vec.nbr_of_deme());

                std::vector<std::array<int, 2>> Q1_locus_j_deme;
                Q1_locus_j_deme.reserve(data_plane_vec.nbr_of_deme());

                for (int deme = 0; deme < data_plane_vec.nbr_of_deme(); ++deme)
                {
                    Q1_locus_i_deme.push_back(calc_Q_inter_indiv_per_chr_per_locus_per_deme(data_plane_vec, chr, *locus_i, deme));

                    Q1_locus_j_deme.push_back(calc_Q_inter_indiv_per_chr_per_locus_per_deme(data_plane_vec, chr, *locus_j, deme));
                }

                double dist_locus;

                if (data_plane_vec.nbr_chr_dist_class() == 0)
                {
                    dist_locus = data_plane_vec.chr_dist_btw_locus(chr, *locus_i, *locus_j);
                }
                else
                {
                    dist_locus = data_plane_vec.chr_dist_class_btw_locus(chr, *locus_i, *locus_j);
                }

                for (int deme_x = 0; deme_x < data_plane_vec.nbr_of_deme(); ++deme_x)
                {
                    for (int deme_y = deme_x + 1; deme_y < data_plane_vec.nbr_of_deme(); ++deme_y)
                    {
                        auto Q1_loc_i_xy = static_cast<double>(Q1_locus_i_deme[deme_x].at(0) + Q1_locus_i_deme[deme_y].at(0)) / (Q1_locus_i_deme[deme_x].at(1) + Q1_locus_i_deme[deme_y].at(1));
                        auto Q1_loc_j_xy = static_cast<double>(Q1_locus_j_deme[deme_x].at(0) + Q1_locus_j_deme[deme_y].at(0)) / (Q1_locus_j_deme[deme_x].at(1) + Q1_locus_j_deme[deme_y].at(1));
                        auto eta = calc_eta_ij_xy(data_plane_vec, chr, *locus_i, Q1_loc_i_xy, Q1_loc_i_xy, *locus_j, Q1_loc_j_xy, Q1_loc_j_xy, deme_x, deme_y);
                        //dist-deme, dist-locus, value eta
                        double dist_deme;
                        if (data_plane_vec.nbr_geo_dist_class() > 1)
                        {
                            dist_deme = data_plane_vec.geo_dist_class_btw_deme(deme_x, deme_y);
                        }
                        else
                        {
                            dist_deme = data_plane_vec.geo_dist_btw_deme(deme_x, deme_y);
                        }
                        result.push_back({{dist_deme, dist_locus, eta.at(0), eta.at(1), eta.at(2)}});
                    }
                }
            }
        }
    }
    return result;
}
//Continous habitat isolation by distance
std::vector<std::array<double, 5>> calc_eta_1_indiv_deme_v(data_plane_vec_c const &data_plane_vec)
{
    std::vector<std::array<double, 5>> result;
    for (int chr = 0; chr < data_plane_vec.nbr_of_chr(); ++chr)
    {
        std::vector<int> const &poly_loc = data_plane_vec.polymorph_locus_list(chr);
        //TODO : Need to be handle up to avoid general throws
        if (poly_loc.size() < 2)
        {
            throw std::logic_error("Less than 2 polymorphics locus in data file. Can't calculate eta.");
        }

        int locus_pair_nbr = (poly_loc.size() * (poly_loc.size() - 1)) / 2;
        int deme_pair_nbr = (data_plane_vec.nbr_of_deme() * (data_plane_vec.nbr_of_deme() - 1)) / 2;

        result.reserve(result.size() + locus_pair_nbr * deme_pair_nbr);

        std::vector<double> Q2_locus;
        Q2_locus.reserve(poly_loc.size());
        std::vector<double> Q0_locus;
        Q0_locus.reserve(poly_loc.size());
        for (int locus : poly_loc)
        {
            auto temp = calc_Q_inter_deme_per_chr_per_locus(data_plane_vec, chr, locus);
            Q2_locus.push_back(static_cast<double>(temp.at(0)) / temp.at(1));
            if (data_plane_vec.get_Ploidy() == 2)
            {
                temp = calc_Q_intra_indiv_per_chr_per_locus(data_plane_vec, chr, locus);
                Q0_locus.push_back(static_cast<double>(temp.at(0)) / temp.at(1));
            }
        }

        for (auto locus_i = poly_loc.begin(); locus_i < poly_loc.end(); ++locus_i)
        {
            for (auto locus_j = locus_i + 1; locus_j < poly_loc.end(); ++locus_j)
            {

                double dist_locus;

                if (data_plane_vec.nbr_chr_dist_class() == 0)
                {
                    dist_locus = data_plane_vec.chr_dist_btw_locus(chr, *locus_i, *locus_j);
                }
                else
                {
                    dist_locus = data_plane_vec.chr_dist_class_btw_locus(chr, *locus_i, *locus_j);
                }

                for (int deme_x = 0; deme_x < data_plane_vec.nbr_of_deme(); ++deme_x)
                {
                    for (int deme_y = deme_x + 1; deme_y < data_plane_vec.nbr_of_deme(); ++deme_y)
                    {
                        std::array<double, 3> eta;
                        if (data_plane_vec.get_Ploidy() == 2)
                        {
                            eta = calc_eta_ij_xy(data_plane_vec, chr, *locus_i, Q2_locus[*locus_i], Q0_locus[*locus_i], *locus_j, Q2_locus[*locus_j], Q0_locus[*locus_j], deme_x, deme_y);
                        }
                        if (data_plane_vec.get_Ploidy() == 1)
                        {
                            eta = calc_eta_ij_xy(data_plane_vec, chr, *locus_i, Q2_locus[*locus_i], Q2_locus[*locus_i], *locus_j, Q2_locus[*locus_j], Q2_locus[*locus_j], deme_x, deme_y);
                        }
                        //dist-deme, dist-locus, value eta
                        double dist_deme;
                        if (data_plane_vec.nbr_geo_dist_class() > 1)
                        {
                            dist_deme = data_plane_vec.geo_dist_class_btw_deme(deme_x, deme_y);
                        }
                        else
                        {
                            dist_deme = data_plane_vec.geo_dist_btw_deme(deme_x, deme_y);
                        }
                        result.push_back({{dist_deme, dist_locus, eta.at(0), eta.at(1), eta.at(2)}});
                    }
                }
            }
        }
    }
    return result;
}

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "stdafx.h"
#include "interpolation.h"
#include "common_tools.hpp"
using namespace alglib;
std::array<double, 3> exp_regr(std::vector<std::array<double, 3>> const &dist_geo_eta_weights)
{
    //a+ b * (1-exp(-b_g * Geo_dist_btw_deme))
    //c[0] a ; c[1] = b ; c[2] = b_g ; x[0] = Geo_dist_btw_deme
    auto function_cx_1_func = [](const real_1d_array &c, const real_1d_array &x, double &func, void *ptr)
    {
        // this callback calculates f(c,x)=c[0]+ c[1] * (1-exp(-c[2] *x)
        // where x = Dist_btw_locus_pb and c are adjustables parameters
        func = c[0] + c[1] * (1 - exp(-c[2] * x[0]));
    };
    auto function_cx_1_grad = [](const real_1d_array &c, const real_1d_array &x, double &func, real_1d_array &grad, void *ptr)
    {
        // this callback calculates f(c,x)=c[0]+ c[1] * (1-exp(-c[2] *x) and gradient G={df/dc[i]}
        // where x = Dist_btw_locus_pb and c are adjustables parameters
        // IMPORTANT: gradient is calculated with respect to C, not to X
        func = c[0] + c[1] * (1 - exp(-c[2] * x[0]));
        grad[0] = 1;
        grad[1] = (1 - exp(-c[2] * x[0]));
        grad[2] = (c[1] * x[0] * exp(-c[2] * x[0]));
    };
    auto function_cx_1_hess = [](const real_1d_array &c, const real_1d_array &x, double &func, real_1d_array &grad, real_2d_array &hess, void *ptr)
    {
        // this callback calculates f(c,x)=c[0]+ c[1] * (1-exp(-c[2] *x), gradient G={df/dc[i]} and Hessian H={d2f/(dc[i]*dc[j])}
        // where x = Dist_btw_locus_pb and c are adjustables parameters
        // IMPORTANT: gradient/Hessian are calculated with respect to C, not to X
        func = c[0] + c[1] * (1 - exp(-c[2] * x[0]));

        grad[0] = 1;
        grad[1] = 1 - exp(-c[2] * x[0]);
        grad[2] = c[1] * x[0] * exp(-c[2] * x[0]);

        // Hessian
        //                     a               b               b_g
        // a    [              0               0               0]
        // b    [              0               0      x*e^(-c*x)]
        // b_g  [              0      x*e^(-c*x) -b*x^2*e^(-c*x)]

        hess[0][0] = 0;
        hess[0][1] = 0;
        hess[0][2] = 0;

        hess[1][0] = 0;
        hess[1][1] = 0;
        hess[1][2] = x[0] * exp(-c[2] * x[0]);

        hess[2][0] = 0;
        hess[2][1] = x[0] * exp(-c[2] * x[0]);
        hess[2][2] = -c[1] * pow(x[0], 2) * exp(-c[2] * x[0]);

        // std::cout << "Gradient" << std::endl;
        // std::cout << "[a] " << grad[0] << std::endl;
        // std::cout << "[b] " << grad[1] << std::endl;
        // std::cout << "[c] " << grad[2] << std::endl;

        // std::cout << "Hessian Matrix" << std::endl;
        // std::cout << "[a,a] " << hess[0][0] << std::endl;
        // std::cout << "[a,b] " << hess[0][1] << std::endl;
        // std::cout << "[a,c] " << hess[0][2] << std::endl;
        // std::cout << "[b,b] " << hess[1][1] << std::endl;
        // std::cout << "[b,c] " << hess[2][1] << std::endl;
        // std::cout << "[c,c] " << hess[2][2] << std::endl;
    };

    // In this example we demonstrate exponential fitting
    // by f(x) = a+ b * (1-exp(-c *x)
    // using function value, gradient and Hessian (with respect to a, b, c)
    std::string param_str = "[0, ";
    double mini = dist_geo_eta_weights[0].at(1);
    double maxi = dist_geo_eta_weights[0].at(1);

    for (auto value : dist_geo_eta_weights)
    {
        maxi = max(value.at(1), maxi);
        mini = min(value.at(1), mini);
    }
    double range = maxi - mini;
    param_str += std::to_string(range);
    param_str += ", ";
    param_str += std::to_string(log(2) / dist_geo_eta_weights.size());
    param_str += "]";

    //std::cout << "Param start : " << param_str << std::endl;
    real_1d_array param = param_str.c_str();

    std::string dist_geo_str = "[[";
    std::string eta_str = "[";
    std::string weights_str = "[";

    for (int i = 0; i < dist_geo_eta_weights.size() - 1; ++i)
    {
        dist_geo_str += std::to_string(dist_geo_eta_weights[i].at(0));
        dist_geo_str += "],[";
        eta_str += std::to_string(dist_geo_eta_weights[i].at(1));
        eta_str += ", ";
        weights_str += std::to_string(dist_geo_eta_weights[i].at(2));
        weights_str += ", ";
    }

    dist_geo_str += std::to_string(dist_geo_eta_weights[dist_geo_eta_weights.size() - 1].at(0));
    dist_geo_str += "]]";
    eta_str += std::to_string(dist_geo_eta_weights[dist_geo_eta_weights.size() - 1].at(1));
    eta_str += "]";
    weights_str += std::to_string(dist_geo_eta_weights[dist_geo_eta_weights.size() - 1].at(2));
    weights_str += "]";

    real_2d_array dist_geo = dist_geo_str.c_str();
    real_1d_array eta = eta_str.c_str();
    real_1d_array weights = weights_str.c_str();

    double epsx = 0.000001;
    ae_int_t maxits = 0;
    ae_int_t info;
    lsfitstate state;
    lsfitreport rep;

    // Fitting with weights
    lsfitcreatewfgh(dist_geo, eta, weights, param, state);
    lsfitsetcond(state, epsx, maxits);
    alglib::lsfitfit(state, function_cx_1_func, function_cx_1_grad, function_cx_1_hess);
    lsfitresults(state, info, param, rep);
    // printf("%d\n", int(info));

    //return a, b, b_g
    return std::array<double, 3>{param[0], param[1], param[2]};
}