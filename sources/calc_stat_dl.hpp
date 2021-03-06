#pragma once

#include "data_plane_vec.hpp"

double calc_phi_ij_xy(std::array<std::array<int, 4>, 2> const &locus_chr);
//Phi for all indiv, calc without missing data
std::vector<double> calc_phi_ij(data_plane_vec_c const &data_plane_vec, int ploidy);

//<locus_pair_nbr<value eta>>
std::array<double, 3> calc_eta_ij_per_deme_pair_xy(data_plane_vec_c const &data_plane_vec, int chr, int locus_i, double Q2_loc_i, double Q1_loc_i, int locus_j, double Q2_loc_j, double Q1_loc_j, int deme_x, int deme_y);
//<locus_pair_nbr,<dist-deme, dist-locus, value eta>>
std::vector<std::array<double, 5>> calc_eta_ij_per_deme_pair(data_plane_vec_c const &data_plane_vec, int chr, int locus_i, double Q2_loc_i, int locus_j, double Q2_loc_j);
//<pair of deme * pair of locus,<dist-deme, dist-locus, value eta>>
std::vector<std::array<double, 5>> calc_eta_diploide(data_plane_vec_c const &data_plane_vec);
//Haploid version with > 1 indiv/deme
std::vector<std::array<double, 5>> calc_eta_haploid(data_plane_vec_c const &data_plane_vec);
//Continous habitat isolation by distance
std::vector<std::array<double, 5>> calc_eta_1_indiv_deme_v(data_plane_vec_c const &data_plane_vec);
//Exp regresssion
std::array<double, 3> exp_regr(std::vector<std::array<double, 3>> const &dist_geo_eta_weights);