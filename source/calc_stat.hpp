#pragma once
#include <array>

#include "data_plane_vec.hpp"

//TODO : vecteur de mot clé pour pouvoir comparer un ensemble d'attribut (indiv, pop, habitat) => fonction de comparaison
double calc_Q_intra_indiv(data_plane_vec_c const &data_plane_vec);
double calc_Q_inter_indiv_intra_pop(data_plane_vec_c const &data_plane_vec);
double calc_Q_inter_pop(data_plane_vec_c const &data_plane_vec);

std::vector<std::array<double, 2>> calc_qr_loc_by_loc(data_plane_vec_c const &data_plane_vec, int dist_max, int locus);
std::vector<std::array<double, 2>> calc_qr_all_loc(data_plane_vec_c const &data_plane_vec, int dist_max);
// std::vector<std::array<double, 2>> ar_ln_dist_qr(data_plane_vec_c const &data_plane_vec, int dist_max);
//std::vector<std::array<double, 3>> ar_by_locus_by_pair(data_plane_vec_c const &data_plane_vec, int locus);
std::vector<std::array<double, 2>> ar_by_pair(data_plane_vec_c const &data_plane_vec);
std::vector<std::array<double, 2>> er_by_pair(data_plane_vec_c const &data_plane_vec);
//esti Fis, Fst, Fit
std::array<std::array<double, 2>, 2> Fstat_by_loc_with_probid(data_plane_vec_c const &data_plane_vec, int locus);
std::array<std::array<double, 2>, 2> Fstat_by_loc_with_indic(data_plane_vec_c const &data_plane_vec, int locus);
std::array<double, 2> Fstat_genepop(data_plane_vec_c const &data_plane_vec);