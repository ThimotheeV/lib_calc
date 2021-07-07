#pragma once
#include <array>

#include "data_plane_vec.hpp"

//WARNING : All calc are based in assumption than no empty locus are provided in the data set
//TODO : vecteur de mot clé pour pouvoir comparer un ensemble d'attribut (indiv, deme, habitat) => fonction de comparaison
double calc_Q_intra_indiv_per_chr_per_locus(data_plane_vec_c const &data_plane_vec, int chr, int locus);
double calc_Q_intra_indiv(data_plane_vec_c const &data_plane_vec);

double calc_Q_inter_indiv_per_chr_per_locus_per_deme(data_plane_vec_c const &data_plane_vec, int chr, int locus, int deme);
double calc_Q_inter_indiv_per_chr_per_locus(data_plane_vec_c const &data_plane_vec, int chr, int locus);
double calc_Q_inter_indiv_intra_deme(data_plane_vec_c const &data_plane_vec);

double calc_Q_inter_deme_per_chr_per_locus(data_plane_vec_c const &data_plane_vec, int chr, int locus);
double calc_Q_inter_deme(data_plane_vec_c const &data_plane_vec);

double calc_Hnei_per_chr_per_loc(data_plane_vec_c const &data_plane_vec, int chr, int locus);
double calc_Hnei(data_plane_vec_c const &data_plane_vec);
//WARNING : for microsat only
double calc_Var_per_chr_per_loc(data_plane_vec_c const &data_plane_vec, int chr, int locus);
//calc_Q_intra_indiv => calc_Qwi_frac
double calc_Hobs_per_chr_per_loc(data_plane_vec_c const &data_plane_vec, int chr, int locus);
double calc_Hobs(data_plane_vec_c const &data_plane_vec);
double calc_MGW_per_chr_per_loc(data_plane_vec_c const &data_plane_vec, int chr, int locus);
double calc_MGW(data_plane_vec_c const &data_plane_vec);

// std::vector<qr_num, qr_denum>
std::vector<double> calc_qr_per_chr_by_loc(data_plane_vec_c const &data_plane_vec, int chr, int locus);
// std::vector<qr>
std::vector<double> calc_qr(data_plane_vec_c const &data_plane_vec);

std::vector<std::array<double, 2>> ar_by_pair(data_plane_vec_c const &data_plane_vec);
std::vector<std::array<double, 2>> er_by_pair(data_plane_vec_c const &data_plane_vec);
//esti Fis, Fst, Fit
//WARNING : Not valid if missing value or indiv_per_deme < 2 but more efficient for big data
std::array<std::array<double, 2>, 2> Fstat_per_chr_by_loc_with_probid(data_plane_vec_c const &data_plane_vec, int chr, int locus);
//WARNING : Validity if missing value but less efficient for big data. Not valid if indiv_per_deme < 2
std::array<std::array<double, 2>, 2> Fstat_per_chr_by_loc_with_indic(data_plane_vec_c const &data_plane_vec, int chr, int locus);
std::array<double, 2> Fstat_genepop(data_plane_vec_c const &data_plane_vec, bool missing_data);

//Return <state, <frequence, nbr of locus>>
std::map<int, std::map<int, double>> calc_AFS(data_plane_vec_c const &data_plane_vec, int limit_min_gene_per_locus);